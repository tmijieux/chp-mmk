module grad
  implicit none

contains
  subroutine grad_conj(Nx,Ny,B,Cx,Cy,X,rhs)
    implicit none
    integer,intent(in)::Nx, Ny
    real*8,intent(inout),dimension(Nx*Ny)::X
    real*8,intent(in),dimension(Nx*Ny)::rhs
    real*8,intent(in)::B,Cx,Cy
    real*8,dimension(Nx*Ny)::w,d,R,Rnext,Ax
    integer::k
    integer,parameter::niter=20000000
    real*8,parameter::eps=1e-5
    real*8::beta,alpha

    call prodAx(B,Cx,Cy,X,Nx,Ny,Ax)
    R=Ax-rhs
    d=R
    k=0
    do while(norme(R)>eps .and. k<=niter)
       call prodAx(B,Cx,Cy,d,Nx,Ny,w)
       alpha=dot_product(d,R)/dot_product(d,w)
       X=X-alpha*d
       Rnext=R-alpha*w
       beta=norme(Rnext)**2/norme(R)**2
       d=Rnext+beta*d
       R=Rnext
       k=k+1
    end do
  end subroutine grad_conj

  subroutine CalculRHS(Cx,Cy,h,g,RHS)
    implicit none
    real*8, intent(in) :: Cx,Cy
    real*8, dimension(:),intent(in) :: h,g
    real*8, dimension(:),intent(inout) :: RHS
    integer :: Nx,Ny,i

    Nx=size(g)/2
    Ny=size(h)/2
    write(*,*) Nx, Ny, size(RHS)

    do i=1, size(RHS)
       if (i<=Nx) then
          RHS(i) = RHS(i) - g(i)*Cy
       elseif (i>(Nx*Ny-Nx)) then
          RHS(i) = RHS(i) - g(i-(Nx*Ny-2*Nx))*Cy
       end if

       if (mod(i,Nx) == 1) then
          RHS(i) = RHS(i) - h(1+(i/Nx))*Cx
       elseif (mod(i, Nx) == 0) then
          RHS(i) = RHS(i) - h(Ny +(i/Nx))*Cx
       end if
    enddo

  end subroutine CalculRHS

  subroutine prodAx(B,Cx,Cy,X,Nx,Ny,AX)
    implicit none
    integer,intent(in)::Nx,Ny
    real*8,intent(in)::B,Cx,Cy
    real*8,dimension(Nx*Ny),intent(in)::X
    real*8,dimension(Nx*Ny),intent(out)::AX
    integer::i

    AX(1)=B*X(1)+Cx*X(2)+Cy*X(Nx+1)
    do i=2,Nx-1
       AX(i)=B*X(i)+Cx*(X(i-1)+X(i+1))+Cy*X(Nx+i)
    end do
    AX(Nx)=B*X(Nx)+Cx*X(Nx-1)+Cy*X(2*Nx)

    do i=Nx+1,Nx*Ny-Nx
       if(mod(i,Nx)==1)then
          AX(i)=Cy*(X(i-Nx)+X(i+Nx))+Cx*X(i+1)+B*X(i)
       else if(mod(i,Nx)==0)then
          AX(i)=Cy*(X(i-Nx)+X(i+Nx))+Cx*X(i-1)+B*X(i)
       else
          AX(i)=Cy*(X(i-Nx)+X(i+Nx))+Cx*(X(i+1)+X(i-1))+B*X(i)
       end if
    end do

    AX(Nx*Ny-Nx+1)=B*X(Nx*Ny-Nx+1)+Cx*X(Nx*Ny-Nx+2)+Cy*X(Nx*Ny-2*Nx+1)
    do i=Nx*Ny-Nx+2,Nx*Ny-1
       AX(i)=B*X(i)+Cx*(X(i-1)+X(i+1))+Cy*X(i-Nx)
    end do
    AX(Nx*Ny)=B*X(Nx*Ny)+Cx*X(Nx*Ny-1)+Cy*X(Nx*Ny-Nx)

  end subroutine prodAx

  real*8 function norme(X)
    implicit none
    real*8,dimension(:)::X
    integer::i
    norme=0.
    do i=1,size(X)
       norme=norme+X(i)*X(i)
    end do
    norme=sqrt(norme)
  end function norme

end module grad
