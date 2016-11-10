module fonctions
  implicit none

contains

  !g1
  function g1(X)
    implicit none
    real*8,dimension(:)::X
    real*8,dimension(2*size(X))::g1
    g1=0.
  end function g1

  !h1
  function h1(Y)
    implicit none
    real*8,dimension(:)::Y
    real*8,dimension(2*size(Y))::h1
    h1=0.
  end function h1

  function f1(X,Y)
    implicit none
    real*8,dimension(:)::X,Y
    real*8,dimension(size(X)*size(Y))::f1
    integer::i,j,Nx,Ny

    Nx=size(X)
    Ny=size(Y)

    do i=1,Nx
       do j=1,Ny
          f1((j-1)*Nx+i)=2*(Y(j)-Y(j)**2+X(i)-X(i)**2)
       end do
    end do

  end function f1

  !g2
  function g2(X,Ly)
    implicit none
    real*8::Ly
    real*8,dimension(:)::X
    real*8,dimension(2*size(X))::g2
    integer::i

    do i=1, size(X)
       g2(i) = sin(X(i)) + cos(0.)
    end do

    do i=1, size(X)
       g2(i+size(X)) = sin(X(i)) + cos(Ly)
    end do
  end function g2

  !h2
  function h2(Y,Lx)
    implicit none
    real*8::Lx
    real*8,dimension(:)::Y
    real*8,dimension(2*size(Y))::h2
    integer::i

    do i=1,size(Y)
       h2(i) = sin(0.) + cos(Y(i))
    end do

    do i=1, size(Y)
       h2(i+size(Y)) = sin(Lx) + cos(Y(i))
    end do

  end function h2

  !f2
  function f2(X,Y)
    implicit none
    real*8,dimension(:)::X,Y
    real*8,dimension(size(X)*size(Y))::f2
    integer::i,j,Nx,Ny

    Nx=size(X)
    Ny=size(Y)

    do i=1,Ny
       do j=1,Nx
          f2((i-1)*Nx+j) = sin(X(j))+cos(Y(i))
       end do
    end do

  end function f2


  !g3
  function g3(X)
    implicit none
    real*8,dimension(:)::X
    real*8,dimension(2*size(X))::g3
    g3=0.
  end function g3

  !h3
  function h3(Y)
    implicit none
    real*8,dimension(:)::Y
    real*8,dimension(2*size(Y))::h3
    h3=1.0

  end function h3

  !f3
  function f3(X,Y,Lx,Ly,t)
    implicit none
    real*8::t,Lx,Ly,pi
    real*8,dimension(:)::X,Y
    real*8,dimension(size(X)*size(Y))::f3
    integer::i,j,Nx,Ny
    pi=4*atan(1.)
    Nx=size(X)
    Ny=size(Y)

    do i=1,Nx
       do j=1,Ny
          f3((j-1)*Nx+i)=exp(-(X(i)-0.5*Lx)**2)*exp(-(Y(i)-0.5*Ly)**2)*cos(pi*t/2)
       end do
    end do

  end function f3

  subroutine Init(X,Y,dx,dy,Nx,Ny)
    implicit none
    integer,intent(in)::Nx,Ny
    real*8,intent(in)::dx,dy
    integer::i
    real*8,dimension(:),intent(out)::X
    real*8,dimension(:),intent(out)::Y

    do i=1,Nx
       X(i)=i*dx
    end do
    do i=1,Ny
       Y(i)=i*dy
    end do

  end subroutine Init

  subroutine ecriture(fichier,X,Y,U)
    implicit none
    integer,intent(in)::fichier
    real*8,dimension(:),intent(in)::X,Y,U
    integer::Nx,Ny,i,j

    Nx=size(X)
    Ny=size(Y)

    do i=1,Nx
       do j=1,Ny
          write(fichier,*) X(i),Y(j),U((j-1)*Nx+i)
       end do
       write(fichier,*) " "
    end do
  end subroutine ecriture

end module fonctions
