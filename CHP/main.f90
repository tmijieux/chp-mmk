program chp
  use fonctions
  use grad
  implicit none

  character(len=20)::parametres='param.txt'
  integer::choixT
  real*8::start,finish
  
  write(*,*) "1.Stationnaire 2.Instationnaire"
  read(*,*) choixT
  if(choixT==1)then
     call stationnaire(parametres)
  else if(choixT==2)then
     call instationnaire(parametres)
  else
     write(*,*) "erreur"
  end if

  call cpu_time(finish)
  write(*,*) finish-start,"seconds"

contains

  subroutine stationnaire(parametres)
    implicit none
    character(len=20),intent(in)::parametres
    integer::Nx,Ny,choixF,i,j
    real*8::Lx,Ly,D,dx,dy,B,Cx,Cy
    real*8,dimension(:),allocatable::X,Y,U,Uexacte,RHS,f,g,h
    
    call cpu_time(start)

    open(unit=10,file=parametres)
    read(10,*)Nx,Ny,Lx,Ly,D
    close(10)

    allocate(U(Nx*Ny))
    allocate(Uexacte(Nx*Ny))
    allocate(RHS(Nx*Ny))
    allocate(f(Nx*Ny))
    allocate(g(2*Nx))
    allocate(h(2*Ny))
    allocate(X(Nx))
    allocate(Y(Ny))

    D=1.
    dx=Lx/(1+Nx)
    dy=Ly/(1+Ny)

    call Init(X,Y,dx,dy,Nx,Ny)

    B=2*D/dx**2+2*D/dy**2 
    Cx=-D/dx**2
    Cy=-D/dy**2
    U=0.

    write(*,*) "fonctions 1 ou fonctions 2 ?"
    read(*,*) choixF
    if(choixF==1)then
       f=f1(X,Y)
       g=g1(X)
       h=h1(Y)
       do i=1,Nx
          do j=1,Ny
             Uexacte((j-1)*Nx+i)=X(i)*(Lx-X(i))*Y(j)*(Ly-Y(j))
          end do
       end do
    else if(choixF==2)then
       f=f2(X,Y)
       g=g2(X,Ly)
       h=h2(Y,Lx)
       do i=1,Nx
          do j=1,Ny
             Uexacte((j-1)*Nx+i)=sin(X(i))+cos(Y(j))
          end do
       end do
    else
       write(*,*)"erreur"
    end if

    RHS=U+f
    call CalculRHS(Cx,Cy,h,g,RHS)
    call grad_conj(Nx,Ny,B,Cx,Cy,U,RHS)

    open(unit=11,file='numerique.dat')
    call ecriture(11,X,Y,U)
    close(11)
    print*,U
    open(unit=12,file='exacte.dat')
    call ecriture(12,X,Y,Uexacte)
    close(12)

    deallocate(U)
    deallocate(Uexacte)
    deallocate(RHS)
    deallocate(f)
    deallocate(g)
    deallocate(h)
    deallocate(X)
    deallocate(Y)

  end subroutine stationnaire

  subroutine instationnaire(parametres)
    implicit none
    character(len=20),intent(in)::parametres
    integer::Nx,Ny,Nit,i,step
    real*8::Lx,Ly,D,dt,dx,dy,B,Cx,Cy,t,Tmax
    real*8,dimension(:),allocatable::U0,U,RHS,f,g,h,X,Y,U4,U8
    character(len=20)::filename

    open(unit=10,file=parametres)
    read(10,*)Nx,Ny,Lx,Ly,D
    close(10)


    allocate(U0(Nx*Ny))
    allocate(U(Nx*Ny))
    allocate(RHS(Nx*Ny))
    allocate(f(Nx*Ny))
    allocate(g(2*Nx))
    allocate(h(2*Ny))
    allocate(X(Nx))
    allocate(Y(Ny))
    allocate(U4(Nx*Ny))
    allocate(U8(Nx*Ny))

    Tmax=10
    Nit=2000
    D=1.
    dx=Lx/(1+Nx)
    dy=Ly/(1+Ny)
    dt=Tmax/Nit

    call Init(X,Y,dx,dy,Nx,Ny)
    
    !B=1+2*D*dt/dx**2+2*D*dt/dy**2
    !Cx=-D*dt/dx**2
    !Cy=-D*dt/dy**2
    B=1/dt+2*D/dx**2+2*D/dy**2
    Cx=-D/dx**2
    Cy=-D/dy**2

    U0=0.
    t=0.
    do i=1,Nit
       t=t+dt
       f=f3(X,Y,Lx,Ly,t)
       !RHS=U0+dt*f
       RHS=U0/dt+f
       call CalculRHS(Cx,Cy,h,g,RHS)
       call grad_conj(Nx,Ny,B,Cx,Cy,U,RHS)
       U0=U

       !gnuplot : do for [i=1:2000]{ splot "sol/sol".i.".dat" w pm3d ; pause 0.05}
       step=1
       if(mod(i,step)==0)then
          if(i/step<10)then
             write(filename,"(A7,I1,A4)")"sol/sol",i/step,".dat"
          else if(i/step<100)then
             write(filename,"(A7,I2,A4)")"sol/sol",i/step,".dat"
          else if(i/step<1000)then
             write(filename,"(A7,I3,A4)")"sol/sol",i/step,".dat"
          else
             write(filename,"(A7,I4,A4)")"sol/sol",i/step,".dat"
          end if
          open(unit=i+11,file=filename)
          call ecriture(i+11,X,Y,U)
          close(i+11)
       end if

    end do
    write(*,*)"nombre de fichiers",Nit/step
    print*,U
    open(unit=11,file='numeriqueInsta.dat')
    call ecriture(11,X,Y,U)
    close(11)

    deallocate(U0)
    deallocate(U)
    deallocate(RHS)
    deallocate(f)
    deallocate(g)
    deallocate(h)
    deallocate(X)
    deallocate(Y)
    deallocate(U4)
    deallocate(U8)

  end subroutine instationnaire

end program chp
