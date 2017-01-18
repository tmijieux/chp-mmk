program test_grad
  use grad
  implicit none
  integer,parameter::Nx=10,Ny=10
  real*8,dimension(Nx*Ny)  :: A
  real*8,dimension(2*Nx)   :: g
  real*8,dimension(2*Ny)   :: h
  real*8                   :: Cx,Cy
  integer                  :: i

  Cx=0.1
  Cy=0.1

  do i = 1, Nx*Ny
     A(i) = 0.0
  end do

  do i=1,2*Nx
     g(i) = 3.0
  end do
  do i=1,2*Ny
     h(i) = 1.0
  end do

  call CalculRHS(Cx, Cy, h, g, A)
  print*,A

end program test_grad
