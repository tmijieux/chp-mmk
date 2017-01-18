program test
  use grad
  implicit none
  integer,parameter ::  Nx=10,Ny=10
  real*8,dimension(Nx*Ny) :: X, AX
  real*8  :: B,Cx,Cy,dx,dy,Lx=1.,Ly=1.
  integer :: i

  dx = Lx/(1+Nx)
  dy = Ly/(1+Ny)
  B = 3.0
  Cx = 1.0
  Cy = 1.0

  do i = 1,size(X)
     X(i) = 1.0 * i
  end do

  call prodAx(B,Cx,Cy,X,Nx,Ny,AX)
  print*,AX

end program test
