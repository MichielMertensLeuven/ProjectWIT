program main
  implicit none
  integer(kind = 4) nx,ny,i
  real(kind = 8) :: startmean, diffmean, rhodiff, rhodefault
  real(kind = 8), dimension(nx*ny) :: rho, grad
  
  rhodefault = 0.4D+00
  nx = 20
  ny = 20
  rho = rhodefault
  rhodiff = 0.01D+00
  call heat(startmean, rho, nx, ny)
    
  do i = 1,nx*ny
    rho(i) += rhodiff
    call heat(diffmean, rho, nx, ny)
    grad(i) = (startmean-diffmean)/rhodiff
    rho(i) = rhodefault
  enddo  
  
  
end program main  
