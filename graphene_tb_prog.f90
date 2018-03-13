program graphene_tb_prog
use energy
    implicit none

real, allocatable :: kx(:), ky(:)
real, allocatable :: energy_band(:,:,:)
integer :: i,j,n
real :: dk

100 format(8(f15.3))
150 format(8(A15))

open(unit=15, file='graphene_bands.dat', status='replace')
write(15,150) 'kx', 'ky', 'energy_plus', 'energy_minus'

n=100
allocate(kx(n), ky(n))
allocate(energy_band(n,n,2))
dk=2*pi/((n-1)*a)

do i = 1, n

    kx(i) = -pi/a + dk*(i-1)

    do j = 1, n
        ky(j) = -pi/a + dk*(j-1)
        energy_band(i,j,1:2) = energy_tb(kx(i), ky(j))
    write(15,100) kx(i)*a/pi, ky(j)*a/pi, energy_band(i,j,1), energy_band(i,j,2)
    end do
write(15,*) ' '
end do

deallocate(kx,ky)
deallocate(energy_band)

close(15)

end program graphene_tb_prog
