program graphene_tb_prog
use energy
use matrix_tb
use matrix_etb

    implicit none

real, allocatable :: kx(:), ky(:)
real, allocatable :: energy_band_tb(:,:,:), energy_band_etb(:,:,:)
integer :: i,j,n
real :: dk
complex :: matr_H_tb(2,2), matr_S_tb(2,2)
complex :: matr_H_etb(8,8), matr_S_etb(8,8)

!for chegv TB (N=2)
integer :: info !INFO is INTEGER
complex :: work(20)
! The length of the array WORK.  LWORK >= max(1,2*N-1).
!          For optimal efficiency, LWORK >= (NB+1)*N
real :: rwork(4)
! RWORK is REAL array, dimension (max(1, 3*N-2))

!for chegv ETB (N=8)
complex :: work_etb(100)
real :: rwork_etb(22)


100 format(20(f15.3))
150 format(20(A15))

open(unit=15, file='graphene_bands_tb.dat', status='replace')
write(15,150) 'kx', 'ky', 'energy_plus', 'energy_minus'

open(unit=16, file='graphene_bands_etb.dat', status='replace')
write(16,150) 'kx', 'ky'

n=100
allocate(kx(n), ky(n))
allocate(energy_band_tb(n,n,2))
allocate(energy_band_etb(n,n,8))

dk=2*pi/((n-1)*a)

! TB ******************************************
do i = 1, n

    kx(i) = -pi/a + dk*(i-1)

    do j = 1, n
        ky(j) = -pi/a + dk*(j-1)
        matr_H_tb = matrix_H(kx(i),ky(j))
        matr_S_tb = matrix_S(kx(i),ky(j))

        call chegv(1, 'N', 'U', 2, matr_H_tb, 2, matr_S_tb, 2, energy_band_tb(i,j,1:2), work, 20 , rwork, info)

        !energy_band(i,j,1:2) = energy_tb(kx(i), ky(j))
    write(15,100) kx(i)*a/pi, ky(j)*a/pi, energy_band_tb(i,j,1), energy_band_tb(i,j,2)
    end do

write(15,*) ' '
end do

! ETB ****************************************
do i = 1, n

    kx(i) = -pi/a + dk*(i-1)

    do j = 1, n
        ky(j) = -pi/a + dk*(j-1)
        matr_H_etb = matrix_H_etb(kx(i),ky(j))
        matr_S_etb = matrix_S_etb(kx(i),ky(j))

        call chegv(1, 'N', 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_etb(i,j,1:8), work_etb, 100 , rwork_etb, info)

        !energy_band(i,j,1:2) = energy_tb(kx(i), ky(j))
    write(16,100) kx(i)*a/pi, ky(j)*a/pi, energy_band_etb(i,j,1:8)
    end do

write(16,*) ' '
end do

deallocate(kx,ky)
deallocate(energy_band_tb)
deallocate(energy_band_etb)

close(15)
close(16)

end program graphene_tb_prog
