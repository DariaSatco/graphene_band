program graphene_tb_prog
use energy
use matrix_tb
use matrix_etb

    implicit none

real, allocatable :: kx(:), ky(:)
real, allocatable :: energy_band_tb(:,:,:), energy_band_etb(:,:,:), energy_band_etb_plot(:,:)
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

open(unit=17, file='graphene_bands_etb_plot.dat', status='replace')
write(17,150) 'n', 'kx', 'ky'


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

!ETB for plotting in linear space Gamma -> K -> M -> Gamma ********************
allocate(energy_band_etb_plot(n,8))

! Gamma -> K ======================

do i = 1,n

    kx(i) = 0. + (i-1)*2*pi/(3*a*n)
    ky(i) = 0. + kx(i)*2/sqrt(3.)

    matr_H_etb = matrix_H_etb(kx(i),ky(i))
    matr_S_etb = matrix_S_etb(kx(i),ky(i))

    call chegv(1, 'N', 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_etb_plot(i,1:8), work_etb, 100 , rwork_etb, info)

    write(17,100) real(i), kx(i)*a/pi, ky(i)*a/pi, energy_band_etb_plot(i,1:8)
end do

! K -> M ===========================

do i = 1,n

    kx(i) = 2*pi/(3*a)
    ky(i) = 4*pi/(3*sqrt(3.)*a) - (i-1)*4*pi/(3*sqrt(3.)*a*n)

    matr_H_etb = matrix_H_etb(kx(i),ky(i))
    matr_S_etb = matrix_S_etb(kx(i),ky(i))

    call chegv(1, 'N', 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_etb_plot(i,1:8), work_etb, 100 , rwork_etb, info)

    write(17,100) real(i+n), kx(i)*a/pi, ky(i)*a/pi, energy_band_etb_plot(i,1:8)
end do

! M -> Gamma ==========================
do i = 1,n

    kx(i) = 2*pi/(3*a) - (i-1)*2*pi/(3*a*n)
    ky(i) = 0.

    matr_H_etb = matrix_H_etb(kx(i),ky(i))
    matr_S_etb = matrix_S_etb(kx(i),ky(i))

    call chegv(1, 'N', 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_etb_plot(i,1:8), work_etb, 100 , rwork_etb, info)

    write(17,100) real(i+2*n), kx(i)*a/pi, ky(i)*a/pi, energy_band_etb_plot(i,1:8)
end do


deallocate(kx,ky)
deallocate(energy_band_tb)
deallocate(energy_band_etb)

close(15)
close(16)
close(17)

end program graphene_tb_prog
