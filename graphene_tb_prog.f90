program graphene_tb_prog
use energy
use matrix_tb
use matrix_etb

    implicit none

real, allocatable :: kx(:), ky(:)
real, allocatable :: energy_band_tb(:,:,:), energy_band_etb(:,:,:), energy_band_cnt(:,:,:)
integer :: i,j,n, nc, nt, kc, kt
real :: dk, kfin
real(8):: kcx,kcy,ktx,kty
real:: kxn, kyn
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

! create output files
open(unit=15, file='graphene_bands_tb.dat', status='replace')
write(15,150) 'kx', 'ky', 'energy_plus', 'energy_minus'

open(unit=16, file='graphene_bands_etb.dat', status='replace')
write(16,150) 'kx', 'ky'

open(unit=17, file='graphene_bands_etb_plot.dat', status='replace')
write(17,150) 'k', 'kx', 'ky'

open(unit=18, file='graphene_bands_tb_plot.dat', status='replace')
write(18,150) 'k', 'kx', 'ky'

! n -number of k-points
n=100
allocate(kx(n), ky(n))
allocate(energy_band_tb(n,n,2))
allocate(energy_band_etb(n,n,8))

! step for k-vector change
dk=2*pi/((n-1)*a)

! TB energy dispersion for graphene ******************************************
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

! ETB energy dispersion for graphene ****************************************
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

!ETB and TB for plotting in linear space of characteristic path Gamma -> K -> M -> Gamma ********************

! Gamma -> K =================================================

do i = 1,n

    kx(i) = 0. + (i-1)*2*pi/(3*a*(n-1))
    ky(i) = 0. + kx(i)/sqrt(3.)

! TB ===================================

    matr_H_tb = matrix_H(kx(i),ky(i))
    matr_S_tb = matrix_S(kx(i),ky(i))

    call chegv(1, 'N', 'U', 2, matr_H_tb, 2, matr_S_tb, 2, energy_band_tb(i,i,1:2), work, 20 , rwork, info)

    write(18,100) sqrt(kx(i)**2+ky(i)**2)*a/pi, kx(i)*a/pi, ky(i)*a/pi, energy_band_tb(i,i,1:2)

! ETB ==================================

    matr_H_etb = matrix_H_etb(kx(i),ky(i))
    matr_S_etb = matrix_S_etb(kx(i),ky(i))

    call chegv(1, 'N', 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_etb(i,i,1:8), work_etb, 100 , rwork_etb, info)

    write(17,100) sqrt(kx(i)**2+ky(i)**2)*a/pi, kx(i)*a/pi, ky(i)*a/pi, energy_band_etb(i,i,1:8)
end do

! we calculate the final length of k-vec to continue from its end point in the next loop
kfin = sqrt(kx(n)**2+ky(n)**2)*a/pi

! K -> M ===================================================

do i = 1,n

    kx(i) = 2*pi/(3*a)
    ky(i) = 2*pi/(3*sqrt(3.)*a) - i*2*pi/(3*sqrt(3.)*a*n)

! TB =================================

    matr_H_tb = matrix_H(kx(i),ky(i))
    matr_S_tb = matrix_S(kx(i),ky(i))

    call chegv(1, 'N', 'U', 2, matr_H_tb, 2, matr_S_tb, 2, energy_band_tb(i,i,1:2), work, 20 , rwork, info)

    write(18,100) kfin + (2*pi/(3*sqrt(3.)*a)-ky(i))*a/pi, kx(i)*a/pi, ky(i)*a/pi, energy_band_tb(i,i,1:2)

! ETB ==================================

    matr_H_etb = matrix_H_etb(kx(i),ky(i))
    matr_S_etb = matrix_S_etb(kx(i),ky(i))

    call chegv(1, 'N', 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_etb(i,i,1:8), work_etb, 100 , rwork_etb, info)

    write(17,100) kfin + (2*pi/(3*sqrt(3.)*a)-ky(i))*a/pi, kx(i)*a/pi, ky(i)*a/pi, energy_band_etb(i,i,1:8)
end do

kfin = kfin + 2*pi/(3*sqrt(3.)*a)*a/pi

! M -> Gamma =====================================================
do i = 1,n

    kx(i) = 2*pi/(3*a) - i*2*pi/(3*a*n)
    ky(i) = 0.

! TB ================================

    matr_H_tb = matrix_H(kx(i),ky(i))
    matr_S_tb = matrix_S(kx(i),ky(i))

    call chegv(1, 'N', 'U', 2, matr_H_tb, 2, matr_S_tb, 2, energy_band_tb(i,i,1:2), work, 20 , rwork, info)

    write(18,100) kfin + (2*pi/(3*a)-kx(i))*a/pi, kx(i)*a/pi, ky(i)*a/pi, energy_band_tb(i,i,1:2)

! ETB ===================================

    matr_H_etb = matrix_H_etb(kx(i),ky(i))
    matr_S_etb = matrix_S_etb(kx(i),ky(i))

    call chegv(1, 'N', 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_etb(i,i,1:8), work_etb, 100 , rwork_etb, info)

    write(17,100) kfin + (2*pi/(3*a)-kx(i))*a/pi, kx(i)*a/pi, ky(i)*a/pi, energy_band_etb(i,i,1:8)
end do

deallocate(kx,ky)
deallocate(energy_band_tb)
deallocate(energy_band_etb)

close(15)
close(16)
close(17)
close(18)

! ETB for carbon nanotibes =====================================
call tube_geom(kcx,kcy,ktx,kty,nc,nt)
allocate(energy_band_cnt(nc,nt,8))

open(unit=19, file='cnt_bands_etb.dat', status='replace')
write(19,150) 'k', 'kx', 'ky', 'energy'

do kc = 1 - nc/2, nc/2
    do kt = 1 - nt/2, nt/2
    kxn = kc*kcx + kt*ktx
    kyn = kc*kcy + kt*kty

    matr_H_etb = matrix_H_etb(kxn,kyn)
    matr_S_etb = matrix_S_etb(kxn,kyn)

    call chegv(1, 'N', 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_cnt(kc,kt,1:8), work_etb, 100 , rwork_etb, info)

    write(19,100) kt*sqrt(ktx**2+kty**2)*sqrt(3.)*a, kxn*a, kyn*a, energy_band_cnt(kc,kt,1:8)
    end do
write(19,*) ' '
end do

close(19)
deallocate(energy_band_cnt)

end program graphene_tb_prog

!-------------------------------------------------------------------------------

    include '/Users/dariasatco/Documents/eclipse_projects/georgiifortran/tbdftsnt/entube.f'
    include '/Users/dariasatco/Documents/eclipse_projects/georgiifortran/tbdftsnt/tbtube.f'
    include '/Users/dariasatco/Documents/eclipse_projects/georgiifortran/tbdftsnt/tbpot.f'
    include '/Users/dariasatco/Documents/eclipse_projects/georgiifortran/cntgrlib/gcd.f'
    include '/Users/dariasatco/Documents/eclipse_projects/georgiifortran/cntgrlib/tubpar.f'

!-------------------------------------------------------------------------------

!subroutine to generate geometry of cnt reciprocal space
    subroutine tube_geom(kcx,kcy,ktx,kty,nc,nt)

    implicit none

    real*8 l
    parameter (l=+0.3d+02)
    real*8 pi,sr3
    parameter (pi=+0.3141592653589793d+01,sr3=+0.1732050807568877d+01)!sr3 means sqrt(3)
    real*8 acc,auc
    parameter (acc=+0.142d-9,auc=acc*sr3) !acc = interatomic distance, auc = sqrt(3)*acc
    integer:: n,m,nc,nt
    real*8 t,d,theta,tab,phiab,dab,a1t,a1phi,a2t,a2phi,kcx,kcy,ktx,kty

! tube indices (n,m)
    print*, 'Enter tube indices (n,m):'
    print*, 'n:'
    read*, n
    print*, 'm:'
    read*, m
    !n=7
    !m=7
    call tubpar (n,m,l,t,d,theta,tab,phiab,dab,a1t,a1phi,a2t,a2phi,nc,nt,kcx,kcy,ktx,kty)

    kcx=kcx/auc
    kcy=kcy/auc
    ktx=ktx/auc
    kty=kty/auc

    end subroutine tube_geom

