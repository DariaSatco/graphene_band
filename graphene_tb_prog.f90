program graphene_tb_prog

use energy
! energy.f90 includes TB analytical expression for graphene energy dispersion
use matrix_tb
! matrix_tb.f90 includes hamiltonian and overlap 2x2 matrices, simple TB model
use matrix_etb
! matrix_etb.f90 includes hamiltonian and overlap 8x8 matrices, 4-orbital TB
use vec_mult
use ordering

    implicit none

real, allocatable :: kx(:), ky(:)
! kx, ky - reciprocal vectors
real, allocatable :: energy_band_tb(:,:,:), energy_band_etb(:,:,:), energy_band_cnt(:,:,:)
! energy_band_tb - array with energies within TB model
! energy_band_et - array with energies within 4-orbital TB model
! energy_band_cnt - array with CNT energies, zone-folding approach

complex, allocatable :: wf_cnt(:,:,:,:)
! wf_cnt - array of wave function coefficients
complex, allocatable :: overlap(:,:,:,:)

integer :: i,j,l,m !counters
integer:: np !number of kx, ky points
integer:: nc, nt, kc, kt !see tubpar.f
real :: dk, kfin
! dk - step in reciprocal space
! kfin - additional variable to build path Gamma - K - M - Gamma
real(8):: kcx,kcy,ktx,kty,ttube !see tubpar.f, ttube <-> t

integer, parameter:: d_tb = 2
complex :: matr_H_tb(d_tb,d_tb), matr_S_tb(d_tb,d_tb)

integer, parameter:: d_etb = 8
complex :: matr_H_etb(d_etb,d_etb), matr_S_etb(d_etb,d_etb)

!variables for chegv (see LAPACK)
character*1:: jobz
integer :: info
complex, allocatable :: work(:), work_etb(:)
integer:: lwork, lwork_etb
real, allocatable:: rwork(:), rwork_etb(:)

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

open(unit=25, file='cnt_bands_etb.dat', status='replace')
write(25,150) 'k', 'kx', 'ky', 'energy'

open(unit=26, file='cnt_wf_etb.dat', status='replace')
write(26,150) 'k', 'kx', 'ky', 'wf1', 'wf2', 'wf3', 'wf4', 'wf5', 'wf6', 'wf7', 'wf8'

! n -number of k-points
np=100
allocate(kx(np), ky(np))
allocate(energy_band_tb(np,np,d_tb))
allocate(energy_band_etb(np,np,d_etb))

! step for k-vector change
dk=2*pi/((np-1)*a)

! chegv parameters
! jobz controls the computation of eigenvectors
! = 'N':  Compute eigenvalues only;
! = 'V':  Compute eigenvalues and eigenvectors.

jobz = 'N'
lwork = (d_tb+1)*d_tb + 10
allocate(rwork(max(1,3*d_tb-2)))
allocate(work(max(1,lwork)))

lwork_etb = (d_etb+1)*d_etb + 10
allocate(rwork_etb(max(1,3*d_etb-2)))
allocate(work_etb(max(1,lwork_etb)))

! TB energy dispersion for graphene ******************************************
do i = 1, np

    kx(i) = -pi/a + dk*(i-1)

    do j = 1, np
        ky(j) = -pi/a + dk*(j-1)
        matr_H_tb = matrix_H(kx(i),ky(j))
        matr_S_tb = matrix_S(kx(i),ky(j))

        call chegv(1, jobz, 'U', 2, matr_H_tb, 2, matr_S_tb, 2, energy_band_tb(i,j,1:2), work, lwork , rwork, info)

    write(15,100) kx(i)*a/pi, ky(j)*a/pi, energy_band_tb(i,j,1:2)
    end do

write(15,*) ' '
end do

! ETB energy dispersion for graphene ****************************************
do i = 1, np

    kx(i) = -pi/a + dk*(i-1)

    do j = 1, np
        ky(j) = -pi/a + dk*(j-1)
        matr_H_etb = matrix_H_etb(kx(i),ky(j))
        matr_S_etb = matrix_S_etb(kx(i),ky(j))

        call chegv(1, jobz, 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_etb(i,j,1:8), work_etb, lwork_etb , rwork_etb, info)

    write(16,100) kx(i)*a/pi, ky(j)*a/pi, energy_band_etb(i,j,1:8)
    end do

write(16,*) ' '
end do

!ETB and TB for plotting in linear space of characteristic path Gamma -> K -> M -> Gamma ********************

! Gamma -> K =================================================

do i = 1,np

    kx(i) = 0. + (i-1)*2*pi/(3*a*(np-1))
    ky(i) = 0. + kx(i)/sqrt(3.)

! TB ===================================

    matr_H_tb = matrix_H(kx(i),ky(i))
    matr_S_tb = matrix_S(kx(i),ky(i))

    call chegv(1, jobz, 'U', 2, matr_H_tb, 2, matr_S_tb, 2, energy_band_tb(i,i,1:2), work, lwork, rwork, info)

    write(18,100) sqrt(kx(i)**2+ky(i)**2)*a/pi, kx(i)*a/pi, ky(i)*a/pi, energy_band_tb(i,i,1:2)

! ETB ==================================

    matr_H_etb = matrix_H_etb(kx(i),ky(i))
    matr_S_etb = matrix_S_etb(kx(i),ky(i))

    call chegv(1, jobz, 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_etb(i,i,1:8), work_etb, lwork_etb, rwork_etb, info)

    write(17,100) sqrt(kx(i)**2+ky(i)**2)*a/pi, kx(i)*a/pi, ky(i)*a/pi, energy_band_etb(i,i,1:8)
end do

! we calculate the final length of k-vec to continue from its end point in the next loop
kfin = sqrt(kx(np)**2+ky(np)**2)*a/pi

! K -> M ===================================================

do i = 1,np

    kx(i) = 2*pi/(3*a)
    ky(i) = 2*pi/(3*sqrt(3.)*a) - i*2*pi/(3*sqrt(3.)*a*np)

! TB =================================

    matr_H_tb = matrix_H(kx(i),ky(i))
    matr_S_tb = matrix_S(kx(i),ky(i))

    call chegv(1, jobz, 'U', 2, matr_H_tb, 2, matr_S_tb, 2, energy_band_tb(i,i,1:2), work, lwork, rwork, info)

    write(18,100) kfin + (2*pi/(3*sqrt(3.)*a)-ky(i))*a/pi, kx(i)*a/pi, ky(i)*a/pi, energy_band_tb(i,i,1:2)

! ETB ==================================

    matr_H_etb = matrix_H_etb(kx(i),ky(i))
    matr_S_etb = matrix_S_etb(kx(i),ky(i))

    call chegv(1, jobz, 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_etb(i,i,1:8), work_etb, lwork_etb, rwork_etb, info)

    write(17,100) kfin + (2*pi/(3*sqrt(3.)*a)-ky(i))*a/pi, kx(i)*a/pi, ky(i)*a/pi, energy_band_etb(i,i,1:8)
end do

kfin = kfin + 2*pi/(3*sqrt(3.)*a)*a/pi

! M -> Gamma =====================================================
do i = 1,np

    kx(i) = 2*pi/(3*a) - i*2*pi/(3*a*np)
    ky(i) = 0.

! TB ================================

    matr_H_tb = matrix_H(kx(i),ky(i))
    matr_S_tb = matrix_S(kx(i),ky(i))

    call chegv(1, jobz, 'U', 2, matr_H_tb, 2, matr_S_tb, 2, energy_band_tb(i,i,1:2), work, lwork, rwork, info)

    write(18,100) kfin + (2*pi/(3*a)-kx(i))*a/pi, kx(i)*a/pi, ky(i)*a/pi, energy_band_tb(i,i,1:2)

! ETB ===================================

    matr_H_etb = matrix_H_etb(kx(i),ky(i))
    matr_S_etb = matrix_S_etb(kx(i),ky(i))

    call chegv(1, jobz, 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_etb(i,i,1:8), work_etb, lwork_etb, rwork_etb, info)

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
    call tube_geom(ttube,kcx,kcy,ktx,kty,nc,nt)

allocate(energy_band_cnt(nc,nt+1,8))
allocate(wf_cnt(nc,nt+1,8,8))
allocate(kx(nt+1),ky(nt+1))
allocate(overlap(nc,nt+1,8,8))

jobz='V'

do kc = 1 - nc/2, nc/2
    do kt = - nt/2, nt/2
        kx(kt) = kc*kcx + kt*ktx
        ky(kt) = kc*kcy + kt*kty

        matr_H_etb = matrix_H_etb(kx(kt),ky(kt))
        matr_S_etb = matrix_S_etb(kx(kt),ky(kt))
        overlap(kc,kt,1:8,1:8) = matr_S_etb

        ! solve generalized eigenvalue problem
        call chegv(1, jobz, 'U', 8, matr_H_etb, 8, matr_S_etb, 8, energy_band_cnt(kc,kt,1:8), work_etb, lwork_etb, rwork_etb, info)

        ! save eigen vectors
        if (jobz.eq.'V') then
            wf_cnt(kc,kt,1:8,1:8) = matr_H_etb
        endif

    enddo
enddo

    call order(nc, nt, energy_band_cnt, wf_cnt)

do kc=1-nc/2, nc/2
    do kt=-nt/2, nt/2
        ! write energy bands
        write(25,100) kt*sqrt(ktx**2 + kty**2)*ttube, kt*ktx*ttube, kt*kty*ttube, energy_band_cnt(kc,kt,1:8)

        ! write wave functions coefficients
        do i=1,8
           write(26,100) kt*sqrt(ktx**2 + kty**2)*ttube, kx(kt)*ttube, ky(kt)*ttube, wf_cnt(kc,kt,i,1:8)
        end do
    enddo
    write(25,*) ' '
    write(26,*) ' '
enddo


close(25)
close(26)
deallocate(energy_band_cnt)
deallocate(wf_cnt)
deallocate(kx,ky)
deallocate(work, work_etb)
deallocate(rwork, rwork_etb)
deallocate(overlap)

end program graphene_tb_prog

!-------------------------------------------------------------------------------

    include '/Users/dariasatco/Documents/eclipse_projects/georgiifortran/tbdftsnt/entube.f'
    include '/Users/dariasatco/Documents/eclipse_projects/georgiifortran/tbdftsnt/tbtube.f'
    include '/Users/dariasatco/Documents/eclipse_projects/georgiifortran/tbdftsnt/tbpot.f'
    include '/Users/dariasatco/Documents/eclipse_projects/georgiifortran/cntgrlib/gcd.f'
    include '/Users/dariasatco/Documents/eclipse_projects/georgiifortran/cntgrlib/tubpar.f'

!-------------------------------------------------------------------------------

!subroutine to generate geometry of cnt reciprocal space
    subroutine tube_geom(t,kcx,kcy,ktx,kty,nc,nt)

    implicit none

    real*8 l
    parameter (l=+0.3d+02)
    real*8 pi,sr3
    parameter (pi=+0.3141592653589793d+01,sr3=+0.1732050807568877d+01)!sr3 means sqrt(3)
    real*8 acc,auc
    parameter (acc=+0.142d-9,auc=acc*sr3) !acc = interatomic distance, auc = sqrt(3)*acc
    integer:: n,m,nc,nt
    real*8 t,d,theta,tab,phiab,dab,a1t,a1phi,a2t,a2phi,kcx,kcy,ktx,kty
    integer:: mult

! tube indices (n,m)
    print*, 'Enter tube indices (n,m):'
    print*, 'n:'
    read*, n
    print*, 'm:'
    read*, m

    call tubpar (n,m,l,t,d,theta,tab,phiab,dab,a1t,a1phi,a2t,a2phi,nc,nt,kcx,kcy,ktx,kty)

    !change the units of k to the m^(-1), initially the were in 1/auc
    kcx=kcx/auc
    kcy=kcy/auc
    !enlarge the number of mesh points
    mult=10
    nt=nt*mult
    ktx=ktx/(mult*auc)
    kty=kty/(mult*auc)
    !change the units of t to m (it was in nm)
    t=t*1.e-9
    if (mod(nt,2).ne.0) print*, 'nt/2 is odd!'
    print*, 'nc=', nc
    print*, 'nt=', nt

    return
    end subroutine tube_geom

