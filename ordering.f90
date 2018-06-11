module ordering
    implicit none

contains

subroutine order(nc, nt, energy_band_cnt, wf_cnt)

! input
integer:: nc, nt
real:: energy_band_cnt(nc,nt+1,8), energy_temp(nc,nt+1,8)
complex :: wf_cnt(nc,nt+1,8,8), wf_temp(nc,nt+1,8,8)


complex:: deriv1_wf(8,8), deriv2_wf(8,8)
real:: deriv1_en(8), deriv2_en(8)
integer:: kc, kt, kt1, kc1
integer:: i,j,l,m


energy_temp(:,:,:) = 0.
wf_temp(:,:,:,:) = 0.

do kc = 1 - nc/2, nc/2
    do kc1 = kc, nc/2
        do kt = 1 - nt/2, nt/2 - 1

            do j=0,1

            do l = 1+4*j, 3+4*j
                do m = l+1, 4+4*j
                    ! check nearest energy points
                    if ( abs(energy_band_cnt(kc,kt,l) - energy_band_cnt(kc1,kt,m)) .le. 3.e-1 ) then
                        print*, kc, kt, kc1

                        ! calculate derivative for wf coefficients
                        do i=1,8
                            deriv1_wf(i,l) =  wf_cnt(kc,kt,i,l) - wf_cnt(kc,kt-1,i,l)
                            deriv2_wf(i,l) = wf_cnt(kc,kt+1,i,l) - wf_cnt(kc,kt,i,l)

                            deriv1_wf(i,m) =  wf_cnt(kc1,kt,i,m) - wf_cnt(kc1,kt-1,i,m)
                            deriv2_wf(i,m) = wf_cnt(kc1,kt+1,i,m) - wf_cnt(kc1,kt,i,m)
                        enddo

!                        print*, "l: ", "deriv1=", sum(abs(deriv1_wf(1:8,l))), "deriv2=", sum(abs(deriv2_wf(1:8,l))), &
!                        "ratio=", sum(abs(deriv1_wf(1:8,l)-deriv2_wf(1:8,l)))/max( sum(abs(deriv1_wf(1:8,l))), &
!                        sum(abs(deriv2_wf(1:8,l))) )
!
!                        print*, "m: ", "deriv1=", sum(abs(deriv1_wf(1:8,m))), "deriv2=", sum(abs(deriv2_wf(1:8,m))), &
!                        "ratio=", sum(abs(deriv1_wf(1:8,m)-deriv2_wf(1:8,m)))/max( sum(abs(deriv1_wf(1:8,m))), &
!                        sum(abs(deriv2_wf(1:8,m))) )

                        ! check continuity of wf
                        if (sum(abs(deriv1_wf(1:8,l)-deriv2_wf(1:8,l)))/max( sum(abs(deriv1_wf(1:8,l))), &
                        sum(abs(deriv2_wf(1:8,l))) ).ge. 0.2 &
                        .and. sum(abs(deriv1_wf(1:8,m)-deriv2_wf(1:8,m)))/max( sum(abs(deriv1_wf(1:8,l))),&
                        sum(abs(deriv2_wf(1:8,l))) ) .ge. 0.2) then

                        ! if it is not continuous, change vector places

                            ! check the position of interchange start
                            if (sum(abs(deriv1_wf(1:8,l))).ge.sum(abs(deriv2_wf(1:8,l)))) then
                            kt1 = kt
                            else
                            kt1 = kt+1
                            endif

                            ! create a copy of wf vector values
                            wf_temp(kc,kt1:nt/2,1:8,l) =  wf_cnt(kc,kt1:nt/2,1:8,l)
                            wf_temp(kc1,kt1:nt/2,1:8,m) =  wf_cnt(kc1,kt1:nt/2,1:8,m)

                           ! create a copy of energy values
                            energy_temp(kc,kt1:nt/2,l) = energy_band_cnt(kc,kt1:nt/2,l)
                            energy_temp(kc1,kt1:nt/2,m) = energy_band_cnt(kc1,kt1:nt/2,m)

                           ! interchange m-th and l-th wf vectors
                            wf_cnt(kc,kt1:nt/2,1:8,l) = wf_temp(kc1,kt1:nt/2,1:8,m)
                            wf_cnt(kc1,kt1:nt/2,1:8,m) = wf_temp(kc,kt1:nt/2,1:8,l)


                        ! calculate derivative for wf coefficients after interchange
                            do i=1,8
                                deriv1_wf(i,l) =  wf_cnt(kc,kt,i,l) - wf_cnt(kc,kt-1,i,l)
                                deriv2_wf(i,l) = wf_cnt(kc,kt+1,i,l) - wf_cnt(kc,kt,i,l)

                                deriv1_wf(i,m) =  wf_cnt(kc,kt,i,m) - wf_cnt(kc,kt-1,i,m)
                                deriv2_wf(i,m) = wf_cnt(kc,kt+1,i,m) - wf_cnt(kc,kt,i,m)
                            enddo

!                            print*, "after interchange"
!                            print*, "deriv1=", sum(abs(deriv1_wf(1:8,l))), "deriv2=", sum(abs(deriv2_wf(1:8,l))), &
!                            "ratio=", sum(abs(deriv1_wf(1:8,l)-deriv2_wf(1:8,l)))/max( sum(abs(deriv1_wf(1:8,l))), &
!                            sum(abs(deriv2_wf(1:8,l))) )

                        ! check continuity of wf again
                            if (sum(abs(deriv1_wf(1:8,l)-deriv2_wf(1:8,l)))/max( sum(abs(deriv1_wf(1:8,l))), &
                            sum(abs(deriv2_wf(1:8,l))) ).ge. 0.2 &
                            .or. sum(abs(deriv1_wf(1:8,m)-deriv2_wf(1:8,m)))/max( sum(abs(deriv1_wf(1:8,l))),&
                            sum(abs(deriv2_wf(1:8,l))) ) .ge. 0.2) then

                        ! if they are not continuous cancle permutation
                                print*, "wf is not continuous again!"

                                wf_cnt(kc,kt1:nt/2,1:8,l) = wf_temp(kc,kt1:nt/2,1:8,l)
                                wf_cnt(kc1,kt1:nt/2,1:8,m) = wf_temp(kc1,kt1:nt/2,1:8,m)

                                else
                        ! if they are continuous perform energies interchange
                                energy_band_cnt(kc,kt1:nt/2,l) = energy_temp(kc1,kt1:nt/2,m)
                                energy_band_cnt(kc1,kt1:nt/2,m) = energy_temp(kc,kt1:nt/2,l)
                                print*, 'changed'

                            endif

                       ! calculate derivatives for energies
                       deriv1_en(l) = energy_band_cnt(kc,kt,l) - energy_band_cnt(kc,kt-1,l)
                       deriv2_en(l) = energy_band_cnt(kc,kt+1,l) - energy_band_cnt(kc,kt,l)

                       deriv1_en(m) = energy_band_cnt(kc1,kt,m) - energy_band_cnt(kc1,kt-1,m)
                       deriv2_en(m) = energy_band_cnt(kc1,kt+1,m) - energy_band_cnt(kc1,kt,m)

                       if (abs(deriv1_en(l)-deriv2_en(l)).ge. 0.15 .and. &
                       abs(deriv1_en(m)-deriv2_en(m)).ge. 0.15) then

                           print*, "energy is not continuous!", "l = ", l,  "m = ", m
!                           print*, "l: ", deriv1_en(l), deriv2_en(l), "m: ", deriv1_en(m), deriv2_en(m)

                       ! if energy is not continuous make interchange of l-th and m-th vectors
                           wf_cnt(kc,kt1:nt/2,1:8,l) = wf_temp(kc,kt1:nt/2,1:8,l)
                           wf_cnt(kc1,kt1:nt/2,1:8,m) = wf_temp(kc1,kt1:nt/2,1:8,m)

                           energy_band_cnt(kc,kt1:nt/2,l) = energy_temp(kc1,kt1:nt/2,m)
                           energy_band_cnt(kc1,kt1:nt/2,m) = energy_temp(kc,kt1:nt/2,l)

!                           ! once again
!                           deriv1_en(l) = energy_band_cnt(kc,kt,l) - energy_band_cnt(kc,kt-1,l)
!                           deriv2_en(l) = energy_band_cnt(kc,kt+1,l) - energy_band_cnt(kc,kt,l)
!
!                           deriv1_en(m) = energy_band_cnt(kc1,kt,m) - energy_band_cnt(kc1,kt-1,m)
!                           deriv2_en(m) = energy_band_cnt(kc1,kt+1,m) - energy_band_cnt(kc1,kt,m)
!
!                           print*, "after reverse"
!                           print*, "l: ", deriv1_en(l), deriv2_en(l), "m: ", deriv1_en(m), deriv2_en(m)

                       endif
                    endif
                  endif
                enddo
            enddo
            enddo
        enddo
     enddo
enddo

return

end subroutine


end module ordering
