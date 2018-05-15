module matrix_etb
use energy
    implicit none
    !Slater-Koster parameters
    real, parameter:: Es = -7.3 !eV on-site energy for s-electrons
    real, parameter:: sssigma = -4.3 !eV
    real, parameter:: spsigma = 4.98 !eV
    real, parameter:: ppsigma = 6.38 !eV
    real, parameter:: pppi = -2.8 !eV

contains

    function matrix_H_etb(kx,ky)

    real:: kx, ky
    complex:: matrix_H_etb(8,8)
    complex:: g
    integer:: i,j
    complex:: M11(4,4), M12(4,4)
    complex:: AsBs, AsBx, AsBy, AxBx, AxBy, AyBy, AzBz

    g = exp(-Im_i*a*kx) + 2*exp(Im_i*a*kx/2)*cos(a*sqrt(3.)*ky/2)

    M11(1,1) = Es
    M11(1,2) = 0.
    M11(1,3) = 0.
    M11(1,4) = 0.

    do i=2,4
        do j=1,4
        M11(i,j)=0.
        end do
    end do

    AsBs = sssigma*g
    AsBx = spsigma*(-exp(-Im_i*a*kx) + exp(Im_i*a*kx/2)*cos(a*sqrt(3.)*ky/2))
    AsBy = spsigma*Im_i*sqrt(3.)*exp(Im_i*a*kx/2)*sin(a*sqrt(3.)*ky/2)
    AxBx = ppsigma*exp(-Im_i*a*kx) + (1./2*ppsigma + 3./2*pppi)*exp(Im_i*a*kx/2)*cos(a*sqrt(3.)*ky/2)
    AxBy = sqrt(3.)/2*(ppsigma - pppi)*Im_i*exp(Im_i*a*kx/2)*sin(a*sqrt(3.)*ky/2)
    AyBy = pppi*exp(-Im_i*a*kx) + (3./2*ppsigma + 1./2*pppi)*exp(Im_i*a*kx/2)*cos(a*sqrt(3.)*ky/2)
    AzBz = pppi*g

    M12(1,1) = AsBs
    M12(1,2) = AsBx
    M12(1,3) = AsBy
    M12(1,4) = 0.

    M12(2,1) = -conjg(AsBx)
    M12(2,2) = AxBx
    M12(2,3) = AxBy
    M12(2,4) = 0.

    M12(3,1) = -conjg(AsBy)
    M12(3,2) = -conjg(AxBy)
    M12(3,3) = AyBy
    M12(3,4) = 0.

    M12(4,1) = 0.
    M12(4,2) = 0.
    M12(4,3) = 0.
    M12(4,4) = AzBz

    do i=1,4
        do j=1,4

        matrix_H_etb(i,j)=M11(i,j)
        matrix_H_etb(i,j+4)=M12(i,j)
        matrix_H_etb(i+4,j)=conjg(M12(j,i))
        matrix_H_etb(i+4,j+4)=M11(i,j)

        end do
    end do

    end function

!--------------------------------------

    function matrix_S_etb(kx,ky)

    complex:: matrix_S_etb(8,8)
    real:: kx, ky
    integer:: i
    complex:: g

    g = exp(-Im_i*a*kx) + 2*exp(Im_i*a*kx/2)*cos(a*sqrt(3.)*ky/2)
    matrix_S_etb(:,:) = 0.

    do i=1,8
        matrix_S_etb(i,i) = 1.
    end do

    end function


end module matrix_etb
