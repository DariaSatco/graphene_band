module matrix_tb
use energy
    implicit none

contains

    function matrix_H(kx, ky)

        real:: kx, ky
        complex :: g
        complex:: matrix_H(2,2)

        g = exp(-Im_i*a*kx) + 2*exp(Im_i*a*kx/2)*cos(a*sqrt(3.)*ky/2)

        matrix_H(1,1) = cmplx(e_p + delta, 0.)
        matrix_H(2,2) = cmplx(e_p - delta, 0.)
        matrix_H(1,2) = t*g
        matrix_H(2,1) = t*conjg(g)

    end function

!*******************************************

    function matrix_S(kx, ky)

        real:: kx, ky
        complex :: g
        complex:: matrix_S(2,2)

        g = exp(-Im_i*a*kx) + 2*exp(Im_i*a*kx/2)*cos(a*sqrt(3.)*ky/2)

        matrix_S(1,1) = cmplx(1., 0.)
        matrix_S(2,2) = cmplx(1.,0.)
        matrix_S(1,2) = s*g
        matrix_S(2,1) = s*conjg(g)

    end function

end module matrix_tb
