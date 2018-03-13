module energy
    implicit none

    real, parameter :: t = -3.033   !hopping integral
    real, parameter :: s = 0.129    !overlap
    real, parameter :: a = 1.42*1e-10   !C-C bonding length
    real, parameter :: pi = 3.14159
    real, parameter :: e_p = 0.0    !pi-electron energy
    real, parameter :: delta = 1.0
    complex, parameter :: Im_i = cmplx(0.0, 1.0)    !imaginary unity

contains

    function energy_tb(kx, ky)

        real:: kx, ky, g
        real:: energy_tb(2)

        g=sqrt(1. + 4*cos(a*sqrt(3.)*ky/2)**2 + 4*cos(a*sqrt(3.)*ky/2)*cos(a*3*kx/2))

        energy_tb(1)=t*g/(1+s*g)
        energy_tb(2)=-t*g/(1-s*g)

    end function

end module energy
