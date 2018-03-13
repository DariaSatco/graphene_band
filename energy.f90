module energy
    implicit none

    real, parameter :: t = -3.033
    real, parameter :: s = 0.129
    real, parameter :: a=1.42*1e-10
    real, parameter :: pi=3.14

contains

    function energy_tb(kx, ky)

        real:: kx, ky
        real:: g
        real:: energy_tb(2)

        g=sqrt(1. + 4*cos(a*sqrt(3.)*ky/2)**2 + 4*cos(a*sqrt(3.)*ky/2)*cos(a*3*kx/2))
        energy_tb(1)=t*g/(1+s*g)
        energy_tb(2)=-t*g/(1-s*g)

    end function

end module energy
