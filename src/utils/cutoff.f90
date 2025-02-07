module cutoffmodule
    use iso_fortran_env
    implicit none

    real(real64), parameter, private :: PI     = 3.14159265358979d0
    real(real64), parameter, private :: PI_INV = 1.0d0/PI
    real(real64), parameter, private :: PI2    = 2.0d0*PI
    real(real64), parameter, private :: EPS    = 1.0d-12

    contains

    real(real64) pure function cutoff_fc(Rij,Rc) result(fc)
        implicit none
        real(real64), intent(in) :: Rij, Rc

        if (Rij >= Rc) then
            fc  = 0.0d0
        else
            fc  =  0.5d0*(cos(PI/Rc*Rij) + 1.0d0)
        end if
    end function cutoff_fc

    real(real64) pure function cutoff_fc_d1(Rij,Rc) result(dfc)
        implicit none
        real(real64), intent(in) :: Rij, Rc
        real(real64) :: a

        if (Rij >= Rc) then
            dfc  = 0.0d0
        else
            a = PI/Rc
            dfc = -0.5d0*a*sin(a*Rij)
        end if
    end function cutoff_fc_d1
end module cutoffmodule