module mod_networks
    implicit none

    type, abstract,public :: Network
        contains
        procedure(testsub), deferred :: test
    end type Network

    abstract interface
        function testsub(this) result(res)
            import :: Network
            class(Network), intent(in) :: this
            real :: res
        end function testsub
    end interface


end module mod_networks