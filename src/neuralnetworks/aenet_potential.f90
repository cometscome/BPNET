module module_aenet_potential
    implicit none

    type, public:: aenet_potential
        integer                                       :: nlayers
        integer,          dimension(:),   allocatable :: nnodes
        integer                                       :: nnodes_max

        integer,          dimension(:),   allocatable :: f_a

        double precision, dimension(:),   allocatable :: W
        integer,          dimension(:),   allocatable :: iw
        integer                                       :: Wsize

        integer,          dimension(:),   allocatable :: iv
        integer                                       :: nvalues
    end type 

    interface aenet_potential
        procedure::setup_aenet_potential
    end interface 

    contains


    type(aenet_potential) function setup_aenet_potential(u) result(net)
        implicit  none
        integer,intent(in) ::u

        read(u) net%nlayers
        read(u) net%nnodes_max
        read(u) net%Wsize
        read(u) net%nvalues
    
        allocate(net%nnodes(net%nlayers), &
                 net%f_a(2:net%nlayers),  &
                 net%iw(1:net%nlayers),   &
                 net%iv(1:net%nlayers),   &
                 net%W(net%Wsize))
    
        read(u) net%nnodes(:)
        read(u) net%f_a(:)
        read(u) net%iw(:)
        read(u) net%iv(:)
        read(u) net%W(:)

        
    end function setup_aenet_potential

!   0 : linear function f(x) = x                                     !
!   1 : hyperbolic tangent, y in [-1:1]                              !
!   2 : sigmoid,            y in [ 0:1]                              !
!   3 : modified tanh,      y in [-1.7159:1.7159]  f(+/-1) = +/-1    !
!   4 : tanh & linear twisting term    
    character(len=128) function choose_activationfunction_fromaenet(fa) result(name)
        integer,intent(in)::fa
        if (fa == 0) then
            name = "linear"
        elseif(fa == 1) then
            name = "tanh"
        elseif(fa == 2) then
            name = "sigmoid"
        elseif(fa == 3) then
            name = "modifiedtanh"
        elseif(fa == 4) then
            name = "twist"
        endif
    end


end module module_aenet_potential