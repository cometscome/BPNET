module MLP
    use iso_fortran_env
    use mod_networks,only:Network
    implicit none

    type, public,extends(Network) :: MLPnet
        integer                                       :: nlayers
        integer,          dimension(:),   allocatable :: nnodes
        integer                                       :: nnodes_max

        integer,          dimension(:),   allocatable :: f_a

        double precision, dimension(:),   allocatable :: W
        integer,          dimension(:),   allocatable :: iw
        integer                                       :: Wsize

        integer,          dimension(:),   allocatable :: iv
        integer                                       :: nvalues
    contains
        procedure :: test => MLP_test
    end type

    interface MLPnet
        procedure::setup_MLP
    end interface 

    contains



    function MLP_test(this) result(res)
        class(MLPnet), intent(in):: this
        real :: res
        res = 0d0
    end function MLP_test   
    

    type(MLPnet) function setup_MLP(u) result(net)
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

        
    end function
    


end module MLP