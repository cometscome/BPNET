module mod_bppotential
    use iso_fortran_env
    use bfio,only:FILEPATHLEN
    use fingerprint_basis,only:FingerPrintBasis
    use chebyshevbasis_fp,only:FP_ChebyshevBasis
    use module_multi_layers,only:Multi_layers
    use module_aenet_potential,only:aenet_potential
    use module_multi_layers,only:Multi_layers
    use ioutils,only:to_lowercase
   ! use fingerprints,only:initialize_fingerprint_basisset_fromparam
    !use mod_networks,only:Network
    !use MLP,only:MLPnet

    implicit none

    integer, parameter :: NENV_MAX = 2
    
    type, public :: BPpotential
        character(len=FILEPATHLEN) :: potential_filename
        character(len=2) :: typeName
        !class(Network),allocatable::net


        integer                                             :: neval
        character(len=1024)                                 :: description
        character(len=2)                              :: atomtype
        integer                                             :: nenv
        character(len=2), dimension(:),   allocatable :: envtypes
   
        integer                                             :: ntypes_global
        integer,                dimension(:),   allocatable :: gtype
        integer,                dimension(:),   allocatable :: ltype
   
        double precision                                    :: Rc_min
        double precision                                    :: Rc_max
   
        character(len=100)                                  :: sftype
        integer                                             :: nsf
        integer,                dimension(:),   allocatable :: sf
        integer                                             :: nsfparam
        double precision,       dimension(:,:), allocatable :: sfparam
        integer,                dimension(:,:), allocatable :: sfenv
   
        double precision,       dimension(:),   allocatable :: sfval_min
        double precision,       dimension(:),   allocatable :: sfval_max
        double precision,       dimension(:),   allocatable :: sfval_avg
        double precision,       dimension(:),   allocatable :: sfval_cov

        character(len=FILEPATHLEN)                            :: training_filename
        logical                                           :: normalized
        double precision                                  :: E_scale, shift

        integer                                           :: nTypes
        integer                                           :: nAtomsTot
        character(len=2), dimension(:), allocatable :: typeNames
        double precision,       dimension(:), allocatable :: E_atom
        integer                                           :: nStrucs
   
        double precision                                  :: E_min, E_max, E_av

        type(Multi_layers)::layers
  

        class(FingerPrintBasis), allocatable  ::  fingerprint_basis


        contains 
        procedure :: display => BPpotential_display
        procedure :: forward => BPpotential_forward
        procedure :: backward => BPpotential_backward
        procedure:: set_global_types
        procedure::  normalize => normalize_fingerprints
        procedure::deallocate
    end type

    interface BPpotential
        module procedure init_BPpotential
    end interface BPpotential


    contains 

    type(BPpotential) function init_BPpotential(typeName,potential_filename) result(BPpot)
        implicit none
        character(len=*),intent(in) ::potential_filename
        character(len=*),intent(in) ::typeName
        integer::u
        type(Multi_layers)::net_MLP
        type(aenet_potential)::potential_params
        !type(MLPnet):: net_MLP

        real(real64),allocatable::xin(:),xout(:)
        real(real64),allocatable::dEdx(:),dEdE(:)

        real(real64),allocatable::xinp(:),xoutp(:)
        integer::i
        real(real64)::eta
        integer                      :: ntypes_global
        real(real64)::scale
        

        BPpot%potential_filename = potential_filename
        BPpot%typeName = typeName

        !write(*,*) BPpot%typeName,BPpot%potential_filename
    
        open(newunit= u, file=trim(potential_filename), status='old', action='read', &
            form='unformatted')

        potential_params = aenet_potential(u)
        
        read(u) BPpot%description
        read(u) BPpot%atomtype
        read(u) BPpot%nenv
        allocate(BPpot%envtypes(BPpot%nenv))
        read(u) BPpot%envtypes
        read(u) BPpot%Rc_min
        read(u) BPpot%Rc_max
        read(u) BPpot%sftype
        read(u) BPpot%nsf
        read(u) BPpot%nsfparam
        allocate(BPpot%sf(BPpot%nsf),                    &
                    BPpot%sfparam(BPpot%nsfparam, BPpot%nsf), &
                    BPpot%sfenv(NENV_MAX,BPpot%nsf),        &
                    BPpot%sfval_min(BPpot%nsf),             &
                    BPpot%sfval_max(BPpot%nsf),             &
                    BPpot%sfval_avg(BPpot%nsf),             &
                    BPpot%sfval_cov(BPpot%nsf)              )
        read(u) BPpot%sf(:)
        read(u) BPpot%sfparam(:,:)
        read(u) BPpot%sfenv(:,:)
        read(u) BPpot%neval
        read(u) BPpot%sfval_min(:)
        read(u) BPpot%sfval_max(:)
        read(u) BPpot%sfval_avg(:)
        read(u) BPpot%sfval_cov(:)

        !write(*,*) BPpot%sfval_cov(:)


        read(u) BPpot%training_filename
        read(u) BPpot%normalized
        read(u) scale
        BPpot%E_scale = 1d0/scale
        read(u) BPpot%shift
        read(u) BPpot%nTypes
        allocate(BPpot%typeNames(BPpot%nTypes), BPpot%E_atom(BPpot%nTypes))
        read(u) BPpot%typeNames(1:BPpot%nTypes)
        read(u) BPpot%E_atom(1:BPpot%nTypes)
        read(u) BPpot%nAtomsTot
        read(u) BPpot%nStrucs
        read(u) BPpot%E_min, BPpot%E_max, BPpot%E_av

        write(*,*) BPpot%training_filename
        write(*,*)  BPpot%E_min, BPpot%E_max, BPpot%E_av
        close(u)

        BPpot%layers = Multi_layers(potential_params)

        ntypes_global = size(BPpot%typeNames(:))
        
        allocate( BPpot%gtype(BPpot%nenv), &
                    BPpot%ltype(ntypes_global))
        call BPpot%set_global_types(ntypes_global, BPpot%typeNames)

        select case(trim(to_lowercase(BPpot%sftype)))
        case('chebyshev')
           allocate(BPpot%fingerprint_basis,&
            source=FP_ChebyshevBasis(BPpot%sfparam,&
                    BPpot%nenv,&
                    BPpot%envtypes))
        case('behler2011')
            write(*,*) "behler2011 is not implemeted yet!!"
        case default
            write(0,*) "Error: Unknown basis function type : ", trim(BPpot%sftype)
            stop
        end select

        !write(*,*) BPpot%fingerprint_basis%get_num_parameters()


        !call test()

        contains 
        
        subroutine test()
            implicit none
            call BPpot%layers%display()
            allocate(xin(1:potential_params%nnodes(1)),source=0d0)
            allocate(xout(1:potential_params%nnodes(potential_params%nlayers)),source=0d0)
            call random_number(xin)
            call BPpot%layers%forward(xin,xout)
            !
            write(*,*) xout
    
            allocate(dEdx(1:potential_params%nnodes(1)),source=0d0)
            allocate(dEdE(1:potential_params%nnodes(potential_params%nlayers)),source=0d0)
            dEdE(1) =1d0
            call BPpot%layers%backward(dEdE,dEdx)
    
    
            allocate(xinp(1:potential_params%nnodes(1)),source=0d0)
            allocate(xoutp(1:potential_params%nnodes(potential_params%nlayers)),source=0d0)
            eta = 1d-7
            do i=1,potential_params%nnodes(1)
                xinp(1:potential_params%nnodes(1)) = xin(1:potential_params%nnodes(1))
                xinp(i) = xinp(i) + eta
                !write(*,*) "before",xin
                !write(*,*) "after ",xinp
                call BPpot%layers%forward(xinp,xoutp)
                !write(*,*) xoutp,xout
                write(*,*) i,(xoutp-xout)/eta,dEdx(i),(xoutp-xout)/eta-dEdx(i)
            end do
            
        end subroutine 
        !write(*,*) dEdx

    end function init_BPpotential

    subroutine deallocate(self)
        implicit none
        class(BPpotential),intent(inout)::self

        if (allocated(self%envtypes)) deallocate(self%envtypes)
        

        if (allocated(self%sf)) deallocate(self%sf)
        if (allocated(self%sfparam)) deallocate(self%sfparam)
        if (allocated(self%sfenv)) deallocate(self%sfenv)
        if (allocated(self%sfval_min)) deallocate(self%sfval_min)
        if (allocated(self%sfval_max)) deallocate(self%sfval_max)
        if (allocated(self%sfval_avg)) deallocate(self%sfval_avg)
        if (allocated(self%sfval_cov)) deallocate(self%sfval_cov)
        if (allocated(self%typeNames)) deallocate(self%typeNames)
        if (allocated(self%E_atom)) deallocate(self%E_atom)
        if (allocated(self%gtype)) deallocate(self%gtype)
        if (allocated(self%ltype)) deallocate(self%ltype)
        if (allocated(self%fingerprint_basis)) deallocate(self%fingerprint_basis)

    end subroutine


    subroutine normalize_fingerprints(self,values,sfderiv_i,sfderiv_j,n)
        implicit none
        class(BPpotential),intent(in)::self
        real(real64),intent(inout)::values(:)
        real(real64), dimension(:,:),   optional, intent(inout) :: sfderiv_i
        real(real64), dimension(:,:,:), optional, intent(inout) :: sfderiv_j
        integer,                            optional, intent(in)    :: n
        integer::i
        real(real64) ::scale, shift, s
        logical          :: do_deriv

        if (present(sfderiv_i) .and. present(sfderiv_j) .and. present(n)) then
            do_deriv = .true.
         else
            do_deriv = .false.
         end if

        do i=1,self%nsf
            shift = self%sfval_avg(i)
            s = sqrt(self%sfval_cov(i) - shift*shift)
            scale = 1.0d0/s

            values(i) = scale*(values(i)-shift)
            if (do_deriv) then
                sfderiv_i(1:3,i) = scale*sfderiv_i(1:3,i)
                sfderiv_j(1:3,i,1:n) = scale*sfderiv_j(1:3,i,1:n)
            end if
        end do

    end subroutine

    subroutine BPpotential_display(self)
        implicit none
        class(BPpotential),intent(in)::self
        call self%layers%display()
    end subroutine

    subroutine BPpotential_forward(self,xin,xout)
        use iso_fortran_env
        implicit none
        class(BPpotential),intent(inout)::self
        real(real64),intent(in)::xin(:)
        real(real64),intent(out)::xout(:)
        call self%layers%forward(xin,xout)
    end subroutine

    subroutine BPpotential_backward(self,dLdy,dLdx)
        use iso_fortran_env
        class(BPpotential), intent(inout) :: self
        real(real64),intent(in)::dLdy(:)
        real(real64),intent(out)::dLdx(:)
        call self%layers%backward(dLdy,dLdx)
    end subroutine

    subroutine set_global_types(self, ntypes_global, global_types)

        implicit none
        class(BPpotential) ::self
        integer,                                    intent(in)    :: ntypes_global
        character(len=*), dimension(ntypes_global), intent(in)    :: global_types
    
        integer :: i, j
    
        ! atom type indices that have no corresponding entry are set to 0
    
        do i = 1, ntypes_global
           self%ltype(i) = 0
           env : do j = 1, self%nenv
              if (trim(global_types(i)) == trim(self%envtypes(j))) then
                 self%ltype(i) = j
                 exit env
              end if
           end do env
        end do
    
        ! reverse direction, because we do not know, if the sets of types
        ! are identical
        do i = 1, self%nenv
           self%gtype(i) = 0
           global : do j = 1, ntypes_global
              if (trim(global_types(j)) == trim(self%envtypes(i))) then
                 self%gtype(i) = j
                 exit global
              end if
           end do global
        end do
    
    end subroutine set_global_types
end module mod_bppotential