module Multibasis_fp
    use iso_fortran_env
    use fingerprint_basis,only:FingerPrintBasis
    use chebyshevbasis_fp,only:FP_ChebyshevBasis
    use splinebasis_fp,only:FP_SplineBasis
    use LJbasis_fp,only:FP_LJBasis

    implicit none

    type :: FPWrapper
        class(FingerPrintBasis), allocatable :: basisset
    end type

    type, extends(FingerPrintBasis) :: FP_MultiBasis
        type(FPWrapper),allocatable::basis(:)
        integer                                       :: num_of_different_species
        integer                                       :: num_of_parameters
        integer::num_kinds

        !real(real64), dimension(:), allocatable :: temp_values
        !real(real64), dimension(:,:), allocatable :: temp_deriv_i
        !real(real64), dimension(:,:,:), allocatable :: temp_deriv_j

        contains
        procedure ::  get_num_parameters =>  get_num_parameters_multi
        procedure :: evaluate => evaluate_values_of_fingerprints_multi

    end type

    interface FP_MultiBasis
        procedure::setup_basis_multi
    end interface

    contains

    function get_num_parameters_multi(self) result(num)
        class(FP_MultiBasis), intent(in):: self
        integer :: num
        num = self%num_of_parameters
    end function get_num_parameters_multi

    subroutine evaluate_values_of_fingerprints_multi(self,itype0,cartesian_coordinates, num_of_neighbor_atoms, &
                                                        coordinates_of_neighbor_atoms, types_of_neighbor_atoms, ltype,&
                                                        values,&
                                                        deriv_i,deriv_j,scaled)
        implicit none
        class(FP_MultiBasis), intent(inout) :: self
        integer,                                      intent(in)    :: itype0
        real(real64), dimension(3),               intent(in)    :: cartesian_coordinates
        integer,                                      intent(in)    :: num_of_neighbor_atoms
        real(real64), dimension(3,num_of_neighbor_atoms),             intent(in)    :: coordinates_of_neighbor_atoms
        integer,          dimension(num_of_neighbor_atoms),               intent(in)    :: types_of_neighbor_atoms
        integer,                dimension(:),   intent(in) :: ltype
        real(real64), dimension(:),               intent(out)   :: values
        real(real64), dimension(:,:),   optional, intent(out)   :: deriv_i
        real(real64), dimension(:,:,:), optional, intent(out)   :: deriv_j
        logical,                            optional, intent(in)    :: scaled 
        integer::ikind
        integer::num_params_i
        integer::istart
        logical ::do_deriv
        


        if(present(deriv_i) .and. present(deriv_j)) then
            do_deriv = .true.
            deriv_i(1:3,1:self%num_of_parameters) = 0d0
            deriv_j(1:3,1:self%num_of_parameters,1:num_of_neighbor_atoms) = 0d0
        else
            do_deriv = .false.
        endif


        istart = 1
        do ikind=1,self%num_kinds
            num_params_i = self%basis(ikind)%basisset%get_num_parameters()

            if(do_deriv) then
                call self%basis(ikind)%basisset%evaluate(itype0,cartesian_coordinates, num_of_neighbor_atoms, &
                                                            coordinates_of_neighbor_atoms, types_of_neighbor_atoms, ltype,&
                                                            values(istart:istart + num_params_i-1),&
                                                            deriv_i(1:3,istart:istart + num_params_i-1),&
                                                            deriv_j(1:3,istart:istart + num_params_i-1,1:num_of_neighbor_atoms),&
                                                            scaled)
            else
                call self%basis(ikind)%basisset%evaluate(itype0,cartesian_coordinates, num_of_neighbor_atoms, &
                                                        coordinates_of_neighbor_atoms, types_of_neighbor_atoms, ltype,&
                                                        values(istart:istart + num_params_i-1))
            end if
            !write(*,*) values(istart:istart + num_params_i-1)
            istart = istart + num_params_i
            
        end do

        !write(*,*) "multi basis",values
        !stop
    end 

    type(FP_MultiBasis) function setup_basis_multi(fingerprint_parameters,num_of_different_species,&
            species_of_surrounding_atom,num_of_maxparameters    ) result(fpbasis)
        implicit none
        real(real64),dimension(:,:),intent(in)::fingerprint_parameters
        integer,intent(in)::num_of_different_species,num_of_maxparameters
        character(len=2), dimension(:),intent(in)::species_of_surrounding_atom
        integer::ikind
        integer::ibasis
        type(FP_ChebyshevBasis)::chebyshev_basis 
        type(FP_SplineBasis)::spline_basis
        type(FP_LJBasis)::LJ_basis
        integer::istart,numparams_i
        !integer::num_of_parameters

        fpbasis%num_kinds =  nint(fingerprint_parameters(2,1))
        allocate(fpbasis%basis(fpbasis%num_kinds))
        fpbasis%num_of_parameters = 0
        fpbasis%num_of_different_species = num_of_different_species

        do ikind = 1,fpbasis%num_kinds
            istart = 3+(ikind-1)*(num_of_maxparameters+2)
            ibasis = nint(fingerprint_parameters(istart,1) )
            numparams_i = nint(fingerprint_parameters(istart+1,1) )
            if (ibasis == 1) then !chebyshev
                chebyshev_basis = FP_ChebyshevBasis(fingerprint_parameters(istart+2:istart+num_of_maxparameters+1,:),&
                    num_of_different_species,&
                    species_of_surrounding_atom)

                allocate(fpbasis%basis(ikind)%basisset,&
                        source= chebyshev_basis &
                        )
            elseif (ibasis == 2) then !spline
                spline_basis = FP_SplineBasis(fingerprint_parameters(istart+2:istart+num_of_maxparameters+1,:),&
                    num_of_different_species,&
                    species_of_surrounding_atom)

                allocate(fpbasis%basis(ikind)%basisset,&
                    source= spline_basis &
                    )
            elseif (ibasis == 3) then !LJ
                LJ_basis = FP_LJBasis(fingerprint_parameters(istart+2:istart+num_of_maxparameters+1,:),&
                    num_of_different_species,&
                    species_of_surrounding_atom)

                allocate(fpbasis%basis(ikind)%basisset,&
                    source= LJ_basis &
                    )
            else
                write(*,*) "ibasis = ",ibasis," is not supoorted in multibasis."
                stop 
            end if
            fpbasis%num_of_parameters = fpbasis%num_of_parameters  + fpbasis%basis(ikind)%basisset%get_num_parameters()
        end do

        !allocate(fpbasis%temp_values(fpbasis%num_of_parameters),source=0d0)
        !allocate(fpbasis%temp_deriv_i(1:3,fpbasis%num_of_parameters),source=0d0)
        !allocate(fpbasis%temp_deriv_j(1:3,fpbasis%num_of_parameters,1:1000),source=0d0)

    end

end module 