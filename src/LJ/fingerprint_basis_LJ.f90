module LJbasis_fp
    use iso_fortran_env
    use fingerprint_basis,only:FingerPrintBasis
    use chebyshev,only:chebyshev_polynomial,chebyshev_polynomial_d1
    use cutoffmodule,only:cutoff_fc,cutoff_fc_d1

    implicit none
    real(real64), parameter, private :: PI     = 3.14159265358979d0
    real(real64), parameter, private :: PI_INV = 1.0d0/PI
    real(real64), parameter, private :: PI2    = 2.0d0*PI
    real(real64), parameter, private :: EPS    = 1.0d-12


    type, extends(FingerPrintBasis) :: FP_LJBasis
        real(real64)                                       :: radial_Rc
        integer                                       :: num_of_different_species
        integer                                       :: num_of_parameters

        character(len=2), dimension(:),   allocatable :: species_of_surrounding_atom !atomtypes
        integer,          dimension(:),   allocatable :: typeid
        real(real64), dimension(:),   allocatable          :: typespin


        contains
        !procedure :: area => circle_area
        !procedure :: test => chebyshev_test
        procedure ::  get_num_parameters =>  get_num_parameters_LJ
        procedure :: evaluate => evaluate_values_of_fingerprints_LJ

    end type FP_LJBasis

    interface FP_LJBasis
        procedure::setup_basis_LJ
    end interface

    contains


    function get_num_parameters_LJ(self) result(num)
        class(FP_LJBasis), intent(in):: self
        integer :: num
        num = self%num_of_parameters
    end function get_num_parameters_LJ

    subroutine evaluate_values_of_fingerprints_LJ(self,itype0,cartesian_coordinates, num_of_neighbor_atoms, &
                                                        coordinates_of_neighbor_atoms, types_of_neighbor_atoms, ltype,&
                                                        values,&
                                                        deriv_i,deriv_j,scaled)
        implicit none
        class(FP_LJBasis), intent(inout) :: self
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
        logical ::do_deriv
        integer::i,j,k
        real(real64)::s_j,s_k
        real(real64)::R_ij(1:3),R_ik(1:3),R_jk(1:3)
        real(real64)::d_ij,d_ik,d_jk
        integer::istart,iend

        integer               :: type0_loc
        integer, dimension(num_of_neighbor_atoms) :: type1_loc
        

        !write(*,*) "num_of_neighbor_atoms: ",num_of_neighbor_atoms

        if(present(deriv_i) .and. present(deriv_j)) then
            do_deriv = .true.
            deriv_i(1:3,1:self%num_of_parameters) = 0d0
            deriv_j(1:3,1:self%num_of_parameters,1:num_of_neighbor_atoms) = 0d0
        else
            do_deriv = .false.
        endif

        ! convert global atom type IDs to setup local IDs
        type0_loc = ltype(itype0)
        do i = 1, num_of_neighbor_atoms
           type1_loc(i) = ltype(types_of_neighbor_atoms(i))
        end do

        for_j:&
        do j=1,num_of_neighbor_atoms
            R_ij = coordinates_of_neighbor_atoms(1:3,j) - cartesian_coordinates(1:3)
            d_ij = sqrt(dot_product(R_ij, R_ij))

            if ((d_ij <= self%radial_Rc) .and. (d_ij > EPS)) then
                values(1+2*(type1_loc(j)-1)) = values(1+2*(type1_loc(j)-1))   + 1d0/(d_ij**6)
                values(2+2*(type1_loc(j)-1)) = values(2+2*(type1_loc(j)-1))  + 1d0/(d_ij**12)
            end if
        end do for_j


    end subroutine evaluate_values_of_fingerprints_LJ




    type(FP_LJBasis) function setup_basis_LJ(fingerprint_parameters,num_of_different_species,&
                    species_of_surrounding_atom    ) result(fpbasis)
        implicit none
        real(real64),dimension(:,:),intent(in)::fingerprint_parameters
        integer,intent(in)::num_of_different_species
        character(len=2), dimension(:),intent(in)::species_of_surrounding_atom
        integer::i,s


        fpbasis%radial_Rc = fingerprint_parameters(1,1)
        !fpbasis%angular_Rc =  fingerprint_parameters(3,1)
        fpbasis%num_of_different_species = num_of_different_species

  
        fpbasis%num_of_parameters = 2*num_of_different_species


    
        allocate(fpbasis%species_of_surrounding_atom(num_of_different_species),      &
             fpbasis%typeid(num_of_different_species),          &
             fpbasis%typespin(num_of_different_species))

        fpbasis%species_of_surrounding_atom = species_of_surrounding_atom

        do i = 1, num_of_different_species
            fpbasis%typeid(i) = i
        end do


    end function setup_basis_LJ
    
end module LJbasis_fp