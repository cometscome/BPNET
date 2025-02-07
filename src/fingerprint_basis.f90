module fingerprint_basis
    use iso_fortran_env
    implicit none


    type, abstract,public :: FingerPrintBasis
        contains
        !procedure(testsub), deferred :: test
        procedure(func_get_num_parameters), deferred :: get_num_parameters
        !procedure(area_interface), deferred :: area
        procedure(evaluate_values_of_fingerprints), deferred :: evaluate
    end type FingerPrintBasis

    abstract interface
        !function testsub(this) result(res)
        !    import :: FingerprintBasis
        !    class(FingerprintBasis), intent(in) :: this
        !    real :: res
        !end function testsub

        function func_get_num_parameters(self) result(num)
            import :: FingerprintBasis
            class(FingerprintBasis), intent(in) :: self
            integer :: num
        end function func_get_num_parameters

        subroutine evaluate_values_of_fingerprints(self,itype0,cartesian_coordinates, num_of_neighbor_atoms, &
                    coordinates_of_neighbor_atoms, types_of_neighbor_atoms, ltype,&
                    values,&
                    deriv_i,deriv_j,scaled)
            use iso_fortran_env,only:real64
            import ::FingerPrintBasis
            class(FingerprintBasis), intent(inout) :: self
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
        end subroutine evaluate_values_of_fingerprints
    end interface




end module fingerprint_basis