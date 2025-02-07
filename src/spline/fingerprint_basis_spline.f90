module splinebasis_fp
    use iso_fortran_env
    use fingerprint_basis,only:FingerPrintBasis
    use chebyshev,only:chebyshev_polynomial,chebyshev_polynomial_d1
    use cutoffmodule,only:cutoff_fc,cutoff_fc_d1

    implicit none
    real(real64), parameter, private :: PI     = 3.14159265358979d0
    real(real64), parameter, private :: PI_INV = 1.0d0/PI
    real(real64), parameter, private :: PI2    = 2.0d0*PI
    real(real64), parameter, private :: EPS    = 1.0d-12


    type, extends(FingerPrintBasis) :: FP_SplineBasis
        integer                                       :: radial_points
        real(real64)                                       :: radial_Rc
        integer                                       :: radial_N
        integer                                       :: angular_points
        real(real64)                                       :: angular_Rc
        integer                                       :: angular_N
        integer                                       :: num_of_different_species
        real(real64), dimension(:),   allocatable :: radial_knots
        real(real64), dimension(:),   allocatable :: angular_knots

        integer                                       :: num_of_parameters

        integer                                       :: radial_1_initialindex, radial_1_finalindex
        integer                                       :: radial_2_initialindex, radial_2_finalindex
        integer                                       :: angular_1_initialindex, angular_1_finalindex
        integer                                       :: angular_2_initialindex, angular_2_finalindex
        logical                                       :: multi
        character(len=2), dimension(:),   allocatable :: species_of_surrounding_atom !atomtypes
        integer,          dimension(:),   allocatable :: typeid
        real(real64), dimension(:),   allocatable          :: typespin

        real(real64), dimension(:), allocatable :: temp_values
        real(real64), dimension(:,:), allocatable :: temp_deriv_i
        real(real64), dimension(:,:), allocatable :: temp_deriv_j
        real(real64), dimension(:,:), allocatable :: temp_deriv_k

        contains
        !procedure :: area => circle_area
        !procedure :: test => chebyshev_test
        procedure ::  get_num_parameters =>  get_num_parameters_spline
        procedure :: evaluate => evaluate_values_of_fingerprints_spline
        procedure :: evaluate_radial => evaluate_radial_spline
        procedure :: evaluate_angular => evaluate_angular_spline
    end type FP_SplineBasis

    interface FP_SplineBasis
        procedure::setup_basis_spline
    end interface

    contains

    !function chebyshev_test(this) result(res)
    !    class(FP_ChebyshevBasis), intent(in):: this
    !    real :: res
    !end function chebyshev_test

    function get_num_parameters_spline(self) result(num)
        class(FP_SplineBasis), intent(in):: self
        integer :: num
        num = self%num_of_parameters
    end function get_num_parameters_spline


    type(FP_SplineBasis) function setup_basis_spline(fingerprint_parameters,num_of_different_species,&
                    species_of_surrounding_atom    ) result(fpbasis)
        use bspline,only:d => bspline_order
        use bspline,only:make_knotsvector
        implicit none
        real(real64),dimension(:,:),intent(in)::fingerprint_parameters
        integer,intent(in)::num_of_different_species
        character(len=2), dimension(:),intent(in)::species_of_surrounding_atom
        integer::i,s


        fpbasis%radial_Rc = fingerprint_parameters(1,1)
        fpbasis%radial_points =  nint(fingerprint_parameters(2,1))
        fpbasis%angular_Rc =  fingerprint_parameters(3,1)

        fpbasis%angular_points = nint(fingerprint_parameters(4,1))
        fpbasis%num_of_different_species = num_of_different_species

        fpbasis%radial_N = fpbasis%radial_points +2*d -d -1
        !fpbasis%radial_N = fpbasis%radial_order + 1
        allocate(fpbasis%radial_knots(fpbasis%radial_points +2*d))
        call make_knotsvector(fpbasis%radial_knots,d,fpbasis%radial_points ,fpbasis%radial_Rc)
    
        fpbasis%angular_N= fpbasis%angular_points +2*d -d -1
        !fpbasis%angular_N = fpbasis%angular_points + 1
        !fpbasis%num_values = max(fpbasis%r_N, fpbasis%angular_N)
        allocate(fpbasis%angular_knots(fpbasis%angular_points +2*d))
        
        call make_knotsvector(fpbasis%angular_knots,d,fpbasis%angular_points ,2*fpbasis%angular_Rc)

        fpbasis%num_of_parameters = fpbasis%radial_N + fpbasis%angular_N
        fpbasis%radial_1_initialindex = 1
        fpbasis%radial_1_finalindex = fpbasis%radial_1_initialindex + fpbasis%radial_N - 1
        fpbasis%angular_1_initialindex = fpbasis%radial_1_finalindex + 1
        fpbasis%angular_1_finalindex = fpbasis%angular_1_initialindex + fpbasis%angular_N - 1
        fpbasis%radial_2_initialindex = fpbasis%angular_1_finalindex + 1
        fpbasis%radial_2_finalindex = fpbasis%radial_2_initialindex + fpbasis%radial_N - 1
        fpbasis%angular_2_initialindex = fpbasis%radial_2_finalindex + 1
        fpbasis%angular_2_finalindex = fpbasis%angular_2_initialindex + fpbasis%angular_N - 1

        if (fpbasis%num_of_different_species > 1) then
            fpbasis%multi = .true.
            fpbasis%num_of_parameters = 2*fpbasis%num_of_parameters
        else
            fpbasis%multi = .false.
        end if
    
        allocate(fpbasis%species_of_surrounding_atom(num_of_different_species),      &
             fpbasis%typeid(num_of_different_species),          &
             fpbasis%typespin(num_of_different_species))

        fpbasis%species_of_surrounding_atom = species_of_surrounding_atom

        do i = 1, num_of_different_species
            fpbasis%typeid(i) = i
        end do

        s = -num_of_different_species/2
        do i = 1, num_of_different_species
            if ((s == 0) .and. (mod(num_of_different_species, 2) == 0)) s = s + 1
            fpbasis%typespin(i) = dble(s)
            s = s + 1
        end do


        allocate(fpbasis%temp_values(fpbasis%num_of_parameters),source=0d0)
        allocate(fpbasis%temp_deriv_i(1:3,fpbasis%num_of_parameters),source=0d0)
        allocate(fpbasis%temp_deriv_j(1:3,fpbasis%num_of_parameters),source=0d0)
        allocate(fpbasis%temp_deriv_k(1:3,fpbasis%num_of_parameters),source=0d0)


    end function setup_basis_spline    

    subroutine evaluate_values_of_fingerprints_spline(self,itype0,cartesian_coordinates, num_of_neighbor_atoms, &
                                                        coordinates_of_neighbor_atoms, types_of_neighbor_atoms, ltype,&
                                                        values,&
                                                        deriv_i,deriv_j,scaled)
        implicit none
        class(FP_SplineBasis), intent(inout) :: self
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

       
        

        values(1:self%num_of_parameters) = 0d0
        s_j = 1d0

        !write(*,*) " values", values
        

        for_j:&
        do j=1,num_of_neighbor_atoms
            R_ij = coordinates_of_neighbor_atoms(1:3,j) - cartesian_coordinates(1:3)
            !d_ij = sqrt(R_ij(1)*R_ij(1)+R_ij(2)*R_ij(2)+R_ij(3)*R_ij(3))
            d_ij = sqrt(dot_product(R_ij, R_ij))
            
            !write(*,*) j,num_of_neighbor_atoms

            if ((d_ij <= self%radial_Rc) .and. (d_ij > EPS)) then
                istart = self%radial_1_initialindex
                iend = self%radial_1_finalindex

                if (do_deriv) then
                    call self%evaluate_radial(R_ij, d_ij, self%temp_values, &
                                deriv_i=self%temp_deriv_i, deriv_j=self%temp_deriv_j)
                    deriv_i(1:3, istart:iend) = deriv_i(1:3, istart:iend) + self%temp_deriv_i(1:3, 1:self%radial_N)

                    deriv_j(1:3, istart:iend, j) = deriv_j(1:3, istart:iend, j) + self%temp_deriv_j(1:3, 1:self%radial_N)
                else
                    call self%evaluate_radial(R_ij, d_ij, self%temp_values)
                end if
                values(istart:iend) = values(istart:iend) + self%temp_values(1:self%radial_N)

                istart = self%radial_2_initialindex
                iend = self%radial_2_finalindex

                if (self%multi) then
                    s_j = self%typespin(self%typeid(type1_loc(j)))

                    values(istart:iend) = values(istart:iend) + s_j*self%temp_values(1:self%radial_N)
                    if (do_deriv) then
                        deriv_i(1:3, istart:iend) = deriv_i(1:3, istart:iend) &
                                            + s_j*self%temp_deriv_i(1:3, 1:self%radial_N)
                        deriv_j(1:3, istart:iend, j) = deriv_j(1:3,istart:iend, j) &
                                            + s_j*self%temp_deriv_j(1:3, 1:self%radial_N)
                    endif
                endif
            endif

            if (d_ij > self%angular_Rc) cycle for_j
            for_k: &
            do k=j+1,num_of_neighbor_atoms
                
                R_ik = coordinates_of_neighbor_atoms(1:3,k) - cartesian_coordinates(1:3)
                d_ik = sqrt(dot_product(R_ik, R_ik))
                
                !d_ik = sqrt(R_ik(1)*R_ik(1)+R_ik(2)*R_ik(2)+R_ik(3)*R_ik(3))
                if ((d_ik > self%angular_Rc) .or. (d_ik < EPS)) cycle for_k
                R_jk = coordinates_of_neighbor_atoms(1:3,k)- coordinates_of_neighbor_atoms(1:3, j) 
                d_jk = sqrt(dot_product(R_jk, R_jk))

                !write(*,*) R_ij, R_ik
                !write(*,*) "d",d_ij, d_ik, cos_ijk
                !cos_ijk = dot_product(R_ij, R_ik)/(d_ij*d_ik)

                istart = self%angular_1_initialindex
                iend = self%angular_1_finalindex

                if (do_deriv) then
                    
                    call self%evaluate_angular(R_ij, R_ik, R_jk,d_ij, d_ik, d_jk, self%temp_values, &
                                deriv_i=self%temp_deriv_i, deriv_j=self%temp_deriv_j,deriv_k=self%temp_deriv_k)
                    !write(*,*) j,k,self%temp_values(1:self%angular_N)

                    deriv_i(1:3, istart:iend) = deriv_i(1:3, istart:iend) + self%temp_deriv_i(1:3, 1:self%angular_N)
                    deriv_j(1:3, istart:iend, j) = deriv_j(1:3, istart:iend, j) + self%temp_deriv_j(1:3, 1:self%angular_N)
                    deriv_j(1:3, istart:iend, k) = deriv_j(1:3, istart:iend, k) + self%temp_deriv_k(1:3, 1:self%angular_N)


                else
                    call self%evaluate_angular(R_ij, R_ik, R_jk,d_ij, d_ik, d_jk, self%temp_values)
                endif
                !write(*,*) j,k,values(istart:iend)
                values(istart:iend) = values(istart:iend) + self%temp_values(1:self%angular_N)
                
                !write(*,*) "val",self%temp_values(1:self%angular_N)

                istart = self%angular_2_initialindex
                iend = self%angular_2_finalindex

                if (self%multi) then
                    s_k = self%typespin(self%typeid(type1_loc(k)))
                    values(istart:iend) = values(istart:iend) + s_j*s_k*self%temp_values(1:self%angular_N)
                    if (do_deriv) then
                        deriv_i(1:3, istart:iend) = deriv_i(1:3, istart:iend) &
                                            + s_j*s_k*self%temp_deriv_i(1:3, 1:self%angular_N)
                        deriv_j(1:3, istart:iend, j) = deriv_j(1:3,istart:iend, j) &
                                            + s_j*s_k*self%temp_deriv_j(1:3, 1:self%angular_N)
                        deriv_j(1:3, istart:iend, k) = deriv_j(1:3,istart:iend, k) &
                                            + s_j*s_k*self%temp_deriv_k(1:3, 1:self%angular_N)      
                    endif
                endif

            end do for_k

        end do for_j
    end subroutine evaluate_values_of_fingerprints_spline




    subroutine evaluate_radial_spline(self, R_ij, d_ij, values, deriv_i, deriv_j)
        use bspline,only:bspline_basis_functions,&
                        bspline_basis_functions_deriv,d => bspline_order
        implicit none
        class(FP_SplineBasis), intent(in) :: self
        real(real64), dimension(3),             intent(in)    :: R_ij
        real(real64),                           intent(in)    :: d_ij
        real(real64), dimension(:),             intent(out)   :: values
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_i
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_j

        !real(real64)                     :: w_ij, dw_ij
        real(real64), dimension(self%radial_N) :: f, df
        integer::i,k

        !w_ij = cutoff_fc(d_ij, self%radial_Rc)
        !f = chebyshev_polynomial(d_ij, 0.0d0, self%radial_Rc, self%radial_points)

        
        if (present(deriv_i) .and. present(deriv_j)) then
            !dw_ij = cutoff_fc_d1(d_ij, self%radial_Rc)
            call bspline_basis_functions_deriv(values(1:self%radial_N),df(1:self%radial_N),&
                    d_ij,self%radial_knots,d)
            !df = chebyshev_polynomial_d1(d_ij, 0.0d0, self%radial_Rc, self%radial_points)
            !forall (i=1:self%radial_N)
            do concurrent (k=1:3,i=1:self%radial_N)
               deriv_i(k,i) = -R_ij(k)/d_ij*(df(i))
               !deriv_i(k,i) = -R_ij(k)/d_ij*(dw_ij*f(i) + w_ij*df(i))
            end do
            !end forall
            !values(1:self%radial_N) = w_ij*f(1:self%radial_N)
            deriv_j(1:3,1:self%radial_N) = -deriv_i(1:3,1:self%radial_N)
        else 
            call bspline_basis_functions(values(1:self%radial_N),&
                        d_ij,self%radial_knots,d)
            !values(1:self%radial_N) = w_ij*f(1:self%radial_N)
        end if


    end subroutine  evaluate_radial_spline

    subroutine evaluate_angular_spline(self, R_ij, R_ik, R_jk,d_ij, d_ik, d_jk, values, &
                                            deriv_i, deriv_j, deriv_k)
        use bspline,only:bspline_basis_functions,&
                    bspline_basis_functions_deriv,d => bspline_order
        implicit none
        class(FP_SplineBasis), intent(in) :: self
        real(real64), dimension(3),             intent(in)    :: R_ij, R_ik,R_jk
        real(real64),                           intent(in)    :: d_ij, d_ik
        real(real64),                           intent(in)    :: d_jk!cos_ijk
        real(real64), dimension(:),             intent(out)   :: values
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_i
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_j
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_k

        real(real64)                     :: w_ijk
        real(real64)                     :: fc_j, dfc_j, fc_k, dfc_k
        real(real64), dimension(self%angular_N) :: f, df
        real(real64)                     :: id_ij2, id_ik2, id_ij_ik
        real(real64), dimension(3)       :: dj_d_jk, dk_d_jk!dj_cos_ikj, dk_cos_ikj
        real(real64), dimension(3)       :: dj_w_ijk, dk_w_ijk,di_w_ijk
        integer::k,i
        integer::N

        N = self%angular_N

        fc_j = cutoff_fc(d_ij, self%angular_Rc)
        fc_k = cutoff_fc(d_ik, self%angular_Rc)
        w_ijk = fc_j*fc_k

        !f = chebyshev_polynomial(cos_ijk, -1.0d0, 1.0d0, self%angular_points)
        !values(1:self%angular_N) = w_ijk*f

        if (present(deriv_i) .and. present(deriv_j) .and. present(deriv_k)) then
            call bspline_basis_functions_deriv(f(1:N),df(1:N),d_jk,self%angular_knots,d)

            dfc_j = cutoff_fc_d1(d_ij, self%angular_Rc)
            dfc_k = cutoff_fc_d1(d_ik, self%angular_Rc)
            !df = chebyshev_polynomial_d1(cos_ijk, -1.0d0, 1.0d0, self%angular_points)

            !id_ij2 = 1.0d0/(d_ij*d_ij)
            !id_ik2 = 1.0d0/(d_ik*d_ik)
            !id_ij_ik = 1.0d0/(d_ij*d_ik)

            !dj_cos_ikj = -cos_ijk*R_ij*id_ij2 + R_ik*id_ij_ik
            !dk_cos_ikj = -cos_ijk*R_ik*id_ik2 + R_ij*id_ij_ik

            dj_d_jk = -R_jk/d_jk
            dk_d_jk = R_jk/d_jk

            dj_w_ijk = dfc_j*fc_k*R_ij/d_ij
            dk_w_ijk = fc_j*dfc_k*R_ik/d_ik
            di_w_ijk = -dj_w_ijk  - dk_w_ijk
            !dj_w_ijk = dfc_j*fc_k*R_ij/d_ij
            !dk_w_ijk = fc_j*dfc_k*R_ik/d_ik

            do concurrent(k=1:3,i=1:self%angular_N)
                deriv_j(k,i) = dj_w_ijk(k)*f(i) + w_ijk*df(i)*dj_d_jk(k)
                !deriv_j(k,i) = dj_w_ijk(k)*f(i) + w_ijk*df(i)*dj_cos_ikj(k)
                !deriv_k(k,i) = dk_w_ijk(k)*f(i) + w_ijk*df(i)*dk_cos_ikj(k)
                deriv_k(k,i) = dk_w_ijk(k)*f(i) + w_ijk*df(i)*dk_d_jk(k)
                deriv_i(k,i) = di_w_ijk(k)*f(i)
                !deriv_i(k,i) = -deriv_j(k,i) -deriv_k(k,i)  
            end do
        else
            call bspline_basis_functions(f(1:N),d_jk,self%angular_knots,d)
        end if
        values(1:N) = w_ijk*f

    end subroutine evaluate_angular_spline


    
end module splinebasis_fp