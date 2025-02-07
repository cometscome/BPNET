module chebyshevbasis_fp
    use iso_fortran_env
    use fingerprint_basis,only:FingerPrintBasis
    use chebyshev,only:chebyshev_polynomial,chebyshev_polynomial_d1
    use chebyshev,only:chebyshev_polynomial_vector,chebyshev_polynomial_d1_vector
    use cutoffmodule,only:cutoff_fc,cutoff_fc_d1

    implicit none
    real(real64), parameter, private :: PI     = 3.14159265358979d0
    real(real64), parameter, private :: PI_INV = 1.0d0/PI
    real(real64), parameter, private :: PI2    = 2.0d0*PI
    real(real64), parameter, private :: EPS    = 1.0d-12


    type, extends(FingerPrintBasis) :: FP_ChebyshevBasis
        integer                                       :: radial_order
        real(real64)                                       :: radial_Rc
        integer                                       :: radial_N
        integer                                       :: angular_order
        real(real64)                                       :: angular_Rc
        integer                                       :: angular_N
        integer                                       :: num_of_different_species

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
        logical::newchebyshev

        contains
        !procedure :: area => circle_area
        !procedure :: test => chebyshev_test
        procedure ::  get_num_parameters =>  get_num_parameters_chebyshev
        procedure :: evaluate => evaluate_values_of_fingerprints_chebyshev
        procedure :: evaluate_radial => evaluate_radial_chebyshev
        procedure :: evaluate_angular => evaluate_angular_chebyshev
        procedure :: evaluate_angular_new => evaluate_angular_chebyshev_new
        procedure :: evaluate_angular_modified => evaluate_angular_chebyshev_modified
        !procedure :: evaluate => evaluate_values_of_fingerprints_chebyshev_modified
    end type FP_ChebyshevBasis

    interface FP_ChebyshevBasis
        procedure::setup_basis_chebyshev
    end interface

    contains

    !function chebyshev_test(this) result(res)
    !    class(FP_ChebyshevBasis), intent(in):: this
    !    real :: res
    !end function chebyshev_test

    function get_num_parameters_chebyshev(self) result(num)
        class(FP_ChebyshevBasis), intent(in):: self
        integer :: num
        num = self%num_of_parameters
    end function get_num_parameters_chebyshev

    subroutine evaluate_values_of_fingerprints_chebyshev(self,itype0,cartesian_coordinates, num_of_neighbor_atoms, &
                                                        coordinates_of_neighbor_atoms, types_of_neighbor_atoms, ltype,&
                                                        values,&
                                                        deriv_i,deriv_j,scaled)
        implicit none
        class(FP_ChebyshevBasis), intent(inout) :: self
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
        real(real64)::R_ij(1:3),R_ik(1:3)
        real(real64)::d_ij,d_ik,cos_ijk
        integer::istart,iend
        real(real64)::RjRk

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
                !write(*,*) R_ij, R_ik
                !write(*,*) "d",d_ij, d_ik, cos_ijk
                RjRk = dot_product(R_ij, R_ik)
                cos_ijk =  RjRk/(d_ij*d_ik)

                istart = self%angular_1_initialindex
                iend = self%angular_1_finalindex

                if (do_deriv) then
                    
                    if (self%newchebyshev) then
                        call self%evaluate_angular_new(R_ij, R_ik, d_ij, d_ik, RjRk, self%temp_values, &
                                    deriv_i=self%temp_deriv_i, deriv_j=self%temp_deriv_j,deriv_k=self%temp_deriv_k)
                    else
                        call self%evaluate_angular(R_ij, R_ik, d_ij, d_ik, cos_ijk, self%temp_values, &
                                    deriv_i=self%temp_deriv_i, deriv_j=self%temp_deriv_j,deriv_k=self%temp_deriv_k)
                    endif
                    !write(*,*) j,k,self%temp_values(1:self%angular_N)

                    deriv_i(1:3, istart:iend) = deriv_i(1:3, istart:iend) + self%temp_deriv_i(1:3, 1:self%angular_N)
                    deriv_j(1:3, istart:iend, j) = deriv_j(1:3, istart:iend, j) + self%temp_deriv_j(1:3, 1:self%angular_N)
                    deriv_j(1:3, istart:iend, k) = deriv_j(1:3, istart:iend, k) + self%temp_deriv_k(1:3, 1:self%angular_N)


                else
                    call self%evaluate_angular(R_ij, R_ik, d_ij, d_ik, cos_ijk, self%temp_values)
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
    end subroutine evaluate_values_of_fingerprints_chebyshev

    subroutine evaluate_values_of_fingerprints_chebyshev_modified(self,itype0,cartesian_coordinates, num_of_neighbor_atoms, &
                                                        coordinates_of_neighbor_atoms, types_of_neighbor_atoms, ltype,&
                                                        values,&
                                                        deriv_i,deriv_j,scaled)
        implicit none
        class(FP_ChebyshevBasis), intent(inout) :: self
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
        real(real64)::R_ij(1:3),R_ik(1:3)
        real(real64)::d_ij,d_ik,cos_ijk
        integer::istart,iend
        real(real64)::RjRk

        integer               :: type0_loc
        integer, dimension(num_of_neighbor_atoms) :: type1_loc

        real(real64)::R_ik_set(1:3,num_of_neighbor_atoms)
        real(real64)::d_ik_set(num_of_neighbor_atoms)
        real(real64)::cos_ijk_set(num_of_neighbor_atoms)
        real(real64)::fc_k_set(num_of_neighbor_atoms)
        real(real64)::s_k_set(num_of_neighbor_atoms)
        integer::istart1,istart2
        

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
            do concurrent(k=j+1:num_of_neighbor_atoms)
            !do k=j+1,num_of_neighbor_atoms
                R_ik_set(1:3,k) = coordinates_of_neighbor_atoms(1:3,k) - cartesian_coordinates(1:3)
                d_ik_set(k) = sqrt(dot_product(R_ik, R_ik_set(1:3,k)))
                if ((d_ik_set(k) > self%angular_Rc) .or. (d_ik_set(k) < EPS)) cycle
                RjRk = dot_product(R_ij, R_ik)
                cos_ijk_set(k) =  RjRk/(d_ij*d_ik_set(k))
                fc_k_set(k) = cutoff_fc(d_ik_set(k),self%angular_Rc)

                if (self%multi) then
                    s_k_set(k) = self%typespin(self%typeid(type1_loc(k)))
                end if
            end do

            istart = self%angular_1_initialindex
            iend = self%angular_1_finalindex

            istart1 = self%angular_1_initialindex
            istart2 = self%angular_2_initialindex
        
            if (do_deriv) then
                s_j = self%typespin(self%typeid(type1_loc(j)))
                call self%evaluate_angular_modified(R_ij, R_ik_set, d_ij, d_ik_set, cos_ijk_set, istart1,istart2,&
                                    s_j,s_k_set,fc_k_set,j,&
                                    values, self%multi,&
                                    deriv_i=deriv_i, deriv_j=deriv_j)!,deriv_k=self%temp_deriv_k)
        !R_ij, R_ik_set, d_ij, d_ik_set, cos_ijk_set, istart1,istart2,&
                                    !s_j,s_k_set,fc_k_set,jatom,values, multi
                !deriv_i(1:3, istart:iend) = deriv_i(1:3, istart:iend) + self%temp_deriv_i(1:3, 1:self%angular_N)
                !deriv_j(1:3, istart:iend, j) = deriv_j(1:3, istart:iend, j) + self%temp_deriv_j(1:3, 1:self%angular_N)
                !deriv_j(1:3, istart:iend, k) = deriv_j(1:3, istart:iend, k) + self%temp_deriv_k(1:3, 1:self%angular_N)

            else
                call self%evaluate_angular_modified(R_ij, R_ik_set, d_ij, d_ik_set,cos_ijk_set,&
                                    istart1,istart2, s_j,s_k_set,fc_k_set,j,&
                                    values,self%multi)
            end if
            !values(istart:iend) = values(istart:iend) + self%temp_values(1:self%angular_N)


        end do for_j
    end subroutine evaluate_values_of_fingerprints_chebyshev_modified



    subroutine evaluate_radial_chebyshev(self, R_ij, d_ij, values, deriv_i, deriv_j)
        implicit none
        class(FP_ChebyshevBasis), intent(in) :: self
        real(real64), dimension(3),             intent(in)    :: R_ij
        real(real64),                           intent(in)    :: d_ij
        real(real64), dimension(:),             intent(out)   :: values
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_i
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_j

        real(real64)                     :: w_ij, dw_ij
        real(real64), dimension(self%radial_N) :: f, df
        integer::i,k

        w_ij = cutoff_fc(d_ij, self%radial_Rc)
        f = chebyshev_polynomial(d_ij, 0.0d0, self%radial_Rc, self%radial_order)

        values(1:self%radial_N) = w_ij*f(1:self%radial_N)
        if (present(deriv_i) .and. present(deriv_j)) then
            dw_ij = cutoff_fc_d1(d_ij, self%radial_Rc)
            df = chebyshev_polynomial_d1(d_ij, 0.0d0, self%radial_Rc, self%radial_order)
            do concurrent (k=1:3,i=1:self%radial_N)
            !forall (i=1:self%radial_N)
               deriv_i(k,i) = -R_ij(k)/d_ij*(dw_ij*f(i) + w_ij*df(i))
            !end forall
            end do
            deriv_j(1:3,1:self%radial_N) = -deriv_i(1:3,1:self%radial_N)
         end if


    end subroutine  evaluate_radial_chebyshev

    subroutine evaluate_angular_chebyshev(self, R_ij, R_ik, d_ij, d_ik, cos_ijk, values, &
                                            deriv_i, deriv_j, deriv_k)
        implicit none
        class(FP_ChebyshevBasis), intent(in) :: self
        real(real64), dimension(3),             intent(in)    :: R_ij, R_ik
        real(real64),                           intent(in)    :: d_ij, d_ik
        real(real64),                           intent(in)    :: cos_ijk
        real(real64), dimension(:),             intent(out)   :: values
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_i
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_j
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_k

        real(real64)                     :: w_ijk
        real(real64)                     :: fc_j, dfc_j, fc_k, dfc_k
        real(real64), dimension(self%angular_N) :: f, df
        real(real64)                     :: id_ij2, id_ik2, id_ij_ik
        real(real64), dimension(3)       :: dj_cos_ikj, dk_cos_ikj
        real(real64), dimension(3)       :: dj_w_ijk, dk_w_ijk
        integer::k,i

        fc_j = cutoff_fc(d_ij, self%angular_Rc)
        fc_k = cutoff_fc(d_ik, self%angular_Rc)
        w_ijk = fc_j*fc_k

        f = chebyshev_polynomial(cos_ijk, -1.0d0, 1.0d0, self%angular_order)
        values(1:self%angular_N) = w_ijk*f

        if (present(deriv_i) .and. present(deriv_j) .and. present(deriv_k)) then
            dfc_j = cutoff_fc_d1(d_ij, self%angular_Rc)
            dfc_k = cutoff_fc_d1(d_ik, self%angular_Rc)
            df = chebyshev_polynomial_d1(cos_ijk, -1.0d0, 1.0d0, self%angular_order)

            id_ij2 = 1.0d0/(d_ij*d_ij)
            id_ik2 = 1.0d0/(d_ik*d_ik)
            id_ij_ik = 1.0d0/(d_ij*d_ik)

            dj_cos_ikj = -cos_ijk*R_ij*id_ij2 + R_ik*id_ij_ik
            dk_cos_ikj = -cos_ijk*R_ik*id_ik2 + R_ij*id_ij_ik

            dj_w_ijk = dfc_j*fc_k*R_ij/d_ij
            dk_w_ijk = fc_j*dfc_k*R_ik/d_ik

            do concurrent(k=1:3,i=1:self%angular_N)
                deriv_j(k,i) = dj_w_ijk(k)*f(i) + w_ijk*df(i)*dj_cos_ikj(k)
                deriv_k(k,i) = dk_w_ijk(k)*f(i) + w_ijk*df(i)*dk_cos_ikj(k)
                deriv_i(k,i) = -deriv_j(k,i) -deriv_k(k,i)  
            end do
        end if

    end subroutine evaluate_angular_chebyshev

    subroutine evaluate_angular_chebyshev_modified(self, R_ij, R_ik_set, d_ij, d_ik_set, cos_ijk_set, istart1,istart2,&
                                s_j,s_k_set,fc_k_set,jatom,values, multi,&
                                            deriv_i, deriv_j)
        implicit none
        class(FP_ChebyshevBasis), intent(in) :: self
        real(real64), dimension(3),             intent(in)    :: R_ij
        real(real64),                           intent(in)    :: d_ij, d_ik_set(:)
        real(real64),                           intent(in)    :: fc_k_set(:),s_j
        real(real64),                           intent(in)    :: cos_ijk_set(:),R_ik_set(:,:),s_k_set(:)
        integer,intent(in) ::istart1,istart2,jatom
        real(real64), dimension(:),             intent(out)   :: values
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_i
        real(real64), dimension(:,:,:), optional, intent(out)   :: deriv_j
        logical,intent(in)::multi

        real(real64)                     :: w_ijk,fc_j
        real(real64)                     :: dfc_j, fc_k, dfc_k
        real(real64), dimension(self%angular_N,size(d_ik_set)) :: f, df
        real(real64)                     :: id_ij2, id_ik2, id_ij_ik
        real(real64), dimension(3)       :: dj_cos_ikj, dk_cos_ikj
        real(real64), dimension(3)       :: dj_w_ijk, dk_w_ijk
        real(real64)::t1,t2
        integer::k,i,katom

        fc_j = cutoff_fc(d_ij, self%angular_Rc)
        f = chebyshev_polynomial_vector(cos_ijk_set, -1.0d0, 1.0d0, self%angular_order)

        if (present(deriv_i) .and. present(deriv_j) ) then
            df = chebyshev_polynomial_d1_vector(cos_ijk_set, -1.0d0, 1.0d0, self%angular_order)
        end if

        do concurrent (katom=1:size(d_ik_set))
            fc_k = cutoff_fc(d_ik_set(katom), self%angular_Rc)
            w_ijk = fc_j*fc_k_set(katom)
            values(istart1:istart1 -1+ self%angular_N) = values(1:self%angular_N) + w_ijk*f(:,katom)
            if (multi) then
                values(istart2:istart2 -1+ self%angular_N) = values(1:self%angular_N) + s_j*s_k_set(katom)*w_ijk*f(:,katom)
            end if
        
            if (present(deriv_i) .and. present(deriv_j) ) then
                dfc_j = cutoff_fc_d1(d_ij, self%angular_Rc)
                dfc_k = cutoff_fc_d1(d_ik_set(katom), self%angular_Rc)
                !df = chebyshev_polynomial_d1(cos_ijk, -1.0d0, 1.0d0, self%angular_order)

                id_ij2 = 1.0d0/(d_ij*d_ij)
                id_ik2 = 1.0d0/(d_ik_set(katom)*d_ik_set(katom))
                id_ij_ik = 1.0d0/(d_ij*d_ik_set(katom))

                dj_cos_ikj = -cos_ijk_set(katom)*R_ij*id_ij2 + R_ik_set(:,katom)*id_ij_ik
                dk_cos_ikj = -cos_ijk_set(katom)*R_ik_set(:,katom)*id_ik2 + R_ij*id_ij_ik

                dj_w_ijk = dfc_j*fc_k_set(katom)*R_ij/d_ij
                dk_w_ijk = fc_j*dfc_k*R_ik_set(:,katom)/d_ik_set(katom)

                do concurrent(k=1:3,i=1:self%angular_N)
                    t1 = dj_w_ijk(k)*f(i,katom) + w_ijk*df(i,katom)*dj_cos_ikj(k)
                    t2 = dk_w_ijk(k)*f(i,katom) + w_ijk*df(i,katom)*dk_cos_ikj(k)

                    !deriv_j(k,i) = deriv_j(k,i) + dj_w_ijk(k)*f(i,katom) + w_ijk*df(i,katom)*dj_cos_ikj(k)
                    !deriv_k(k,i) = deriv_k(k,i) + dk_w_ijk(k)*f(i,katom) + w_ijk*df(i,katom)*dk_cos_ikj(k)
                    !deriv_i(k,i) = deriv_i(k,i) + -deriv_j(k,i) -deriv_k(k,i)  


                    deriv_j(k, istart1+i-1, jatom) = deriv_j(k, istart1+i-1, jatom) + t1
                    deriv_j(k, istart1+i-1, katom) = deriv_j(k, istart1+i-1, katom) + t2
                    deriv_i(k, istart1+i-1) = deriv_i(k, istart1+i-1) -t1-t2

                    if (multi) then
                        deriv_j(k, istart2+i-1, jatom) = deriv_j(k, istart2+i-1, jatom) + s_j*s_k_set(katom)*t1
                        deriv_j(k, istart2+i-1, katom) = deriv_j(k, istart2+i-1, katom) + s_j*s_k_set(katom)*t2
                        deriv_i(k, istart2+i-1) = deriv_i(k, istart2+i-1) -s_j*s_k_set(katom)*t1-s_j*s_k_set(katom)*t2
                    end if

                end do
            end if
        end do

    end subroutine evaluate_angular_chebyshev_modified    

    subroutine evaluate_angular_chebyshev_new(self, R_ij, R_ik, d_ij, d_ik, RjRk, values, &
                                            deriv_i, deriv_j, deriv_k)
        implicit none
        class(FP_ChebyshevBasis), intent(in) :: self
        real(real64), dimension(3),             intent(in)    :: R_ij, R_ik
        real(real64),                           intent(in)    :: d_ij, d_ik
        real(real64),                           intent(in)    :: RjRk
        real(real64), dimension(:),             intent(out)   :: values
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_i
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_j
        real(real64), dimension(:,:), optional, intent(out)   :: deriv_k

        real(real64)                     :: w_ijk
        real(real64)                     :: fc_j, dfc_j, fc_k, dfc_k
        real(real64), dimension(self%angular_N) :: f, df
        real(real64)                     :: id_ij2, id_ik2, id_ij_ik
        real(real64), dimension(3)       :: dj_cos_ikj, dk_cos_ikj
        real(real64), dimension(3)       :: dj_w_ijk, dk_w_ijk
        integer::k,i

        fc_j = cutoff_fc(d_ij, self%angular_Rc)
        fc_k = cutoff_fc(d_ik, self%angular_Rc)
        w_ijk = fc_j*fc_k

        !f = chebyshev_polynomial(cos_ijk, -1.0d0, 1.0d0, self%angular_order)
        f = chebyshev_polynomial(RjRk, -self%angular_Rc, self%angular_Rc, self%angular_order)
        values(1:self%angular_N) = w_ijk*f

        if (present(deriv_i) .and. present(deriv_j) .and. present(deriv_k)) then
            do concurrent(k=1:3,i=1:self%angular_N)
                deriv_j(k,i) = 0d0
                deriv_k(k,i) = 0d0
                deriv_i(k,i) = 0d0
            end do
        end if

    end subroutine evaluate_angular_chebyshev_new

    type(FP_ChebyshevBasis) function setup_basis_chebyshev(fingerprint_parameters,num_of_different_species,&
                    species_of_surrounding_atom    ) result(fpbasis)
        implicit none
        real(real64),dimension(:,:),intent(in)::fingerprint_parameters
        integer,intent(in)::num_of_different_species
        character(len=2), dimension(:),intent(in)::species_of_surrounding_atom
        integer::i,s

        !fpbasis%newchebyshev = .true.
        fpbasis%newchebyshev = .false.


        fpbasis%radial_Rc = fingerprint_parameters(1,1)
        fpbasis%radial_order =  nint(fingerprint_parameters(2,1))
        fpbasis%angular_Rc =  fingerprint_parameters(3,1)
        fpbasis%angular_order = nint(fingerprint_parameters(4,1))
        fpbasis%num_of_different_species = num_of_different_species

        fpbasis%radial_N = fpbasis%radial_order + 1
        fpbasis%angular_N = fpbasis%angular_order + 1
        !fpbasis%num_values = max(fpbasis%r_N, fpbasis%angular_N)

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


    end function setup_basis_chebyshev
    
end module chebyshevbasis_fp