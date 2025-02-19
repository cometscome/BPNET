module shared_potentials
    use iso_fortran_env
    use mod_bppotential,only:BPpotential
    use bfio,only:FILEPATHLEN
    use lclist,   only: lcl_init,          &
                        lcl_final,         &
                        lcl_print_info,    &
                        lcl_nmax_nbdist,   &
                        lcl_nbdist_cart
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_char, c_bool, &
                        c_ptr, c_f_pointer, c_null_char
                        
    implicit none
    type(BPpotential), allocatable :: BPpots(:)
    real(real64)::Rc_min,Rc_max
    integer::num_of_maxatoms_in_sphere
    integer::num_max_coefficients
    real(real64), dimension(:),     allocatable :: fingerprint_values
    real(real64), dimension(:,:),   allocatable :: fingerprint_deriv_i
    real(real64), dimension(:,:,:), allocatable :: fingerprint_deriv_j
    integer, parameter :: TYPENAME_LEN = 2
    


    public::initialize_potentials
    public::bpnet_atomic_energy_and_forces_novirial,set_Rcmax_and_Rcmin
    public::deallocate_potentials,reload_potentials


     

    contains 

    subroutine deallocate_potentials()
        implicit none
        integer::ntypes
        integer::i

        ntypes = size(BPpots,1)
        do i=1,ntypes
            call BPpots(i)%deallocate()
        end do

    end subroutine

    subroutine reload_potentials(num_of_Types, typeName, network_filenames)
        integer, intent(in)                                   :: num_of_Types
        character(len=2), dimension(:), intent(in)           :: typeName
        character(len=FILEPATHLEN), dimension(:), intent(in) :: network_filenames

        integer :: i
        !character(kind=c_char, len=2) :: typeName_c(num_of_Types)
        !character(kind=c_char, len=FILEPATHLEN) :: network_filenames_c(num_of_Types)
        character(kind=c_char, len=TYPENAME_LEN*num_of_Types) :: typeName_c
        character(kind=c_char, len=FILEPATHLEN*num_of_Types)  :: network_filenames_c


        do i = 1, num_of_Types
            typeName_c((i-1)*TYPENAME_LEN+1:i*TYPENAME_LEN) = typeName(i)
            network_filenames_c((i-1)*FILEPATHLEN+1:i*FILEPATHLEN) = network_filenames(i)
        end do


        call reload_potentials_c(num_of_Types, typeName_c, network_filenames_c)
    end subroutine reload_potentials

    subroutine reload_potentials_c(num_of_Types, typeName, network_filenames) bind(C, name="reload_potentials")
        integer(c_int), intent(in)                            :: num_of_Types
        character(kind=c_char), intent(in)                   :: typeName(*)
        character(kind=c_char), intent(in)                   :: network_filenames(*)
        character(len=TYPENAME_LEN) :: typeName_f
        character(len=FILEPATHLEN) :: network_filenames_f
        integer::itype,i

        call deallocate_potentials()

        do i=1,num_of_Types
            call copy_c_string_to_fortran(typeName((i-1)*TYPENAME_LEN+1:i*TYPENAME_LEN), TYPENAME_LEN, typeName_f)
            call copy_c_string_to_fortran(network_filenames((i-1)*FILEPATHLEN+1:i*FILEPATHLEN), FILEPATHLEN, network_filenames_f)

            BPpots(i) = BPpotential(typeName_f,network_filenames_f)
            !BPpots(i) = BPpotential(typeName(i),network_filenames(i))
        end do

        call set_Rcmax_and_Rcmin(num_of_types)

        num_max_coefficients  = 0
        do itype = 1, num_of_types
            num_max_coefficients = max(num_max_coefficients,BPpots(itype)%nsf)
        end do

        num_of_maxatoms_in_sphere = lcl_nmax_nbdist(Rc_min, Rc_max)

        if (allocated(fingerprint_values)) deallocate(fingerprint_values)
        if (allocated(fingerprint_deriv_i)) deallocate(fingerprint_deriv_i)
        if (allocated(fingerprint_deriv_j)) deallocate(fingerprint_deriv_j)

        allocate(fingerprint_values(num_max_coefficients), source=0d0)
        allocate(fingerprint_deriv_i(3,num_max_coefficients), source=0d0)
        allocate(fingerprint_deriv_j(3,num_max_coefficients,num_of_maxatoms_in_sphere), source=0d0)

    end subroutine

    subroutine reload_potentials_old(num_of_Types,typeName,network_filenames)
        implicit none
        integer,intent(in)                                             :: num_of_Types
        character(len=2), dimension(:),intent(in) :: typeName
        character(len=FILEPATHLEN), dimension(:),intent(in) :: network_filenames
        integer::itype,i

        call deallocate_potentials()

        do i=1,num_of_Types
            BPpots(i) = BPpotential(typeName(i),network_filenames(i))
        end do

        call set_Rcmax_and_Rcmin(num_of_types)

        num_max_coefficients  = 0
        do itype = 1, num_of_types
            num_max_coefficients = max(num_max_coefficients,BPpots(itype)%nsf)
        end do

        num_of_maxatoms_in_sphere = lcl_nmax_nbdist(Rc_min, Rc_max)

        if (allocated(fingerprint_values)) deallocate(fingerprint_values)
        if (allocated(fingerprint_deriv_i)) deallocate(fingerprint_deriv_i)
        if (allocated(fingerprint_deriv_j)) deallocate(fingerprint_deriv_j)

        allocate(fingerprint_values(num_max_coefficients), source=0d0)
        allocate(fingerprint_deriv_i(3,num_max_coefficients), source=0d0)
        allocate(fingerprint_deriv_j(3,num_max_coefficients,num_of_maxatoms_in_sphere), source=0d0)

    end subroutine

    subroutine initialize_potentials(num_of_Types, typeName, network_filenames)
        integer, intent(in)                                   :: num_of_Types
        character(len=2), dimension(:), intent(in)           :: typeName
        character(len=FILEPATHLEN), dimension(:), intent(in) :: network_filenames

        integer :: i
        !character(kind=c_char, len=2) :: typeName_c(num_of_Types)
        !character(kind=c_char, len=FILEPATHLEN) :: network_filenames_c(num_of_Types)
        character(kind=c_char, len=TYPENAME_LEN*num_of_Types) :: typeName_c
        character(kind=c_char, len=FILEPATHLEN*num_of_Types)  :: network_filenames_c


        do i = 1, num_of_Types
            typeName_c((i-1)*TYPENAME_LEN+1:i*TYPENAME_LEN) = typeName(i)
            network_filenames_c((i-1)*FILEPATHLEN+1:i*FILEPATHLEN) = network_filenames(i)
        end do


        call initialize_potentials_c(num_of_Types, typeName_c, network_filenames_c)
    end subroutine initialize_potentials

    subroutine copy_c_string_to_fortran(c_string, c_len, f_string)
        character(kind=c_char), intent(in) :: c_string(*)
        integer, intent(in)                :: c_len
        character(len=*), intent(out)      :: f_string
        integer :: i

        f_string = ' '  ! 初期化
        do i = 1, min(c_len, len(f_string))
            f_string(i:i) = c_string(i)
        end do
    end subroutine copy_c_string_to_fortran

    subroutine initialize_potentials_c(num_of_Types, typeName, network_filenames) bind(C, name="initialize_potentials")
        integer(c_int), intent(in)                            :: num_of_Types
        character(kind=c_char), intent(in)                   :: typeName(*)
        character(kind=c_char), intent(in)                   :: network_filenames(*)
        character(len=TYPENAME_LEN) :: typeName_f
        character(len=FILEPATHLEN) :: network_filenames_f
        
        integer :: i,itype

        allocate(BPpots(1:num_of_Types))


        do i=1,num_of_Types
            call copy_c_string_to_fortran(typeName((i-1)*TYPENAME_LEN+1:i*TYPENAME_LEN), TYPENAME_LEN, typeName_f)
            call copy_c_string_to_fortran(network_filenames((i-1)*FILEPATHLEN+1:i*FILEPATHLEN), FILEPATHLEN, network_filenames_f)

            BPpots(i) = BPpotential(typeName_f,network_filenames_f)
            !BPpots(i) = BPpotential(typeName(i),network_filenames(i))
        end do

        call set_Rcmax_and_Rcmin(num_of_types)

        num_max_coefficients  = 0
        do itype = 1, num_of_types
            num_max_coefficients = max(num_max_coefficients,BPpots(itype)%nsf)
        end do

        num_of_maxatoms_in_sphere = lcl_nmax_nbdist(Rc_min, Rc_max)

        allocate(fingerprint_values(num_max_coefficients), source=0d0)
        allocate(fingerprint_deriv_i(3,num_max_coefficients), source=0d0)
        allocate(fingerprint_deriv_j(3,num_max_coefficients,num_of_maxatoms_in_sphere), source=0d0)


    end subroutine initialize_potentials_c

    subroutine initialize_potentials_old(num_of_Types,typeName,network_filenames) !,Rc_min, Rc_max)
        implicit none
        integer,intent(in)                                             :: num_of_Types
        character(len=2), dimension(:),intent(in) :: typeName
        character(len=FILEPATHLEN), dimension(:),intent(in) :: network_filenames
        integer::itype,i
        !real(real64)::Rc_min,Rc_max
        

        allocate(BPpots(1:num_of_Types))

        do i=1,num_of_Types
            BPpots(i) = BPpotential(typeName(i),network_filenames(i))
        end do

        call set_Rcmax_and_Rcmin(num_of_types)

        num_max_coefficients  = 0
        do itype = 1, num_of_types
            num_max_coefficients = max(num_max_coefficients,BPpots(itype)%nsf)
        end do

        num_of_maxatoms_in_sphere = lcl_nmax_nbdist(Rc_min, Rc_max)

        allocate(fingerprint_values(num_max_coefficients), source=0d0)
        allocate(fingerprint_deriv_i(3,num_max_coefficients), source=0d0)
        allocate(fingerprint_deriv_j(3,num_max_coefficients,num_of_maxatoms_in_sphere), source=0d0)



    end subroutine


    subroutine set_Rcmax_and_Rcmin(num_of_types)
        integer,intent(in)::num_of_types
        integer::itype
        

        Rc_min =BPpots(1)%Rc_min
        Rc_max =BPpots(1)%Rc_max

        do itype = 1, num_of_types
            Rc_min = min(Rc_min, BPpots(itype)%Rc_min)
            Rc_max = max(Rc_max, BPpots(itype)%Rc_max)
        end do
            
    end subroutine set_Rcmax_and_Rcmin

    subroutine bpnet_atomic_energy_and_forces_novirial( &
        coo_i, type_i, index_i, n_j, coo_j, type_j, index_j, natoms, &
        E_i, F, stat) bind(C)

        use, intrinsic :: iso_c_binding, only: c_double, c_int, c_char, c_bool, &
        c_ptr, c_f_pointer, c_null_char
        implicit none
    
        real(kind=c_double), dimension(3),        intent(in)    :: coo_i
        integer(kind=c_int), value,               intent(in)    :: type_i
        integer(kind=c_int), value,               intent(in)    :: index_i
        integer(kind=c_int), value,               intent(in)    :: n_j
        real(kind=c_double), dimension(3,n_j),    intent(in)    :: coo_j
        integer(kind=c_int), dimension(n_j),      intent(in)    :: type_j
        integer(kind=c_int), dimension(n_j),      intent(in)    :: index_j
        integer(kind=c_int), value,               intent(in)    :: natoms
        real(kind=c_double),                      intent(out)   :: E_i
        real(kind=c_double), dimension(3,natoms), intent(inout) :: F
        integer(kind=c_int),                      intent(out)   :: stat

        
        real(real64)::E_i_arr(1:1)
        real(real64)::dEdE(1:1)
        real(real64),allocatable::dEdc(:)
        integer::j
        real(real64),allocatable::Ftest(:,:)

        INTEGER :: start_clock, end_clock, clock_rate
        REAL(real64) :: elapsed_time
        integer::it


        !write(*,*) "energy"
        !write(*,*) index_i,&
        !coo_i, &
        !n_j, &
        !coo_j, &
        !type_j, &
        !BPpots(index_i)%ltype
        dEdE = 1d0
        allocate(dEdc(1:BPpots(type_i)%nsf))

        
        call BPpots(type_i)%fingerprint_basis%evaluate(type_i,&
            coo_i, &
            n_j, &
            coo_j, &
            type_j, &
            BPpots(type_i)%ltype,&
            values=fingerprint_values(1:BPpots(type_i)%nsf), &
            deriv_i=fingerprint_deriv_i(:,1:BPpots(type_i)%nsf), &
            deriv_j=fingerprint_deriv_j(:,1:BPpots(type_i)%nsf,:))
        


        call BPpots(type_i)%normalize(fingerprint_values(1:BPpots(type_i)%nsf),&
                fingerprint_deriv_i(:,1:BPpots(type_i)%nsf), &
                fingerprint_deriv_j(:,1:BPpots(type_i)%nsf,:),&
                n_j)

        !write(*,*)  "d",sum(fingerprint_deriv_j(:,1:BPpots(type_i)%nsf,:))

        !stop
        
        !write(*,*) "fingerprint_values",fingerprint_values(1:BPpots(type_i)%nsf)
        !stop

        call BPpots(type_i)%forward(fingerprint_values(1:BPpots(type_i)%nsf),E_i_arr)


    
        



        call BPpots(type_i)%backward(dEdE,dEdc)


        !write(*,*) E_i_arr
        !stop

        E_i = BPpots(type_i)%E_scale*E_i_arr(1) + BPpots(type_i)%shift
        !write(*,*) E_i,BPpots(type_i)%E_scale,BPpots(type_i)%shift
        E_i = E_i + BPpots(type_i)%E_atom(type_i)
        !write(*,*) E_i
        !stop
        !write(*,*) dEdc
        !call stp_eval(type_i, coo_i, n_j, coo_j, type_j, &
        !          aenet_pot(type_i)%stp, sfval=sfval, &
        !          sfderiv_i=sfderiv_i, sfderiv_j=sfderiv_j, scaled=.true.)
        F(1:3, index_i) = F(1:3, index_i) -  BPpots(type_i)%E_scale &
                    * matmul(fingerprint_deriv_i(1:3,1:BPpots(type_i)%nsf), dEdc(1:BPpots(type_i)%nsf))
        !write(*,*) F(1:3, index_i)
        

        !allocate(Ftest(3,n_j))
        do j = 1, n_j
            !Ftest(1:3,j) = -  BPpots(index_i)%scale &
            !        * matmul(fingerprint_deriv_j(1:3,1:BPpots(index_i)%nsf,j), dEdc(1:BPpots(index_i)%nsf))

            F(1:3, index_j(j)) = F(1:3, index_j(j)) -  BPpots(type_i)%E_scale &
                            * matmul(fingerprint_deriv_j(1:3,1:BPpots(type_i)%nsf,j), dEdc(1:BPpots(type_i)%nsf))
        end do

        !stop



        !call check_force()
        !stop

        contains 

        subroutine check_force()
            implicit none
            real(real64)::value1(1:BPpots(type_i)%nsf),value2(1:BPpots(type_i)%nsf)
            real(real64),allocatable:: coo_i2(:),coo_j2(:,:)
            real(real64)::eta,E_i_arr2(1:1),E_i2
            integer::ii

            allocate(coo_i2,source= coo_i)
            allocate(coo_j2,source= coo_j)
            eta = 1d-6

            do ii=1,3
                coo_i2 = coo_i
                coo_i2(ii) = coo_i2(ii) +eta
                call BPpots(type_i)%fingerprint_basis%evaluate(type_i,&
                    coo_i2, &
                    n_j, &
                    coo_j, &
                    type_j, &
                    BPpots(type_i)%ltype,&
                    values=value1,&
                    deriv_i=fingerprint_deriv_i(:,1:BPpots(type_i)%nsf), &
                    deriv_j=fingerprint_deriv_j(:,1:BPpots(type_i)%nsf,:))
                
                call BPpots(type_i)%forward(value1,E_i_arr2)
                E_i2 = BPpots(type_i)%E_scale*E_i_arr2(1) + BPpots(type_i)%shift
                E_i2 = E_i2 + BPpots(type_i)%E_atom(type_i)
                write(*,*) (E_i-E_i2)/eta,F(ii,index_i)
            end do
            coo_i2 = coo_i
            do j = 1, n_j
                do ii=1,3
                    coo_j2 = coo_j
                    coo_j2(ii,j) = coo_j2(ii,j) + eta
                    call BPpots(type_i)%fingerprint_basis%evaluate(type_i,&
                        coo_i2, &
                        n_j, &
                        coo_j2, &
                        type_j, &
                        BPpots(type_i)%ltype,&
                        values=value1,&
                        deriv_i=fingerprint_deriv_i(:,1:BPpots(type_i)%nsf), &
                        deriv_j=fingerprint_deriv_j(:,1:BPpots(type_i)%nsf,:))
                
                    call BPpots(type_i)%forward(value1,E_i_arr2)
                    E_i2 = BPpots(type_i)%E_scale*E_i_arr2(1) + BPpots(type_i)%shift
                    E_i2 = E_i2 + BPpots(type_i)%E_atom(index_i)
                    write(*,*) E_i,E_i2, E_i-E_i2
                    write(*,*)j,ii,(E_i2-E_i)/eta,Ftest(ii,j),(E_i-E_i2)/eta-Ftest(ii,j)!,F(ii,index_j(j))
                end do
            end do



            
        end subroutine

    end subroutine


    subroutine bpnet_atomic_energy_and_forces( &
        coo_i, type_i, index_i, n_j, coo_j, type_j, index_j, natoms, &
        E_i, F, S,stat) bind(C)

        use, intrinsic :: iso_c_binding, only: c_double, c_int, c_char, c_bool, &
        c_ptr, c_f_pointer, c_null_char
        implicit none
    
        real(kind=c_double), dimension(3),        intent(in)    :: coo_i
        integer(kind=c_int), value,               intent(in)    :: type_i
        integer(kind=c_int), value,               intent(in)    :: index_i
        integer(kind=c_int), value,               intent(in)    :: n_j
        real(kind=c_double), dimension(3,n_j),    intent(in)    :: coo_j
        integer(kind=c_int), dimension(n_j),      intent(in)    :: type_j
        integer(kind=c_int), dimension(n_j),      intent(in)    :: index_j
        integer(kind=c_int), value,               intent(in)    :: natoms
        real(kind=c_double),                      intent(out)   :: E_i
        real(kind=c_double), dimension(3,natoms), intent(inout) :: F
        real(kind=c_double), dimension(3,3), intent(inout) :: S
        integer(kind=c_int),                      intent(out)   :: stat
        integer                                         :: x1,x2,d

        
        real(real64)::E_i_arr(1:1)
        real(real64)::dEdE(1:1)
        real(real64),allocatable::dEdc(:)
        integer::j
        real(real64),allocatable::Ftest(:,:)

        INTEGER :: start_clock, end_clock, clock_rate
        REAL(real64) :: elapsed_time
        integer::it


        !write(*,*) "energy"
        !write(*,*) index_i,&
        !coo_i, &
        !n_j, &
        !coo_j, &
        !type_j, &
        !BPpots(index_i)%ltype
        dEdE = 1d0
        allocate(dEdc(1:BPpots(type_i)%nsf))

        
        call BPpots(type_i)%fingerprint_basis%evaluate(type_i,&
            coo_i, &
            n_j, &
            coo_j, &
            type_j, &
            BPpots(type_i)%ltype,&
            values=fingerprint_values(1:BPpots(type_i)%nsf), &
            deriv_i=fingerprint_deriv_i(:,1:BPpots(type_i)%nsf), &
            deriv_j=fingerprint_deriv_j(:,1:BPpots(type_i)%nsf,:))
        


        call BPpots(type_i)%normalize(fingerprint_values(1:BPpots(type_i)%nsf),&
                fingerprint_deriv_i(:,1:BPpots(type_i)%nsf), &
                fingerprint_deriv_j(:,1:BPpots(type_i)%nsf,:),&
                n_j)

        !write(*,*)  "d",sum(fingerprint_deriv_j(:,1:BPpots(type_i)%nsf,:))

        !stop
        
        !write(*,*) "fingerprint_values",fingerprint_values(1:BPpots(type_i)%nsf)
        !stop

        call BPpots(type_i)%forward(fingerprint_values(1:BPpots(type_i)%nsf),E_i_arr)


    
        



        call BPpots(type_i)%backward(dEdE,dEdc)


        !write(*,*) E_i_arr
        !stop

        E_i = BPpots(type_i)%E_scale*E_i_arr(1) + BPpots(type_i)%shift
        !write(*,*) E_i,BPpots(type_i)%E_scale,BPpots(type_i)%shift
        E_i = E_i + BPpots(type_i)%E_atom(type_i)
        !write(*,*) E_i
        !stop
        !write(*,*) dEdc
        !call stp_eval(type_i, coo_i, n_j, coo_j, type_j, &
        !          aenet_pot(type_i)%stp, sfval=sfval, &
        !          sfderiv_i=sfderiv_i, sfderiv_j=sfderiv_j, scaled=.true.)
        F(1:3, index_i) = F(1:3, index_i) -  BPpots(type_i)%E_scale &
                    * matmul(fingerprint_deriv_i(1:3,1:BPpots(type_i)%nsf), dEdc(1:BPpots(type_i)%nsf))
        !write(*,*) F(1:3, index_i)
        

        !allocate(Ftest(3,n_j))
        do j = 1, n_j
            !Ftest(1:3,j) = -  BPpots(index_i)%scale &
            !        * matmul(fingerprint_deriv_j(1:3,1:BPpots(index_i)%nsf,j), dEdc(1:BPpots(index_i)%nsf))

            F(1:3, index_j(j)) = F(1:3, index_j(j)) -  BPpots(type_i)%E_scale &
                            * matmul(fingerprint_deriv_j(1:3,1:BPpots(type_i)%nsf,j), dEdc(1:BPpots(type_i)%nsf))
        end do

        !stop
    !stress tensor
        do x1=1,3
            do x2 = x1,3
               do j =1, n_j
                  do d = 1, BPpots(type_i)%nsf
                     S(x2,x1) = S(x2,x1) + (coo_i(x2)-coo_j(x2,j))*BPpots(type_i)%E_scale&
                          * fingerprint_deriv_j(x1,d,j)*dEdc(d)
                  end do
               end do
            end do
        end do


        !call check_force()
        !stop

        contains 

        subroutine check_force()
            implicit none
            real(real64)::value1(1:BPpots(type_i)%nsf),value2(1:BPpots(type_i)%nsf)
            real(real64),allocatable:: coo_i2(:),coo_j2(:,:)
            real(real64)::eta,E_i_arr2(1:1),E_i2
            integer::ii

            allocate(coo_i2,source= coo_i)
            allocate(coo_j2,source= coo_j)
            eta = 1d-6

            do ii=1,3
                coo_i2 = coo_i
                coo_i2(ii) = coo_i2(ii) +eta
                call BPpots(type_i)%fingerprint_basis%evaluate(type_i,&
                    coo_i2, &
                    n_j, &
                    coo_j, &
                    type_j, &
                    BPpots(type_i)%ltype,&
                    values=value1,&
                    deriv_i=fingerprint_deriv_i(:,1:BPpots(type_i)%nsf), &
                    deriv_j=fingerprint_deriv_j(:,1:BPpots(type_i)%nsf,:))
                
                call BPpots(type_i)%forward(value1,E_i_arr2)
                E_i2 = BPpots(type_i)%E_scale*E_i_arr2(1) + BPpots(type_i)%shift
                E_i2 = E_i2 + BPpots(type_i)%E_atom(type_i)
                write(*,*) (E_i-E_i2)/eta,F(ii,index_i)
            end do
            coo_i2 = coo_i
            do j = 1, n_j
                do ii=1,3
                    coo_j2 = coo_j
                    coo_j2(ii,j) = coo_j2(ii,j) + eta
                    call BPpots(type_i)%fingerprint_basis%evaluate(type_i,&
                        coo_i2, &
                        n_j, &
                        coo_j2, &
                        type_j, &
                        BPpots(type_i)%ltype,&
                        values=value1,&
                        deriv_i=fingerprint_deriv_i(:,1:BPpots(type_i)%nsf), &
                        deriv_j=fingerprint_deriv_j(:,1:BPpots(type_i)%nsf,:))
                
                    call BPpots(type_i)%forward(value1,E_i_arr2)
                    E_i2 = BPpots(type_i)%E_scale*E_i_arr2(1) + BPpots(type_i)%shift
                    E_i2 = E_i2 + BPpots(type_i)%E_atom(index_i)
                    write(*,*) E_i,E_i2, E_i-E_i2
                    write(*,*)j,ii,(E_i2-E_i)/eta,Ftest(ii,j),(E_i-E_i2)/eta-Ftest(ii,j)!,F(ii,index_j(j))
                end do
            end do



            
        end subroutine

    end subroutine



end module shared_potentials