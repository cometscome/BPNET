module fingerprints
    use iso_fortran_env
    use ioutils,only:read_next_valid_line,&
            to_lowercase,extract_key_value
    use fingerprint_basis,only:FingerPrintBasis
    use chebyshevbasis_fp,only:FP_ChebyshevBasis
    use splinebasis_fp,only:FP_SplineBasis
    use LJbasis_fp,only:FP_LJBasis
    use Multibasis_fp,only:FP_MultiBasis

    !use chebyshevbasis,only:read_basis_chebyshev
    implicit none
    integer,parameter::num_of_maxparameters = 4 !NSFPARAM 
    ! NENV_MAX    maximum number of types involved in single function  
    integer, parameter :: NENV_MAX = 2

    public::get_Rcmax_and_Rcmin,initialize_fingerprint_basisset

    type, public :: FingerPrint
        real(real64)                                             :: Rc_min
        real(real64)                                             :: Rc_max
        integer                                             :: num_of_coeffs !nsf
        character(len=100)                                  :: fingerprint_type !sftype
        integer                                             :: num_of_different_species !nenv
        character(len=2)                                    :: species_of_central_atom !atomtype 
        character(len=2), dimension(:),   allocatable :: species_of_surrounding_atom !envtype
        character(len=1024)                                 :: description
        integer                                             :: num_of_parameters !nsfparam
        real(real64),       dimension(:,:), allocatable :: fingerprint_parameters !sfparam
        integer                                             :: num_evaluations !neval

        real(real64),       dimension(:),   allocatable :: fingerprint_values_min !sfval_min
        real(real64),       dimension(:),   allocatable :: fingerprint_values_max !sfval_max
        real(real64),       dimension(:),   allocatable :: fingerprint_values_avg !sfval_avg
        real(real64),       dimension(:),   allocatable :: fingerprint_values_cov !sfval_cov   
        
        ! gtype(i)      global atom type ID for local type i              !
        ! ltype(i)      local atom type ID for global type i              !        
        integer,                dimension(:),   allocatable :: gtype
        integer,                dimension(:),   allocatable :: ltype
        integer                                             :: ntypes_global

        !Behler-Parrinello basis
        ! function kind of the particular basis type  
        integer,                dimension(:),   allocatable :: fingerprint_function_kind !sf
        ! i-th environment species for j-th basis function
        integer,                dimension(:,:), allocatable :: species_of_basis_functions!sfenv

        contains
        procedure:: read => read_FingerPrint
        procedure:: allocate_statistics_arrays
        procedure:: set_global_types
        procedure:: print_info
        procedure:: write => write_FingerPrint
        procedure:: count_evaluation
        procedure:: update_statistics => update_statistics_FingerPrint
    end type FingerPrint

    interface FingerPrint
        module procedure init_FingerPrint
    end interface FingerPrint

    interface
        module subroutine read_basis_chebyshev(fp,line)
            type(FingerPrint),intent(inout)::fp
            character(len=*),intent(in)::line
        end subroutine read_basis_chebyshev

        module subroutine print_info_chebyshev(fingerprint_parameters)
            real(real64),intent(in)::fingerprint_parameters(:,:)
        end subroutine

        module subroutine read_chebyshevdata(line,paramdata,num_of_coeffs,Rc)
            character(len=*),intent(in)::line
            real(real64),intent(out)::paramdata(:)
            integer,intent(out)::num_of_coeffs
            real(real64),intent(out)::Rc
        end subroutine

        module subroutine read_basis_spline(fp,line)
            type(FingerPrint),intent(inout)::fp
            character(len=*),intent(in)::line
        end subroutine read_basis_spline

        module subroutine print_info_spline(fingerprint_parameters)
            !type(FingerPrint),intent(in)::fp
            real(real64),intent(in)::fingerprint_parameters(:,:)
        end subroutine

        module subroutine read_splinedata(line,paramdata,num_of_coeffs,Rc)
            character(len=*),intent(in)::line
            real(real64),intent(out)::paramdata(:)
            integer,intent(out)::num_of_coeffs
            real(real64),intent(out)::Rc
        end subroutine

        module subroutine read_basis_LJ(fp,line)
            type(FingerPrint),intent(inout)::fp
            character(len=*),intent(in)::line
        end subroutine read_basis_LJ

        module subroutine print_info_LJ(fingerprint_parameters)
            real(real64),intent(in)::fingerprint_parameters(:,:)
            !type(FingerPrint),intent(in)::fp
        end subroutine

        module subroutine read_LJdata(line,ntypes,paramdata,num_of_coeffs,Rc)
            character(len=*),intent(in)::line
            real(real64),intent(out)::paramdata(:)
            integer,intent(out)::num_of_coeffs
            integer,intent(in)::ntypes
            real(real64),intent(out)::Rc
        end subroutine

        module subroutine read_basis_multi(fp,unit,num_kinds)
            type(FingerPrint),intent(inout)::fp
            integer,intent(in)::unit
            integer,intent(in)::num_kinds
        end subroutine read_basis_multi

        module subroutine print_info_multi(fingerprint_parameters)
            !type(FingerPrint),intent(in)::fp
            real(real64),intent(in)::fingerprint_parameters(:,:)
        end subroutine

        module subroutine read_basis_splinetensor(fp,line)
            type(FingerPrint),intent(inout)::fp
            character(len=*),intent(in)::line
        end subroutine read_basis_splinetensor

        module subroutine print_info_splinetensor(fingerprint_parameters)
            !type(FingerPrint),intent(in)::fp
            real(real64),intent(in)::fingerprint_parameters(:,:)
        end subroutine

        module subroutine read_splinetensordata(line,paramdata,num_of_coeffs,Rc)
            character(len=*),intent(in)::line
            real(real64),intent(out)::paramdata(:)
            integer,intent(out)::num_of_coeffs
            real(real64),intent(out)::Rc
        end subroutine

    end interface

    contains

    type(FingerPrint) function init_FingerPrint() result(fp)
        fp%Rc_min = 0.0d0
        fp%Rc_max = 0.0d0
        fp%num_of_coeffs = 0
        fp%num_of_different_species = 0
        
    end function init_FingerPrint

    subroutine print_info(self)
        implicit none
        class(FingerPrint) ::self
        integer::i


        write(*,*) 'Structural fingerprint (SF) set-up for ', &
            trim(adjustl(self%species_of_central_atom))
        write(*,*)
        call print_descr(self%description)


        if (self%num_of_different_species >0) then
            write(*,'(1x,"environment types: ")', advance='no')
            do i = 1, self%num_of_different_species
                write(*,'(A2,1x)', advance='no') self%species_of_surrounding_atom(i)
            end do
        end if

        write(*,*)
        write(*,*) 'minimal distance : ', self%Rc_min, ' Angstrom'
        write(*,*) 'maximal cut-off  : ', self%Rc_max, ' Angstrom'
        write(*,*) 'size of basis    : ', self%num_of_coeffs
!        write(*,*) 'evaluations      : ', self%
        write(*,*)

        select case(trim(to_lowercase(self%fingerprint_type)))
        case('chebyshev')
           call print_info_chebyshev(self%fingerprint_parameters)
        case('spline')
            call print_info_spline(self%fingerprint_parameters)
        case('lj')
            call print_info_LJ(self%fingerprint_parameters)
        case('multi')
            call print_info_multi(self%fingerprint_parameters)
        end select

    end subroutine print_info

    subroutine count_evaluation(self)
        implicit none
        class(FingerPrint) ::self
        self%num_evaluations = self%num_evaluations + 1
        
    end subroutine count_evaluation


    subroutine update_statistics_FingerPrint(self,values)
        implicit none
        class(FingerPrint) ::self
        real(real64),intent(in)::values(:)
        integer::nsf
        integer::i

        nsf = self%num_of_coeffs

        if (self%num_evaluations == 0) then
            self%fingerprint_values_min(1:nsf) = values(1:nsf)
            self%fingerprint_values_max(1:nsf) = values(1:nsf)
            self%fingerprint_values_avg(1:nsf) = values(1:nsf)
            self%fingerprint_values_cov(1:nsf) = values(1:nsf)*values(1:nsf)
        else
            do i = 1, nsf
                self%fingerprint_values_min(i) = min(self%fingerprint_values_min(i), values(i))
                self%fingerprint_values_max(i) = max(self%fingerprint_values_max(i), values(i))
                self%fingerprint_values_avg(i) = (dble(self%num_evaluations)*self%fingerprint_values_avg(i) &
                                  + values(i))/(dble(self%num_evaluations+1))
                self%fingerprint_values_cov(i) = (dble(self%num_evaluations)*self%fingerprint_values_cov(i) &
                                  + values(i)*values(i))/(dble(self%num_evaluations+1))
            end do
        end if
        
    end subroutine update_statistics_FingerPrint

    subroutine print_descr(descr)
        implicit none
        character(len=*), intent(in) :: descr

        integer :: i1, i2, l

        l = len_trim(descr)

        i1 = 1
        i2 = index(descr,'$$')
        do while(i2 >= i1)
        write(*,*) descr(i1:i2-1)
        i1 = i2 + 2
        i2 = i1 + index(descr(i1:l),'$$') - 1
        end do
        if (len_trim(descr(i1:l)) > 0) then
        write(*,*) trim(descr(i1:l))
        end if
        write(*,*)

    end subroutine print_descr

    subroutine allocate_statistics_arrays(self,n)
        implicit none
        class(FingerPrint) ::self
        integer,intent(in)::n

        allocate( &
            self%fingerprint_values_min(n), &
            self%fingerprint_values_max(n), &
            self%fingerprint_values_avg(n), &
            self%fingerprint_values_cov(n) &
        )

        self%fingerprint_values_min = 0d0
        self%fingerprint_values_max = 0d0
        self%fingerprint_values_avg = 0d0
        self%fingerprint_values_cov = 0d0
        
    end subroutine 

    subroutine write_FingerPrint(self, file, unit)
        implicit none
        class(FingerPrint) ::self
        character(len=*), optional, intent(in) :: file
        integer,          optional, intent(in) :: unit
        integer :: u

        if (present(unit)) then
            u = unit
        else if (present(file)) then
            open(newunit=u, file=trim(file), status='replace', action='write', &
                    form='unformatted')
        else
            write(0,*) "Error: neither unit number nor file name given"
            write(0,*) "       in `save_Setup'."
            return
        end if
        

        write(u) self%description
        write(u) self%species_of_central_atom ! %atomtype
        write(u) self%num_of_different_species !%nenv
        write(u) self%species_of_surrounding_atom !%envtypes(:)
        write(u) self%Rc_min
        write(u) self%Rc_max
        write(u) self%fingerprint_type!sftype
        write(u) self%num_of_coeffs!nsf
        write(u) self%num_of_parameters !nsfparam
        
        write(u) self%fingerprint_function_kind(:)
        !write(*,*) size(self%fingerprint_parameters(:,:) ,1),size(self%fingerprint_parameters(:,:) ,2)
        !stop
        write(u) self%fingerprint_parameters(:,:) !sfparam(:,:)

        write(u) self%species_of_basis_functions(:,:) !sfenv(:,:)
        write(u) self%num_evaluations !neval
        
        write(u) self%fingerprint_values_min ! sfval_min(:)
        write(u) self%fingerprint_values_max !sfval_max(:)
        write(u) self%fingerprint_values_avg !sfval_avg(:)
        write(u) self%fingerprint_values_cov !sfval_cov(:)

        if (.not. present(unit)) then
           close(u)
        end if


    end subroutine write_FingerPrint

    

    subroutine read_FingerPrint(fp,inFile, global_types)
        implicit none
        character(len=*),               intent(in) :: inFile
        character(len=*), dimension(:), intent(in) :: global_types
        class(FingerPrint) ::fp

        integer :: unit, ios
        logical :: eof
        character(len=20)            :: keyword
        integer,           parameter :: nc = 1024
        character(len=nc)            :: line
        integer::i
        integer                      :: ntypes_global
        integer::num_kinds,ikind

        ntypes_global = size(global_types(:))

        fp%Rc_min = 0.0d0
        fp%Rc_max = 0.0d0
        fp%num_evaluations  = 0
        fp%num_of_coeffs = 0

        ! ファイルを開く
        open(newunit=unit, file=trim(adjustl(inFile)), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'Error: Could not open file ',inFile
            stop
        end if

        do
            call read_next_valid_line(unit, line, eof)
            if (eof) exit

            read(line,*) keyword
            select case(to_lowercase(trim(adjustl(keyword))))
            case('atom')
                read(line(5:nc), *) fp%species_of_central_atom
            case('env')
                read(line(4:nc), *) fp%num_of_different_species
                allocate(fp%species_of_surrounding_atom(fp%num_of_different_species))
                do i = 1, fp%num_of_different_species
                    call read_next_valid_line(unit, line, eof)
                    read(line,*) fp%species_of_surrounding_atom(i)
                end do
            case('descr')
                fp%description = ' '
                read(unit, '(A)') line
                do while(index(to_lowercase(trim(adjustl(line))),'end descr') <= 0)
                    if (len_trim(fp%description) < len(fp%description)) then
                        fp%description = trim(fp%description) // trim(line) // "$$"
                    end if
                    read(unit, '(A)') line
                end do
            case('rmin')
                read(line(5:nc), *) fp%Rc_min
            case('symmfunc', 'functions')
                write(*,*) "not implemented yet"
            case('basis')
                if (fp%num_of_coeffs > 0) then
                   write(0,*) "Error: Multiple basis definitions found in " // &
                              "structural fingerprint setup."
                   stop
                end if
                !write(*,*) line
                call extract_key_value(line, 'type', fp%fingerprint_type)

                select case(to_lowercase(fp%fingerprint_type))
                case('chebyshev')
                    call read_next_valid_line(unit, line, eof)
                    !write(*,*) line
                    call read_basis_chebyshev(fp,line)
                !    call read_basis_chebyshev(u_stp, stp, iline)     
                case('spline')
                    call read_next_valid_line(unit, line, eof)    
                    
                    call read_basis_spline(fp,line)
                case('lj')
                    call read_next_valid_line(unit, line, eof)    
                    
                    call read_basis_LJ(fp,line)
                case('multi')
                    call read_next_valid_line(unit, line, eof)  
                    read(line, *) num_kinds
                    call read_basis_multi(fp,unit,num_kinds)
                    !call read_next_valid_line(unit, line, eof)    
                    
                    !call read_basis_LJ(fp,line)
                case default
                    write(0,*) "read_FingerPrint: Error: Unknown basis type: ", trim(fp%fingerprint_type)
                    !deallocate(stp%sf, stp%sfparam)
                    stop
                end select
            case default
                write(0,*) "Error: Unknown keyword : ", trim(keyword)
                write(0,*) "       file            : ", trim(inFile)
                !write(0,*) "       line            : ", trim(io_adjustl(iline))
            end select
            !print *, 'Valid Line:', trim(line)
        end do

        close(unit)

        if (fp%Rc_min <= 1d-16) then
            write(0,*) "Warning: RMIN not set in : ", trim(inFile)
            write(0,*) "         --> setting Rc_min = 1.0"
            fp%Rc_min = 1.0d0
        end if
     
        if (fp%Rc_max <= 1d-16) then
            write(0,*) "Error: Cutoff undetermined not set in : ", trim(inFile)
            stop
        end if

        !write(*,*) global_types
        allocate(fp%gtype(fp%num_of_different_species), &
                fp%ltype(ntypes_global))
        call fp%set_global_types(ntypes_global, global_types)
        fp%ntypes_global = ntypes_global

    end subroutine read_FingerPrint

    subroutine initialize_fingerprint_basisset_fromparam(&
        fingerprint_type,&
        fingerprint_parameters,&
        num_of_different_species,&
        species_of_surrounding_atom ,&
        num_of_types,&
        fingerprint_basis_set)
        implicit none
        character(len=*),intent(in)                                  :: fingerprint_type !sftype
        real(real64),       dimension(:,:,:), intent(in) :: fingerprint_parameters !sfparam
        integer,intent(in)                                             :: num_of_different_species(:) !nenv
        character(len=2), dimension(:,:),  intent(in) :: species_of_surrounding_atom !envtype
        integer,                        intent(in)  :: num_of_types
        !integer,                        intent(in)  :: num_of_maxatoms_in_sphere
        class(FingerPrintBasis),dimension(:), allocatable,intent(out)  :: fingerprint_basis_set
        !integer,intent(out)::num_max_coefficients
        integer::itype
        type(FP_ChebyshevBasis), dimension(:), allocatable :: chebyshev_basis_set
        type(FP_SplineBasis), dimension(:), allocatable :: spline_basis_set
        type(FP_LJBasis), dimension(:), allocatable :: LJ_basis_set
        type(FP_MultiBasis), dimension(:), allocatable :: multi_basis_set
        



        select case(trim(to_lowercase(fingerprint_type)))
        case('chebyshev')
           allocate(chebyshev_basis_set(num_of_types))
           do itype = 1, num_of_types
                chebyshev_basis_set(itype)= FP_ChebyshevBasis(fingerprint_parameters(:,:,itype),&
                    num_of_different_species(itype),&
                    species_of_surrounding_atom(:,itype))
                
               !call setup_basis_chebyshev(stp(itype), sfb(itype))
           end do
           allocate(fingerprint_basis_set,source=chebyshev_basis_set)
        case('spline')
            allocate(spline_basis_set(num_of_types))
            do itype = 1, num_of_types
                spline_basis_set(itype)= FP_SplineBasis(fingerprint_parameters(:,:,itype),&
                     num_of_different_species(itype),&
                     species_of_surrounding_atom(:,itype))
                 
                !call setup_basis_chebyshev(stp(itype), sfb(itype))
            end do
            allocate(fingerprint_basis_set,source=spline_basis_set)
        case('lj')
            allocate(LJ_basis_set(num_of_types))
            do itype = 1, num_of_types
                LJ_basis_set(itype)= FP_LJBasis(fingerprint_parameters(:,:,itype),&
                     num_of_different_species(itype),&
                     species_of_surrounding_atom(:,itype))
                 
                !call setup_basis_chebyshev(stp(itype), sfb(itype))
            end do
            allocate(fingerprint_basis_set,source=LJ_basis_set)
        case('multi')
            allocate(multi_basis_set(num_of_types))
            do itype = 1, num_of_types
                multi_basis_set(itype)= FP_MultiBasis(fingerprint_parameters(:,:,itype),&
                     num_of_different_species(itype),&
                     species_of_surrounding_atom(:,itype),num_of_maxparameters)
                !call setup_basis_chebyshev(stp(itype), sfb(itype))
            end do
            allocate(fingerprint_basis_set,source=multi_basis_set)
            
        case('behler2011')
            write(*,*) "behler2011 is not implemeted yet!!"
        case default
            write(0,*) "Error: Unknown basis function type : ", trim(fingerprint_type)
            stop
        end select

        !num_max_coefficients = get_max_number_of_coefficients(fingerprint_set) 

    end subroutine

    subroutine initialize_fingerprint_basisset(fingerprint_set,num_of_types,&
                fingerprint_basis_set,&
                num_max_coefficients)
        implicit none
        type(FingerPrint),dimension(:),intent(in) ::fingerprint_set
        integer,                        intent(in)  :: num_of_types
        !integer,                        intent(in)  :: num_of_maxatoms_in_sphere
        class(FingerPrintBasis),dimension(:), allocatable,intent(out)  :: fingerprint_basis_set
        integer,intent(out)::num_max_coefficients
        type(FP_ChebyshevBasis), dimension(:), allocatable :: chebyshev_basis_set
        type(FP_SplineBasis), dimension(:), allocatable :: spline_basis_set
        type(FP_LJBasis), dimension(:), allocatable :: LJ_basis_set
        type(FP_MultiBasis), dimension(:), allocatable :: multi_basis_set
        
        character(len=100) :: fingerprint_type
        integer::itype
        !real(real64),       dimension(:,:,:), allocatable :: fingerprint_parameters !sfparam
        !integer,allocatable                                         :: num_of_different_species(:) !nenv
        !character(len=2), dimension(:,:),   allocatable :: species_of_surrounding_atom !envtype
        !integer::n1,n2,n3

        !n1 = size(fingerprint_set(1)%fingerprint_parameters,1)
        !n2 =  size(fingerprint_set(1)%fingerprint_parameters,2)
        !allocate(fingerprint_parameters(n1,n2,num_of_types))
        !allocate(num_of_different_species(num_of_types))
        !n1 = size(fingerprint_set(1)%species_of_surrounding_atom,1)
        !allocate( species_of_surrounding_atom(n1,num_of_types))

        !do itype=1,num_of_types
        !    fingerprint_parameters(:,:,itype) = fingerprint_set(itype)%fingerprint_parameters
        !    num_of_different_species(itype) = fingerprint_set(itype)%num_of_different_species
        !    species_of_surrounding_atom(:,itype) = fingerprint_set(itype)%species_of_surrounding_atom
        !end do
        

        fingerprint_type = trim(fingerprint_set(1)%fingerprint_type)
        do itype = 1, num_of_types
            if (trim(fingerprint_set(itype)%fingerprint_type) /= trim(fingerprint_type)) then
                write(0,*) "Error: Mixing of basis functions of different " &
                                // "types not yet implemented."
                write(0,*) trim(fingerprint_type), trim(fingerprint_set(itype)%fingerprint_type)
                stop
            end if
        end do


        !call initialize_fingerprint_basisset_fromparam(&
        !            fingerprint_type,&
        !            fingerprint_parameters,&
        !!            num_of_different_species,&
        !            species_of_surrounding_atom ,&
        !            num_of_types,&
        !            fingerprint_basis_set,&
        !            num_max_coefficients)
                    
        !return 

        select case(trim(to_lowercase(fingerprint_type)))
        case('chebyshev')
           allocate(chebyshev_basis_set(num_of_types))
           do itype = 1, num_of_types
                chebyshev_basis_set(itype)= FP_ChebyshevBasis(fingerprint_set(itype)%fingerprint_parameters,&
                    fingerprint_set(itype)%num_of_different_species,&
                    fingerprint_set(itype)%species_of_surrounding_atom )
                
               !call setup_basis_chebyshev(stp(itype), sfb(itype))
           end do
           allocate(fingerprint_basis_set,source=chebyshev_basis_set)
        case('spline')
            allocate(spline_basis_set(num_of_types))
            do itype = 1, num_of_types
                 spline_basis_set(itype)= FP_SplineBasis(fingerprint_set(itype)%fingerprint_parameters,&
                     fingerprint_set(itype)%num_of_different_species,&
                     fingerprint_set(itype)%species_of_surrounding_atom )
                 
                !call setup_basis_chebyshev(stp(itype), sfb(itype))
            end do
            allocate(fingerprint_basis_set,source=spline_basis_set)
        case('lj')
            allocate(LJ_basis_set(num_of_types))
            do itype = 1, num_of_types
                LJ_basis_set(itype)= FP_LJBasis(fingerprint_set(itype)%fingerprint_parameters,&
                        fingerprint_set(itype)%num_of_different_species,&
                        fingerprint_set(itype)%species_of_surrounding_atom )
                 
                !call setup_basis_chebyshev(stp(itype), sfb(itype))
            end do
            allocate(fingerprint_basis_set,source=LJ_basis_set)
        case('multi')
            allocate(multi_basis_set(num_of_types))
            do itype = 1, num_of_types
                multi_basis_set(itype)= FP_MultiBasis(fingerprint_set(itype)%fingerprint_parameters,&
                        fingerprint_set(itype)%num_of_different_species,&
                        fingerprint_set(itype)%species_of_surrounding_atom ,num_of_maxparameters)
                !call setup_basis_chebyshev(stp(itype), sfb(itype))
            end do
            allocate(fingerprint_basis_set,source=multi_basis_set)
        case('behler2011')
            write(*,*) "behler2011 is not implemeted yet!!"
        case default
            write(0,*) "Error: Unknown basis function type : ", trim(fingerprint_type)
            stop
        end select

        num_max_coefficients = get_max_number_of_coefficients(fingerprint_set) 

    end subroutine initialize_fingerprint_basisset

    integer function get_max_number_of_coefficients(fingerprint_set) result(num_max_coefficients)
        implicit none
        type(FingerPrint),dimension(:),intent(in) ::fingerprint_set
        integer::num_of_types
        integer::itype

        num_of_types = size(fingerprint_set(:))
        num_max_coefficients  = 0
        do itype = 1, num_of_types
            num_max_coefficients = max(num_max_coefficients,fingerprint_set(itype)%num_of_coeffs)
        end do
    end function

    subroutine set_global_types(self, ntypes_global, global_types)

        implicit none
    
        class(FingerPrint) ::self
        integer,                                    intent(in)    :: ntypes_global
        character(len=*), dimension(ntypes_global), intent(in)    :: global_types
    
        integer :: i, j
    
        ! atom type indices that have no corresponding entry are set to 0
    
        do i = 1, ntypes_global
           self%ltype(i) = 0
           env : do j = 1, self%num_of_different_species
              if (trim(global_types(i)) == trim(self%species_of_surrounding_atom(j))) then
                 self%ltype(i) = j
                 exit env
              end if
           end do env
        end do
    
        ! reverse direction, because we do not know, if the sets of types
        ! are identical
        do i = 1, self%num_of_different_species
           self%gtype(i) = 0
           global : do j = 1, ntypes_global
              if (trim(global_types(j)) == trim(self%species_of_surrounding_atom(i))) then
                 self%gtype(i) = j
                 exit global
              end if
           end do global
        end do
    
    end subroutine set_global_types

    subroutine get_Rcmax_and_Rcmin(fps,num_of_types,Rc_min,Rc_max)
        type(FingerPrint),dimension(:),intent(in) ::fps
        integer,intent(in)::num_of_types
        real(real64),intent(out)::Rc_max,Rc_min
        integer::itype
        

        Rc_min = fps(1)%Rc_min
        Rc_max = fps(1)%Rc_max

        do itype = 1, num_of_types
            Rc_min = min(Rc_min, fps(itype)%Rc_min)
            Rc_max = max(Rc_max, fps(itype)%Rc_max)
        end do
            
    end subroutine get_Rcmax_and_Rcmin



end module fingerprints