submodule (fingerprints) LJbasis
    use iso_fortran_env
    implicit none
    contains

    module subroutine read_LJdata(line,ntypes,paramdata,num_of_coeffs,Rc)
        use bspline,only:d => bspline_order
        character(len=*),intent(in)::line
        real(real64),intent(out)::paramdata(:)
        integer,intent(in)::ntypes
        integer,intent(out)::num_of_coeffs
        integer             :: r_N, a_N!,ntypes
        real(real64)    :: r_Rc, a_Rc
        real(real64),intent(out)::Rc

        call extract_key_value(line, 'radial_Rc', r_Rc)
        !call extract_key_value(line, 'angular_Rc', a_Rc)
        !nTypes = fp%num_of_different_species

        num_of_coeffs = 2*nTypes

        paramdata(1) = r_Rc
        paramdata(2) = 0d0
        paramdata(3) = 0d0
        paramdata(4) = 0d0
        Rc = r_Rc
    end

    module subroutine read_basis_LJ(fp,line)
        implicit none
        type(FingerPrint),intent(inout)::fp
        character(len=*),intent(in)::line
        !integer             :: r_N, a_N
        real(real64)    :: r_Rc, a_Rc
        integer::nTypes
        real(real64)::paramdata(4)


        !r_Rc = 0.0d0
        !a_Rc = 0.0d0
        nTypes = fp%num_of_different_species
        call read_LJdata(line,nTypes,paramdata,fp%num_of_coeffs,fp%Rc_max)
        !r_Rc = paramdata(1)
        !a_Rc = paramdata(3)

        !call extract_key_value(line, 'radial_Rc', r_Rc)
        !call extract_key_value(line, 'angular_Rc', a_Rc)
        !nTypes = fp%num_of_different_species

        !fp%num_of_coeffs = num_of_coeffs!2*nTypes

        call fp%allocate_statistics_arrays(fp%num_of_coeffs)
        fp%num_of_parameters = num_of_maxparameters
        allocate(fp%fingerprint_parameters(num_of_maxparameters,fp%num_of_coeffs),source=0d0)

        allocate(fp%fingerprint_function_kind(fp%num_of_coeffs),source=0)
        allocate(fp%species_of_basis_functions(NENV_MAX,fp%num_of_coeffs),source=0)


        fp%fingerprint_parameters(:,1) = paramdata(:)
        !fp%fingerprint_parameters(1,1) = r_Rc
        !fp%fingerprint_parameters(3,1) = a_Rc

        !fp%Rc_max = r_Rc!max(r_Rc, a_Rc)

    end subroutine read_basis_LJ


    !module subroutine print_info_LJ(fp)
    !    implicit none
    !    type(FingerPrint), intent(in) :: fp
    
    !    double precision :: r_Rc, a_Rc
    !    !integer          :: r_N, a_N
    
    !    r_Rc = fp%fingerprint_parameters(1,1)
    !    a_Rc = fp%fingerprint_parameters(3,1)
    
    !    write(*,*) 'Basis function type LJ'
    !    write(*,*) '[Y. Nagai and M. Okumura (2025)]'
    !    write(*,*)
    !    write(*,*) 'Radial Rc     : ', r_Rc
    !    write(*,*) 'Angular Rc    : ', a_Rc
        !write(*,*) 'Radial order  : ', r_N
        !write(*,*) 'Angular order : ', a_N
    !    write(*,*)
    
    !end subroutine print_info_LJ

    module subroutine print_info_LJ(fingerprint_parameters)
        implicit none
        real(real64), intent(in) :: fingerprint_parameters(:,:)
    
        double precision :: r_Rc, a_Rc
        !integer          :: r_N, a_N
    
        r_Rc = fingerprint_parameters(1,1)
        !a_Rc = fingerprint_parameters(3,1)
    
        write(*,*) 'Basis function type LJ'
        write(*,*) '[Y. Nagai and M. Okumura (2025)]'
        write(*,*)
        write(*,*) 'Radial Rc     : ', r_Rc
        !write(*,*) 'Angular Rc    : ', a_Rc
        !write(*,*) 'Radial order  : ', r_N
        !write(*,*) 'Angular order : ', a_N
        write(*,*)
    
    end subroutine print_info_LJ


end submodule LJbasis