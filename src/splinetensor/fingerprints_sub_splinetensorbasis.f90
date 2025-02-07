submodule (fingerprints) splinetensorbasis
    use iso_fortran_env
    implicit none
    contains


    module subroutine read_splinetensordata(line,paramdata,num_of_coeffs,Rc)
        use bspline,only:d => bspline_order
        character(len=*),intent(in)::line
        real(real64),intent(out)::paramdata(:)
        integer,intent(out)::num_of_coeffs
        real(real64),intent(out)::Rc
        integer             :: r_N, a_N
        real(real64)    :: r_Rc, a_Rc

        call extract_key_value(line, 'radial_Rc', r_Rc)
        call extract_key_value(line, 'radial_N', r_N)
        call extract_key_value(line, 'angular_Rc', a_Rc)
        call extract_key_value(line, 'angular_N', a_N)

        num_of_coeffs = (r_N +2*d -d -1) + (a_N +2*d -d -1)

        paramdata(1) = r_Rc
        paramdata(2) = dble(r_N)
        paramdata(3) = a_Rc
        paramdata(4) = dble(a_N)

        Rc = max(r_Rc, a_Rc)
    end

    module subroutine read_basis_splinetensor(fp,line)
        use bspline,only:d => bspline_order
        implicit none
        type(FingerPrint),intent(inout)::fp
        character(len=*),intent(in)::line
        !integer             :: r_N, a_N
        !real(real64)    :: r_Rc, a_Rc
        real(real64)::paramdata(4)

        call read_splinetensordata(line,paramdata,fp%num_of_coeffs,fp%Rc_max)


        !r_Rc = 0.0d0
        !a_Rc = 0.0d0
        !r_N = 0
        !a_N = 0

        !call extract_key_value(line, 'radial_Rc', r_Rc)
        !call extract_key_value(line, 'radial_N', r_N)
        !call extract_key_value(line, 'angular_Rc', a_Rc)
        !call extract_key_value(line, 'angular_N', a_N)

        !write(*,*) r_Rc,r_N,a_Rc,a_N
        !fp%num_of_coeffs = r_N + a_N + 2
        !fp%num_of_coeffs = num_of_coeffs!(r_N +2*d -d -1) + (a_N +2*d -d -1)

        if (fp%num_of_different_species > 1) then
            fp%num_of_coeffs = 2*fp%num_of_coeffs
        end if

        call fp%allocate_statistics_arrays(fp%num_of_coeffs)
        fp%num_of_parameters = num_of_maxparameters
        allocate(fp%fingerprint_parameters(num_of_maxparameters,fp%num_of_coeffs),source=0d0)

        allocate(fp%fingerprint_function_kind(fp%num_of_coeffs),source=0)
        allocate(fp%species_of_basis_functions(NENV_MAX,fp%num_of_coeffs),source=0)

        fp%fingerprint_parameters(:,1) = paramdata(:)

        !fp%fingerprint_parameters(1,1) = r_Rc
        !fp%fingerprint_parameters(2,1) = dble(r_N)
        !fp%fingerprint_parameters(3,1) = a_Rc
        !fp%fingerprint_parameters(4,1) = dble(a_N)

        
    end subroutine read_basis_splinetensor


    module subroutine print_info_splinetensor(fingerprint_parameters)
        implicit none
        !type(FingerPrint), intent(in) :: fp
        real(real64),intent(in)::fingerprint_parameters(:,:)
    
        double precision :: r_Rc, a_Rc
        integer          :: r_N, a_N
    
        r_Rc = fingerprint_parameters(1,1)
        r_N = nint(fingerprint_parameters(2,1))
        a_Rc = fingerprint_parameters(3,1)
        a_N = nint(fingerprint_parameters(4,1))
    
        write(*,*) 'Basis function type splinetensor'
        write(*,*) '[Y. Nagai and M. Okumura (2025)]'
        write(*,*)
        write(*,*) 'Radial Rc     : ', r_Rc
        write(*,*) 'Angular Rc    : ', a_Rc
        write(*,*) 'Radial points  : ', r_N
        write(*,*) 'Angular points : ', a_N
        write(*,*)
    
    end subroutine print_info_splinetensor


end submodule splinetensorbasis