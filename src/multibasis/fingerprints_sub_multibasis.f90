submodule (fingerprints) multibasis
    use iso_fortran_env
    implicit none
    contains

    module subroutine read_basis_multi(fp,unit,num_kinds)
        implicit none
        type(FingerPrint),intent(inout)::fp
        integer,intent(in)::num_kinds
        integer,intent(in)::unit
        integer,           parameter :: nc = 1024
        character(len=nc)            :: line
        character(len=nc)            :: tline
        !integer             :: r_N, a_N
        !real(real64)    :: r_Rc, a_Rc
        integer::nTypes,ikind,num_of_coeffs_i,num_of_coeffs_total
        logical :: eof
        real(real64)::paramdata(5)
        integer::num_of_maxparameters_total,ibasis
        real(real64),allocatable::paramdata_total(:)
        real(real64)::Rc
        integer::istart
        fp%Rc_max = 0d0

        num_of_maxparameters_total = num_kinds*(num_of_maxparameters+2)+2
        allocate(paramdata_total(num_of_maxparameters_total ))
        num_of_coeffs_total = 0
        paramdata_total(1) = 0d0
        paramdata_total(2) = dble(num_kinds)
                
        do ikind=1,num_kinds
            call read_next_valid_line(unit, tline, eof)
            select case(to_lowercase(tline))
            case('chebyshev')
                ibasis = 1
                call read_next_valid_line(unit, line, eof)
                call read_chebyshevdata(line,paramdata,num_of_coeffs_i ,Rc)
                
                if (fp%num_of_different_species > 1) then
                    num_of_coeffs_i = 2*num_of_coeffs_i
                end if
                fp%Rc_max = max(fp%Rc_max, Rc)
                
            case('spline')
                ibasis = 2
                call read_next_valid_line(unit, line, eof)    
                call read_splinedata(line,paramdata,num_of_coeffs_i ,Rc)
                if (fp%num_of_different_species > 1) then
                    num_of_coeffs_i = 2*num_of_coeffs_i
                end if
                fp%Rc_max = max(fp%Rc_max, Rc)
                !call read_basis_spline(fp,line)
            case('lj')
                ibasis = 3
                call read_next_valid_line(unit, line, eof)   
                call read_LJdata(line,fp%num_of_different_species,paramdata,num_of_coeffs_i,Rc ) 
                fp%Rc_max = max(fp%Rc_max, Rc)
                !call read_basis_LJ(fp,line)
            case('multi')
                write(*,*) "multi is not supported in multi basis"
                stop
                !call read_next_valid_line(unit, line, eof)    
                
                !call read_basis_LJ(fp,line)
            case default
                write(0,*) "read_basis_multi: Error: Unknown basis type: ", trim(tline)
                !deallocate(stp%sf, stp%sfparam)
                stop
            end select
            num_of_coeffs_total  = num_of_coeffs_total + num_of_coeffs_i 
            istart = (ikind-1)*(num_of_maxparameters+2) +3
            paramdata_total(istart) = dble(ibasis)
            paramdata_total(istart+1) = dble(num_of_coeffs_i)
            paramdata_total(istart+2:istart+num_of_maxparameters+1) =&
                    paramdata(1:num_of_maxparameters)
        end do


        fp%num_of_coeffs = num_of_coeffs_total 

        call fp%allocate_statistics_arrays(fp%num_of_coeffs)
        fp%num_of_parameters = num_of_maxparameters_total 
        allocate(fp%fingerprint_parameters(num_of_maxparameters_total ,fp%num_of_coeffs),source=0d0)

        allocate(fp%fingerprint_function_kind(fp%num_of_coeffs),source=0)
        allocate(fp%species_of_basis_functions(NENV_MAX,fp%num_of_coeffs),source=0)

        
        fp%fingerprint_parameters(:,1) = paramdata_total(:)


    end subroutine read_basis_multi


    module subroutine print_info_multi(fingerprint_parameters)
        implicit none
        real(real64), intent(in) :: fingerprint_parameters(:,:)
        integer::num_kinds
        double precision :: r_Rc, a_Rc
        integer::ikind,ibasis
        !integer          :: r_N, a_N

        write(*,*) 'Basis function type multi'
        write(*,*) '[Y. Nagai and M. Okumura (2025)]'
        write(*,*)
        num_kinds =  nint(fingerprint_parameters(1,1))
        write(*,*) 'num. of basis is ',num_kinds
        do ikind = 1,num_kinds
            write(*,*) "------------------------------------------"
            write(*,*) ikind,"-th basis: "
            ibasis =  nint(fingerprint_parameters(2+(ikind-1)*(num_of_maxparameters+1),1) )
            !write(*,*) fingerprint_parameters(3+(ikind-1)*(num_of_maxparameters+1):num_of_maxparameters+1+(ikind-1)*(num_of_maxparameters+1),:)
            if (ibasis == 1) then !chebyshev
                call print_info_chebyshev(fingerprint_parameters(3+(ikind-1)* &
                    (num_of_maxparameters+1):num_of_maxparameters+1+(ikind-1)*(num_of_maxparameters+1),:))
            elseif(ibasis == 2) then !spline
                call print_info_spline(fingerprint_parameters(3+(ikind-1)* &
                    (num_of_maxparameters+1):num_of_maxparameters+1+(ikind-1)*(num_of_maxparameters+1),:))
            elseif(ibasis == 3) then !LJ
                call print_info_LJ(fingerprint_parameters(3+(ikind-1)* &
                    (num_of_maxparameters+1):num_of_maxparameters+1+(ikind-1)*(num_of_maxparameters+1),:))
            end if
            write(*,*) "------------------------------------------"
        end do

    end subroutine print_info_multi


end submodule multibasis