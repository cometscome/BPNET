module sub_generate
    use iso_fortran_env
    use inputdata,only:InputDataStruct
    use fingerprints,only:FingerPrint,&
                get_Rcmax_and_Rcmin,&
                initialize_fingerprint_basisset
    use neighborlist,only:max_packed_spheres
    use fingerprint_basis,only:FingerPrintBasis
    use bfio,only:print_centered_header,FILEPATHLEN
    use mod_dataset,only:dataset
    use structureinfo,only:Structuredata
    use fromaenet,only:geo_type_conv
    use lclist,   only: lcl_init,          &
                        lcl_final,         &
                        lcl_print_info,    &
                        lcl_nmax_nbdist,   &
                        lcl_nbdist_cart
    use bpnet_trainbin2ascii,only:trainbin2ascii_subroutine
    !use chebyshevbasis_fp,only:FP_ChebyshevBasis
    implicit none
    public::generate_subroutine


    contains
    


    subroutine generate_subroutine(inFile,ionum)
        implicit none
        character(len=*),intent(in)    :: inFile
        integer,intent(in)             :: ionum
        type(InputDataStruct)          :: inp
        type(FingerPrint),  dimension(:), allocatable :: fingerprint_set
        class(FingerPrintBasis),dimension(:), allocatable  ::  fingerprint_basis_set
        real(real64), dimension(:),     allocatable :: fingerprint_values
        real(real64), dimension(:,:),   allocatable :: fingerprint_deriv_i
        real(real64), dimension(:,:,:), allocatable :: fingerprint_deriv_j

        integer:: num_of_maxatoms_in_sphere
        integer:: num_max_coefficients
        real(real64)::Rc_min,Rc_max
        integer::itype,itype1
        type(Dataset) :: traindataset
        integer::u_in
        integer :: ios
        character(len=256) :: line
        !character(len=1024)                         :: keyword
        !character(len=FILEPATHLEN),dimension(:),allocatable :: datafilenames
        integer::i
        integer::istruct
        type(Structuredata)::xsfdata
        real(real64)::cohesiveEnergy
        integer::iatom
        integer                                        :: num_of_neighbor_atoms!nnb
        real(real64),  dimension(:,:), allocatable :: coordinates_of_neighbor_atoms!nbcoo
        real(real64),  dimension(:),   allocatable :: distances_of_neighbor_atoms!nbdist
        integer,  dimension(:),   allocatable :: list_of_neighbor_atoms!nblist
        integer,           dimension(:),   allocatable :: types_of_neighbor_atoms!!nbtype
        character(len=1024)::outfilename_ascii
        character(len=1024)::outfilename
        logical ::to_bin,to_ascii

        inp = InputDataStruct()
        call inp%read_generate_in(inFile)

        allocate(fingerprint_set(inp%num_of_Types))
        call load_fingerprint_setups(inp,fingerprint_set)

        call get_Rcmax_and_Rcmin(fingerprint_set,inp%num_of_Types,Rc_min, Rc_max)
        !num_of_maxatoms_in_sphere = max_packed_spheres(Rc_min,Rc_max)
        num_of_maxatoms_in_sphere = lcl_nmax_nbdist(Rc_min, Rc_max)
        allocate(coordinates_of_neighbor_atoms(3,num_of_maxatoms_in_sphere ), &
                distances_of_neighbor_atoms(num_of_maxatoms_in_sphere ), &
                list_of_neighbor_atoms(num_of_maxatoms_in_sphere ),&
                types_of_neighbor_atoms(num_of_maxatoms_in_sphere ))
  
        call initialize_fingerprint_basisset(fingerprint_set,inp%num_of_types,&
                fingerprint_basis_set,&
                num_max_coefficients)

        allocate(fingerprint_values(num_max_coefficients), source=0d0)
        allocate(fingerprint_deriv_i(3,num_max_coefficients), source=0d0)
        allocate(fingerprint_deriv_j(3,num_max_coefficients,num_of_maxatoms_in_sphere), source=0d0)

        call print_centered_header("Generation of training set started")
        !write(*,*) fingerprint_basis_set(1)%get_num_parameters()

        write(ionum,*)
  
        write(ionum,*) 'Number of atom types  : ',inp%num_of_Types
        write(ionum,'(1x,"types                 : ")', advance='no')
        do itype = 1, inp%num_of_Types
            if (mod(itype,7) == 0) write(*,'(29x)')
            write(ionum,'(A5,1x)', advance='no') inp%typeName(itype)
        end do
        write(ionum,*)
        write(ionum,*) "Number of structures  : ", inp%num_of_Strucs
        write(ionum,*)

        call print_centered_header("Structural fingerprint basis set-up")
        write(ionum,*)

        do itype1 = 1, inp%num_of_Types
            call fingerprint_set(itype1)%print_info()
        end do

        !write(*,*) trim(inp%trainingdataset_filename)
        traindataset = Dataset(inp%num_of_Types, inp%typeName, inp%isolated_atomic_energy, &
                inp%num_of_Strucs, trim(inp%trainingdataset_filename))

        outfilename_ascii = trim(inp%trainingdataset_filename)//'.ascii'   

        allocate(traindataset%datafilenames(inp%num_of_Strucs))
        call print_centered_header("Adding structures to the training set")
        write(ionum,*)

        open(newunit=u_in, file=trim(adjustl(inFile)), status='old', action='read')
        rewind(u_in)

        loop1:&
        do
            read(u_in, '(A)', iostat=ios) line
            !write(*,*) line
            if (trim(line) == 'FILES') then
                !write(*,*) "found"
                read(u_in, '(A)', iostat=ios) line
                do i=1,inp%num_of_Strucs
                    read(u_in, '(A)', iostat=ios) line
                    !write(*,*) line
                    traindataset%datafilenames(i) = trim(adjustl(line))
                end do
                exit loop1
            end if
            if (ios /= 0) exit  ! EOFで終了
        end do loop1
        close(u_in)


        ! header for stdout
        write(*,'("#",A6,2x,A6,2x,A6,2x,A15,2x,A)') &
        'N', 'nAtoms', 'nTypes', 'E/atom', 'structure file (xsf)'

        xsfdata = Structuredata() 

        structures: &
        do istruct=1,inp%num_of_Strucs
            !write(*,*) istruct,trim(adjustl(traindataset%datafilenames(istruct)))
            call xsfdata%load(trim(adjustl(traindataset%datafilenames(istruct))))

            if (.not. (xsfdata%hasForces .and. xsfdata%hasEnergy)) then
                write(0,*) ">>>", xsfdata%hasForces, xsfdata%hasEnergy
                write(0,*) "Error: incomplete output data in : ", trim(xsfdata%structurefilename)
                stop
            end if

            if (xsfdata%num_of_types > inp%num_of_Types) then
                write(ionum,*) 'Skipping ', trim(adjustl(xsfdata%structurefilename)), &
                        ': too many atomic species'
                cycle structures
            end if

            if (abs(xsfdata%cohesiveEnergy) /= 0.0d0) then
                cohesiveEnergy = xsfdata%cohesiveEnergy
            else
                ! if only the total energy is available, we have to calculate
                ! the cohesive energy at this point
                cohesiveEnergy = xsfdata%totalEnergy
                do iatom = 1, xsfdata%num_of_atoms
                   itype1 = xsfdata%atomType(iatom)
                   itype1 = geo_type_conv(itype1, xsfdata%num_of_types, xsfdata%atomTypeName, &
                                    inp%num_of_Types, inp%typeName)
                   cohesiveEnergy = cohesiveEnergy - inp%isolated_atomic_energy(itype1)
                end do
            end if

            write(*,'(1x,I6,2x,I6,2x,I6,2x,ES15.8,2x,A)') &
                    istruct, xsfdata%num_of_atoms, xsfdata%num_of_types, cohesiveEnergy/dble(xsfdata%num_of_atoms), &
                    trim(adjustl(xsfdata%structurefilename))

            call lcl_init(Rc_min, Rc_max, xsfdata%latticeVec, &
                    xsfdata%num_of_atoms, xsfdata%atomType, xsfdata%abc_ith_atomic_coordinates, xsfdata%pbc)

            call traindataset%write_structure_info(xsfdata%structurefilename, &
                        xsfdata%num_of_atoms, xsfdata%num_of_types, cohesiveEnergy)

            atoms: &
            do iatom = 1, xsfdata%num_of_atoms
                ! determine the training atom type of atom `iatom' in global
                ! index terms
                itype1 = xsfdata%atomType(iatom)
                itype1 = geo_type_conv(itype1, xsfdata%num_of_types, xsfdata%atomTypeName, &
                                inp%num_of_Types, inp%typeName)
                ! assert that atom type is included in the set-ups:
                if (itype1 == 0) then
                    write(0,*) "Error: not a valid structure    : ", trim(xsfdata%structurefilename)
                    write(0,*) "       Additional species found."
                    stop
                end if

                ! write atom info (species, forces) to training set file:
                call traindataset%write_atom_info(itype1, &
                    xsfdata%calc_Cartesian_coordinates(iatom), &
                    xsfdata%calc_Cartesian_forces(iatom))
  
                 ! get all atoms within cut-off:
                num_of_neighbor_atoms = num_of_maxatoms_in_sphere
                call lcl_nbdist_cart(iatom, num_of_neighbor_atoms, coordinates_of_neighbor_atoms, distances_of_neighbor_atoms, &
                                    r_cut=Rc_max, nblist=list_of_neighbor_atoms, nbtype=types_of_neighbor_atoms)
  
                ! convert atom types to global index:
                do i = 1, num_of_neighbor_atoms
                    types_of_neighbor_atoms(i) = geo_type_conv(types_of_neighbor_atoms(i), &
                                            xsfdata%num_of_types,xsfdata%atomTypeName, &
                                            inp%num_of_Types, inp%typeName)
                    if (types_of_neighbor_atoms(i) == 0) then
                        write(0,*) "Error: atom type not found in setup."
                        stop
                    end if
                end do

                ! evaluate the structural fingerprint basis function set-up:
                call fingerprint_basis_set(itype1)%evaluate(itype1,&
                                            xsfdata%calc_Cartesian_coordinates(iatom), &
                                            num_of_neighbor_atoms, &
                                            coordinates_of_neighbor_atoms, &
                                            types_of_neighbor_atoms, &
                                            fingerprint_set(itype1)%ltype,&
                                            values=fingerprint_values, &
                                            deriv_i=fingerprint_deriv_i, &
                                            deriv_j=fingerprint_deriv_j)

                
                call fingerprint_set(itype1)%update_statistics(fingerprint_values)


                call fingerprint_set(itype1)%count_evaluation()

                ! write basis function values and derivatives
                ! to the training set file:
                call traindataset%write_fingerprint_values_info(&
                        fingerprint_set(itype1)%num_of_coeffs, &
                        fingerprint_values(1:fingerprint_set(itype1)%num_of_coeffs))
                    
            end do atoms

            call lcl_final()
        end do structures

        call traindataset%print_info()

        to_bin = .false.
        to_ascii = .true.
        outFileName = trim(adjustl(inp%trainingdataset_filename))

        call traindataset%write_footer(fps=fingerprint_set(1:inp%num_of_Types))
        
        !
        !call traindataset%write_footer()
        call traindataset%close()

        call print_centered_header("Training set generation done.")

        call trainbin2ascii_subroutine(trim(outFileName), trim(outfilename_ascii),to_bin, to_ascii)

    end subroutine

    subroutine load_fingerprint_setups(inp,fingerprint_set)
        implicit none
  
        type(InputDataStruct),           intent(in)  :: inp
        type(FingerPrint), dimension(:), intent(out) :: fingerprint_set

        integer :: i

        do i = 1, inp%num_of_Types
            fingerprint_set(i) = FingerPrint()
            call fingerprint_set(i)%read(inp%fingerprint_setup_file(i), inp%typeName)
        end do


    end subroutine load_fingerprint_setups

end module