module sub_predict
    use iso_fortran_env
    use inputdata,only:InputDataStruct
    use bfio,only:print_centered_header,FILEPATHLEN
    use structureinfo,only:Structuredata
    use mod_bppotential,only:BPpotential
    use fromaenet,only:geo_type_conv
    use lclist,   only: lcl_init,          &
                        lcl_final,         &
                        lcl_print_info,    &
                        lcl_nmax_nbdist,   &
                        lcl_nbdist_cart    
    use shared_potentials, only:initialize_potentials,bpnet_atomic_energy_and_forces_novirial,&
                        num_of_maxatoms_in_sphere,set_Rcmax_and_Rcmin,Rc_min,Rc_max,BPpots,&
                        deallocate_potentials,reload_potentials
    use fingerprints,only:FingerPrint,&
                        get_Rcmax_and_Rcmin,&
                        initialize_fingerprint_basisset 
    use fingerprint_basis,only:FingerPrintBasis                       

    implicit none

    contains

    subroutine predict_subroutine(inFile,ionum)
!        use module_sharing_potentials,only:BPpots
        implicit none
        character(len=*),intent(in)    :: inFile
        integer,intent(in)             :: ionum
        type(InputDataStruct)          :: inp
        integer::nStrucs
        integer::istruc
        character(len=FILEPATHLEN),allocatable::datafilenames(:)
        integer:: i
        character(len=256) :: line
        integer::u_in
        integer :: ios
        type(Structuredata)::xsfdata
        integer::iatom
        integer                                        :: num_of_neighbor_atoms!nnb
        real(real64),  dimension(:,:), allocatable :: coordinates_of_neighbor_atoms!nbcoo
        real(real64),  dimension(:),   allocatable :: distances_of_neighbor_atoms!nbdist
        integer,  dimension(:),   allocatable :: list_of_neighbor_atoms!nblist
        integer,           dimension(:),   allocatable :: types_of_neighbor_atoms!!nbtype
        !integer:: num_of_maxatoms_in_sphere
        !real(real64)::Rc_min,Rc_max
        integer::itype,itype1
        real(real64), dimension(:,:), allocatable :: forCart
        integer::stat
        real(real64)::E_i
        !type(FingerPrint),  dimension(:), allocatable :: fingerprint_set
        !class(FingerPrintBasis),dimension(:), allocatable  ::  fingerprint_basis_set
        integer:: num_max_coefficients
        real(real64)::Etot,Ecoh
        
        !type(BPpotential),allocatable::BPpots(:)

        inp = InputDataStruct()
        call inp%read_predict_in(inFile)


        call initialize_potentials(inp%num_of_Types,inp%typeName,inp%network_filenames)
        call reload_potentials(inp%num_of_Types,inp%typeName,inp%network_filenames)

        !write(*,*) "num_of_maxatoms_in_sphere ",num_of_maxatoms_in_sphere 
        !stop
        allocate(coordinates_of_neighbor_atoms(3,num_of_maxatoms_in_sphere ), &
        distances_of_neighbor_atoms(num_of_maxatoms_in_sphere ), &
        list_of_neighbor_atoms(num_of_maxatoms_in_sphere ),&
        types_of_neighbor_atoms(num_of_maxatoms_in_sphere ))
        !allocate(BPpots(1:inp%num_of_Types))

        !do i=1,inp%num_of_Types
        !    BPpots(i) = BPpotential(inp%typeName(i),inp%network_filenames(i))
        !end do
        

        if (inp%num_of_Strucs <= 0) then
            write(0,*) "Error: no input structures specified in ", &
                        trim(inFile)
            stop
        else
            nStrucs = inp%num_of_Strucs
        end if

        allocate(datafilenames(nStrucs))

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
                    datafilenames(i) = trim(adjustl(line))
                end do
                exit loop1
            end if
            if (ios /= 0) exit  ! EOFで終了
        end do loop1
        close(u_in)


        xsfdata = Structuredata() 

        

        do istruc = 1, nStrucs
            !write(*,*) trim(adjustl(datafilenames(istruc) ))
            call xsfdata%load(trim(adjustl(datafilenames(istruc))))
            !write(*,*) trim(adjustl(datafilenames(istruc) ))

            allocate(forCart(3,xsfdata%num_of_atoms))
            !write(*,*) "dd"
            !write(*,*) Rc_min, Rc_max, xsfdata%latticeVec
            !write(*,*)  xsfdata%num_of_atoms, xsfdata%atomType, xsfdata%abc_ith_atomic_coordinates, xsfdata%pbc


            call lcl_init(Rc_min, Rc_max, xsfdata%latticeVec, &
                    xsfdata%num_of_atoms, xsfdata%atomType, xsfdata%abc_ith_atomic_coordinates, xsfdata%pbc)

            !write(*,*) "atoms"
            call print_fileinfo(istruc, trim(adjustl(datafilenames(istruc))), &
                xsfdata%latticeVec, xsfdata%num_of_atoms, xsfdata%num_of_types)

            Ecoh = 0.0d0
            Etot = 0.0d0
            atoms: &
            do iatom = 1, xsfdata%num_of_atoms
                ! determine the training atom type of atom `iatom' in global
                ! index terms
                itype1 = xsfdata%atomType(iatom)
                !write(*,*) "atodd"
                itype1 = geo_type_conv(itype1, xsfdata%num_of_types, xsfdata%atomTypeName, &
                                inp%num_of_Types, inp%typeName)
                            !write(*,*) "atomsddd"
                ! assert that atom type is included in the set-ups:
                if (itype1 == 0) then
                    write(0,*) "Error: not a valid structure    : ", trim(xsfdata%structurefilename)
                    write(0,*) "       Additional species found."
                    stop
                end if

                ! get all atoms within cut-off:
                num_of_neighbor_atoms = num_of_maxatoms_in_sphere
                !write(*,*) "numa",num_of_maxatoms_in_sphere,coordinates_of_neighbor_atoms
                call lcl_nbdist_cart(iatom, num_of_neighbor_atoms, coordinates_of_neighbor_atoms, distances_of_neighbor_atoms, &
                                    r_cut=Rc_max, nblist=list_of_neighbor_atoms, nbtype=types_of_neighbor_atoms)
                !write(*,*) "numaa",num_of_neighbor_atoms


                call bpnet_atomic_energy_and_forces_novirial( &
                    xsfdata%calc_Cartesian_coordinates(iatom), &
                    itype1,&
                    iatom, &
                    num_of_neighbor_atoms, &
                    coordinates_of_neighbor_atoms, &
                    types_of_neighbor_atoms,&
                    list_of_neighbor_atoms, &
                    xsfdata%num_of_atoms, &
                    E_i, &
                    forCart, stat)
                !write(*,*) iatom,E_i
                Etot = Etot + E_i
                Ecoh = Ecoh + E_i - BPpots(itype1)%E_atom(itype1)
            end do atoms
            !write(*,*) trim(adjustl(datafilenames(istruc) ))
            write(*,*) "Etot",Etot
            call print_coordinates(istruc, &   
                    xsfdata%latticeVec,&
                    xsfdata%num_of_atoms, &
                    xsfdata%num_of_types, xsfdata%abc_ith_atomic_coordinates, &
                    xsfdata%atomType, xsfdata%atomTypeName, xsfdata%origin, forCart=forCart)
            call print_energy(Ecoh, Etot, forCart)

            call lcl_final()
            deallocate(forCart)
        end do

    end subroutine 


!subroutines from aenet 

  !--------------------------------------------------------------------!
  !                           general output                           !
  !--------------------------------------------------------------------!

    subroutine print_fileinfo(istruc, file, latticeVec, nAtoms, nTypes)

        implicit none
    
        integer,                                         intent(in) :: istruc
        character(len=*),                                intent(in) :: file
        double precision, dimension(3,3),                intent(in) :: latticeVec
        integer,                                         intent(in) :: nAtoms
        integer,                                         intent(in) :: nTypes
    
        write(*,*) 'Structure number  : ', istruc
        write(*,*) 'File name         : ', trim(adjustl(file))
        write(*,*) 'Number of atoms   : ', nAtoms
        write(*,*) 'Number of species : ', nTypes
        write(*,*)
    
        write(*,*) 'Lattice vectors (Angstrom):'
        write(*,*)
        write(*,'(3x,"a = ( ",3(2x,F15.8)," )")') latticeVec(1:3,1)
        write(*,'(3x,"b = ( ",3(2x,F15.8)," )")') latticeVec(1:3,2)
        write(*,'(3x,"c = ( ",3(2x,F15.8)," )")') latticeVec(1:3,3)
        write(*,*)
    
      end subroutine print_fileinfo

      subroutine print_energy(Ecoh, Etot, forCart)

        implicit none
    
        double precision,                           intent(in) :: Ecoh, Etot
        double precision, dimension(:,:), optional, intent(in) :: forCart
    
        double precision, dimension(3) :: F_mav, F_max, F_avg
        double precision               :: F_rms
        integer                        :: imax
    
        write(*,'(1x,"Cohesive energy            :",2x,F20.8," eV")') Ecoh
        write(*,'(1x,"Total energy               :",2x,F20.8," eV")') Etot
        if (present(forCart)) then
           call calc_rms_force(forCart, F_mav, F_max, imax, F_avg, F_rms)
           write(*,'(1x,"Mean force (must be zero)  :",3(2x,F12.6))') F_avg(1:3)
           write(*,'(1x,"Mean absolute force        :",3(2x,F12.6))') F_mav(1:3)
           write(*,'(1x,"Maximum force              :",3(2x,F12.6))') F_max(1:3)
           write(*,'(1x,"RMS force                  :",2x,F12.6)')    F_rms
           write(*,*) 'The maximum force is acting on atom ', &
                imax, '.'
           write(*,*) 'All forces are given in eV/Angstrom.'
        end if
        write(*,*)
    
      end subroutine print_energy

    subroutine print_coordinates(iter, latticeVec, nAtoms, nTypes, cooLatt, &
                               atomType, atomTypeName, origin, forCart, &
                               atomicEnergy)

    implicit none

    integer,                                         intent(in) :: iter
    double precision, dimension(3,3),                intent(in) :: latticeVec
    integer,                                         intent(in) :: nAtoms
    integer,                                         intent(in) :: nTypes
    double precision, dimension(3,nAtoms),           intent(in) :: cooLatt
    integer,          dimension(nAtoms),             intent(in) :: atomType
    character(len=*), dimension(nTypes),             intent(in) :: atomTypeName
    double precision, dimension(3),                  intent(in) :: origin
    double precision, dimension(3,nAtoms), optional, intent(in) :: forCart
    double precision, dimension(nAtoms),   optional, intent(in) :: atomicEnergy

    character(len=80)              :: header
    integer                        :: iat
    character(len=2)               :: symbol
    double precision, dimension(3) :: cooCart

        header = 'Cartesian atomic coordinates'
        if (iter == 0) then
        header = trim(header) // ' (input)'
        else
        header = trim(header) // ' (optimized)'
        end if
        if (present(forCart)) then
        header = trim(header) // ' and corresponding atomic forces'
        end if
        header = trim(header) // ':'

        write(*,*) trim(header)
        write(*,*)
        write(*,'(1x,2x,3(2x,A12))', advance='no') &
            '     x      ', '     y      ', '     z      '
        if (present(forCart)) then
        write(*,'(3(2x,A12))', advance='no') &
                '     Fx     ', '    Fy      ', '    Fz      '
        end if
        if (present(atomicEnergy)) then
        write(*,'(2x,A12)', advance='no') '   E_atom   '
        end if
        write(*,*)
        write(*,'(1x,2x,3(2x,A12))', advance='no') &
            '    (Ang)   ', '    (Ang)   ', '    (Ang)   '
        if (present(forCart)) then
        write(*,'(3(2x,A12))', advance='no') &
                '  (eV/Ang)  ', '  (eV/Ang)  ', '  (eV/Ang)  '
        end if
        if (present(atomicEnergy)) then
        write(*,'(2x,A12)', advance='no') '    (eV)    '
        end if
        write(*,*)
        write(*,'(1x,44("-"))', advance='no')
        if (present(forCart)) write(*,'(42("-"))', advance='no')
        if (present(atomicEnergy)) write(*,'(14("-"))', advance='no')
        write(*,*)
        do iat = 1, nAtoms
        symbol = atomTypeName(atomType(iat))
        cooCart(1:3) = matmul(latticeVec, cooLatt(1:3,iat)) + origin
        write(*,'(1x,A2,3(2x,F12.6))', advance='no') symbol, cooCart(1:3)
        if (present(forCart)) then
            write(*,'(3(2x,F12.6))', advance='no') forCart(1:3,iat)
        end if
        if (present(atomicEnergy)) then
            write(*,'(2x,F12.6)', advance='no') atomicEnergy(iat)
        end if
        write(*,*)
        end do
        write(*,*)

    end subroutine print_coordinates  
    
  !--------------------------------------------------------------------!
  !                       analyze atomic forces                        !
  !--------------------------------------------------------------------!

    subroutine calc_rms_force(forCart, F_mav, F_max, imax, F_avg, F_rms)

        implicit none
    
        double precision, dimension(:,:), optional, intent(in)  :: forCart
        double precision, dimension(3),             intent(out) :: F_mav
        double precision, dimension(3),             intent(out) :: F_max
        integer,                                    intent(out) :: imax
        double precision, dimension(3),             intent(out) :: F_avg
        double precision,                           intent(out) :: F_rms
    
        integer                        :: nAtoms
        double precision               :: F_abs2, F_abs2_max
        integer                        :: iat
    
        nAtoms = size(forCart(1,:))
        F_rms       = 0.0d0
        F_mav(1:3)  = 0.0d0
        F_max(1:3)  = 0.0d0
        F_avg(1:3)  = 0.0d0
        F_abs2      = 0.0d0
        F_abs2_max  = 0.0d0
        do iat = 1, nAtoms
           F_avg(1:3) = F_avg(1:3) + forCart(1:3,iat)
           F_mav(1:3) = F_mav(1:3) + abs(forCart(1:3,iat))
           F_abs2 = sum(forCart(1:3,iat)*forCart(1:3,iat))
           F_rms  = F_rms + F_abs2
           if (F_abs2 > F_abs2_max) then
              F_abs2_max  = F_abs2
              F_max(1:3) = forCart(1:3,iat)
              imax       = iat
           end if
        end do
        F_avg(1:3) = F_avg(1:3)/dble(nAtoms)
        F_mav(1:3) = F_mav(1:3)/dble(nAtoms)
        F_rms = sqrt(F_rms/dble(nAtoms))
    
      end subroutine calc_rms_force    

end module sub_predict