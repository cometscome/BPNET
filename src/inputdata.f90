module inputdata
    use iso_fortran_env
    use bfio,only:FILEPATHLEN
    use ioutils,only:file_exists, find_keyword_in_file

    implicit none

    type, public :: InputDataStruct
        integer                                             :: num_of_Types
        integer                                             :: num_of_Strucs
        character(len=FILEPATHLEN)    :: input_filename
        character(len=2), dimension(:),   allocatable :: typeName
        real(real64),       dimension(:),   allocatable :: isolated_atomic_energy
        character(len=FILEPATHLEN), dimension(:),   allocatable :: fingerprint_setup_file
        character(len=FILEPATHLEN)                              :: trainingdataset_filename

        character(len=FILEPATHLEN), dimension(:),   allocatable :: network_filenames
        contains
        procedure:: read_generate_in
        procedure:: read_predict_in
    end type InputDataStruct

    interface InputDataStruct
        module procedure init_InputDataStruct
    end interface InputDataStruct

    contains

        !   コンストラクタ
    type(InputDataStruct) function init_InputDataStruct() result(inp)
        implicit none
        inp%num_of_Types = 0
        inp%num_of_Strucs = 0
    end function init_InputDataStruct

    subroutine read_predict_in(self,file)
        class(InputDataStruct)::self
        character(len=*), intent(in) :: file
        logical::exists
        integer::u

        call  file_exists(file,exists)
        if (exists .eqv. .false.) then
            write(*,*) "file ",file,"is not found"
        endif
        open(newunit=u, file=trim(adjustl(file)), status='old', action='read')
        
        self%input_filename = file

        call read_types_section_forpredict(u, self%num_of_Types,self%typeName)
        close(u)

        open(newunit=u, file=trim(adjustl(file)), status='old', action='read')
        allocate(self%network_filenames(1:self%num_of_Types))
        call read_networks_section(u,self%num_of_Types,self%typeName,self%network_filenames)
        close(u)
        !write(*,*) self%network_filenames

        open(newunit=u, file=trim(adjustl(file)), status='old', action='read')
        call read_files_section(u,self%num_of_Strucs)
        close(u)
        !write(*,*) self%num_of_Strucs

        return

    end subroutine 

    subroutine read_generate_in(self,file)
        class(InputDataStruct)::self
        character(len=*), intent(in) :: file
        logical::exists
        integer::u

       
        call  file_exists(file,exists)
        if (exists .eqv. .false.) then
            write(*,*) "file ",file,"is not found"
        endif
        open(newunit=u, file=trim(adjustl(file)), status='old', action='read')
        
        self%input_filename = file

        call read_types_section(u, self%num_of_Types,self%typeName,self%isolated_atomic_energy)
        close(u)
        open(newunit=u, file=trim(adjustl(file)), status='old', action='read')
        allocate(self%fingerprint_setup_file(1:self%num_of_Types))
        call read_setups_section(u, self%num_of_Types,self%typeName,self%fingerprint_setup_file)
        close(u)

        open(newunit=u, file=trim(adjustl(file)), status='old', action='read')
        call read_files_section(u,self%num_of_Strucs)
        close(u)

        !write(*,*) self%num_of_Strucs
        !write(*,* )self%typeName
        !write(*,*) self%isolated_atomic_energy
        !write(*,*) self%fingerprint_setup_file

        call find_keyword_in_file(file, 'output', self%trainingdataset_filename)

        return
    end subroutine read_generate_in

    subroutine read_types_section_forpredict(unit, num_types,typeName)
        implicit none
        integer, intent(in) :: unit           ! ファイルユニット番号
        integer, intent(out) :: num_types     ! TYPESの数
        character(len=2), allocatable, intent(out) :: typeName(:)
        !character(len=256), allocatable, intent(out) :: types_data(:) ! TYPESデータ配列

        character(len=256) :: line
        integer :: ios, i
        logical :: found_types

        found_types = .false.

        ! ファイルを1行ずつ読み込む
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! EOFで終了
            !write(*,*) line

            ! typesを見つける
            if (trim(line) == 'TYPES') then
                found_types = .true.
                read(unit, '(A)', iostat=ios) line  ! 次の行を読む（元素数）
                if (ios /= 0) then
                    print *, 'Error: Unexpected EOF after types'
                    stop
                end if
                read(line, *) num_types  ! 元素数を取得
                !print *, 'Number of types:', num_types

                ! 配列を確保
                allocate(typeName(num_types))

                ! typesの次のnum_types行を読み込む
                do i = 1, num_types
                    read(unit, '(A)', iostat=ios) line
                    if (ios /= 0) then
                        print *, 'Error: Unexpected EOF while reading TYPES section'
                        stop
                    end if
                    read(line,*) typeName(i)
                end do

                return  ! typesセクションを読み込んだらサブルーチンを終了
            end if
        end do

        if (.not. found_types) then
            print *, 'Error: TYPES section not found'
            num_types = 0
        end if
    end subroutine read_types_section_forpredict  

    subroutine read_networks_section(unit, num_types,typeName_in,network_filename)
        implicit none
        integer, intent(in) :: unit           ! ファイルユニット番号
        integer, intent(in) :: num_types     ! TYPESの数
        !real(real64),allocatable, intent(out) :: isolated_atomic_energy(:)
        !character(len=256), allocatable, intent(out) :: types_data(:) ! TYPESデータ配列
        character(len=FILEPATHLEN),intent(out) :: network_filename(1:num_types)
        character(len=2),intent(in):: typeName_in(1:num_types)

        character(len=256) :: line
        integer :: ios, i
        logical :: found_setups
        character(len=2):: typeName(1:num_types)


        found_setups = .false.

        ! ファイルを1行ずつ読み込む
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! EOFで終了

            ! networksを見つける
            if (trim(line) == 'NETWORKS') then
                found_setups = .true.
            
                ! networksの次のnum_types行を読み込む
                do i = 1, num_types
                    read(unit, '(A)', iostat=ios) line
                    if (ios /= 0) then
                        print *, 'Error: Unexpected EOF while reading networks section'
                        stop
                    end if
                    read(line,*) typeName(i), network_filename(i)
                    if (trim(typeName(i)) .ne. trim(typeName_in(i))) then
                        write(*,*) "type in TYPES section ",trim(typeName_in(i))
                        write(*,*) "type in NETWORKS section ",trim(typeName(i))
                        print *, 'Error: type name mismatch! please check the input file for generate.x '
                        stop
                    endif
                end do

                return  ! SETUPSセクションを読み込んだらサブルーチンを終了
            end if
        end do

        if (.not. found_setups) then
            print *, 'Error: NETWORKS section not found'
            stop
        end if
    end subroutine read_networks_section 


    subroutine read_types_section(unit, num_types,typeName,isolated_atomic_energy)
        implicit none
        integer, intent(in) :: unit           ! ファイルユニット番号
        integer, intent(out) :: num_types     ! TYPESの数
        character(len=2), allocatable, intent(out) :: typeName(:)
        real(real64),allocatable, intent(out) :: isolated_atomic_energy(:)
        !character(len=256), allocatable, intent(out) :: types_data(:) ! TYPESデータ配列

        character(len=256) :: line
        integer :: ios, i
        logical :: found_types

        found_types = .false.

        ! ファイルを1行ずつ読み込む
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! EOFで終了

            ! TYPESを見つける
            if (trim(line) == 'TYPES') then
                found_types = .true.
                read(unit, '(A)', iostat=ios) line  ! 次の行を読む（元素数）
                if (ios /= 0) then
                    print *, 'Error: Unexpected EOF after TYPES'
                    stop
                end if
                read(line, *) num_types  ! 元素数を取得
                !print *, 'Number of types:', num_types

                ! 配列を確保
                allocate(typeName(num_types))
                allocate(isolated_atomic_energy(num_types))

                ! TYPESの次のnum_types行を読み込む
                do i = 1, num_types
                    read(unit, '(A)', iostat=ios) line
                    if (ios /= 0) then
                        print *, 'Error: Unexpected EOF while reading TYPES section'
                        stop
                    end if
                    read(line,*) typeName(i), isolated_atomic_energy(i)
                end do

                return  ! TYPESセクションを読み込んだらサブルーチンを終了
            end if
        end do

        if (.not. found_types) then
            print *, 'Error: TYPES section not found'
            num_types = 0
        end if
    end subroutine read_types_section  

    subroutine read_setups_section(unit, num_types,typeName_in,fingerprint_setup_file)
        implicit none
        integer, intent(in) :: unit           ! ファイルユニット番号
        integer, intent(in) :: num_types     ! TYPESの数
        character(len=FILEPATHLEN),intent(out) :: fingerprint_setup_file(1:num_types)
        !real(real64),allocatable, intent(out) :: isolated_atomic_energy(:)
        !character(len=256), allocatable, intent(out) :: types_data(:) ! TYPESデータ配列
        character(len=2),intent(in):: typeName_in(1:num_types)
        character(len=256) :: line
        integer :: ios, i
        logical :: found_setups
        character(len=2):: typeName(1:num_types)


        found_setups = .false.

        ! ファイルを1行ずつ読み込む
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! EOFで終了

            ! TYPESを見つける
            if (trim(line) == 'SETUPS') then
                found_setups = .true.
            
                ! SETUPSの次のnum_types行を読み込む
                do i = 1, num_types
                    read(unit, '(A)', iostat=ios) line
                    if (ios /= 0) then
                        print *, 'Error: Unexpected EOF while reading SETUPS section'
                        stop
                    end if
                    read(line,*) typeName(i), fingerprint_setup_file(i)
                    if (trim(typeName(i)) .ne. trim(typeName_in(i))) then
                        write(*,*) "type in TYPES section ",trim(typeName_in(i))
                        write(*,*) "type in SETUPS section ",trim(typeName(i))
                        print *, 'Error: type name mismatch! please check the input file for generate.x '
                        stop
                    endif
                end do

                return  ! SETUPSセクションを読み込んだらサブルーチンを終了
            end if
        end do

        if (.not. found_setups) then
            print *, 'Error: SETUPS section not found'
            stop
        end if
    end subroutine read_setups_section 


    subroutine read_files_section(unit,numfile)
        implicit none
        integer, intent(in) :: unit           ! ファイルユニット番号
        integer,intent(out) ::numfile
        !integer, intent(in) :: num_types     ! TYPESの数
        !character(len=FILEPATHLEN),intent(out) :: fingerprint_setup_file(1:num_types)
        !real(real64),allocatable, intent(out) :: isolated_atomic_energy(:)
        !character(len=256), allocatable, intent(out) :: types_data(:) ! TYPESデータ配列
        !character(len=2),intent(in):: typeName_in(1:num_types)
        character(len=256) :: line
        integer :: ios, i
        logical :: found_files
        !character(len=2):: typeName(1:num_types)


        found_files = .false.

        ! ファイルを1行ずつ読み込む
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! EOFで終了

            ! TYPESを見つける
            if (trim(line) == 'FILES') then
                found_files = .true.
                read(unit, '(A)', iostat=ios) line
                read(line,*) numfile
            
                return  ! FILESセクションを読み込んだらサブルーチンを終了
            end if
        end do

        if (.not. found_files) then
            print *, 'Error: FILES section not found'
            stop
        end if
    end subroutine read_files_section     
end module inputdata