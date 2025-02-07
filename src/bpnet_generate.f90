program bpnet_generate
    use iso_fortran_env
    use bfio,only:FILEPATHLEN,print_centered_header,&
            print_mit_license
    use ioutils,only:file_exists
    use sub_generate,only:generate_subroutine

    implicit none
    integer::ionum
    character(len=FILEPATHLEN)   :: inFile

    
    ionum = 139
    call initialize(inFile)
    call generate_subroutine(inFile,ionum)
    !call print_centered_header("test")

    contains


    subroutine initialize(inFile)
        implicit none
        character(len=*), intent(out) :: inFile
        integer :: nargs
        logical :: fexists

        call print_centered_header("bpnet_generate.x - training set generation")!, char='=')
        call print_mit_license()


        nargs = command_argument_count()
        if (nargs < 1) then
            write(0,*) "Error: No input file provided."
            call print_usage()
            stop
        end if
    
        call get_command_argument(1, value=inFile)
        inquire(file=trim(inFile), exist=fexists)
        if (.not. fexists) then
            write(0,*) "Error: File not found: ", trim(inFile)
            call print_usage()
            stop
        end if

    end subroutine initialize

    subroutine print_usage()

        implicit none

        write(ionum,*)
        write(ionum,*) "bpnet_generate.x -- Generate training sets for use with `train.x'"
        write(ionum,'(1x,70("-"))')
        write(ionum,*) 'Usage: bpnet_generate.x <input-file>'
        write(ionum,*)
        write(ionum,*) 'See the documentation or the source code for a description of the '
        write(ionum,*) 'input file format.'
        write(ionum,*)

    end subroutine print_usage


end program