

program nnASCII2bin
    use bpnet_nnASCII2bin
    implicit none

    character(len=1024) :: infile, outfile
    logical             :: to_bin, to_ascii

    call initialize(infile, outfile, to_bin, to_ascii)
    call nnASCII2bin_subroutine(infile, outfile, to_bin, to_ascii)

    contains 

    subroutine initialize(infile, outfile, to_bin, to_ascii)

        implicit none
    
        character(len=*), intent(out) :: infile, outfile
        integer :: iarg, nargs
        character(len=100) :: arg
        logical, intent(out) :: to_bin, to_ascii
    
        nargs = command_argument_count()
        if (nargs < 1) then
           write(0,*) "Error: No input file provided."
           call finalize()
           stop
        end if
    
        infile = ' '
        outfile = ' '
    
        to_bin = .true.
        to_ascii = .false.
        
        iarg = 1
        do while(iarg <= nargs)
           call get_command_argument(iarg, value=arg)
           select case(trim(arg))
           case('--to-binary')
              to_bin = .false.
              to_ascii = .true.
           case('--to-ascii')
              to_ascii = .false.
              to_bin = .true.
           case default
              if (len_trim(infile) == 0) then
                 infile = trim(arg)
              else if (len_trim(outfile) == 0) then
                 outfile = trim(arg)
              else
                 write(0,*) 'Error: Unknown argument: ', trim(arg)
                 call finalize()
                 stop
              end if
           end select
           iarg = iarg + 1
        end do
    
        if ((len(infile) == 0) .or. (len(outfile) == 0))then
           write(0,*) 'Error: No input file specified.'
           call finalize()
           stop
        end if
    
    end subroutine initialize

    subroutine finalize()

        implicit none
    
    end subroutine finalize
end program