module mod_dataset
    use iso_fortran_env
    use bfio,only:FILEPATHLEN,print_centered_header
    use fingerprints,only:FingerPrint
    implicit none

    type,public::Dataset
        character(len=FILEPATHLEN)                  :: filename
        integer                                     :: num_of_Types
        integer                                     :: num_of_Strucs
        character(len=2), dimension(:), allocatable :: typeName
        real(real64),       dimension(:),   allocatable  :: isolated_atomic_energy
        logical                                           :: normalized
        real(real64)                                  :: scale, shift
        integer:: current_file_position !iStruc
        integer::file_iounit
        integer::num_totalatoms
        character(len=FILEPATHLEN),dimension(:),allocatable :: datafilenames
        real(real64)                                  :: E_min, E_max, E_av


        contains
        procedure::write_header => dataset_write_header
        procedure::write_structure_info => dataset_write_structure_info
        procedure::write_atom_info => dataset_write_atom_info
        procedure::write_fingerprint_values_info => dataset_write_fingerprint_values_info
        procedure::print_info => dataset_print_info
        procedure::write_footer => dataset_write_footer
        procedure::close => dataset_close
    end type Dataset

    interface Dataset
        module procedure init_Dataset
    end interface Dataset

    contains

    type(Dataset) function init_Dataset(nTypes, typeName, E_atom, nStrucs, file, scale, &
                shift) result(ts)
        implicit none

        integer,                             intent(in) :: nTypes
        character(len=*), dimension(nTypes), intent(in) :: typeName
        real(real64), dimension(nTypes), intent(in)          :: E_atom
        integer,                             intent(in) :: nStrucs
        character(len=*),                    intent(in) :: file
        real(real64), optional,                   intent(in) :: scale
        real(real64), optional,                   intent(in) :: shift

        logical :: fexists

        inquire(file=trim(adjustl(file)), exist=fexists)
        if (fexists) then
           write(0,*) 'Error: file already exists: ', trim(adjustl(file))
           stop
        end if

        allocate(ts%typeName(nTypes), ts%isolated_atomic_energy(nTypes))

        ts%filename               = trim(adjustl(file))
        ts%num_of_Types             = nTypes
        ts%typeName(1:nTypes) = typeName(1:nTypes)
        ts%isolated_atomic_energy(1:nTypes)   = E_atom(1:nTypes)
        ts%num_of_Strucs            = nStrucs
        ts%current_file_position             = 0

        if (present(scale) .and. present(shift)) then
            ts%normalized = .true.
            ts%scale = scale
            ts%shift = shift
        else
            ts%normalized = .false.
            ts%scale = 1.0d0
            ts%shift = 0.0d0
        end if

        open(newunit=ts%file_iounit, file=trim(ts%filename), status='new', action='write', &
                form='unformatted')

        ts%num_totalatoms = 0

        call ts%write_header()

    end function init_Dataset

    subroutine dataset_write_header(self)
        implicit none
        class(Dataset)::self

        write(self%file_iounit) self%num_of_Types
        write(self%file_iounit) self%num_of_Strucs
        write(self%file_iounit) self%typeName(:)
        write(self%file_iounit) self%isolated_atomic_energy(:)
        write(self%file_iounit) self%normalized
        write(self%file_iounit) self%scale
        write(self%file_iounit) self%shift

    end subroutine dataset_write_header

    subroutine dataset_close(self)
        implicit none
        class(Dataset)::self
        close(self%file_iounit)
    end subroutine 

    subroutine dataset_write_footer(self,fps)
        implicit none
        class(Dataset)::self
        type(FingerPrint), dimension(:), optional, intent(in)    :: fps
        integer :: itype, nTypes
        logical :: has_setups


        if (self%current_file_position < self%num_of_Strucs) then
            write(0,*) "Warning: writing footer to incomplete training set file."
        end if
        

        write(self%file_iounit) self%num_totalatoms
        write(self%file_iounit) self%E_av, self%E_min, self%E_max
        

        if (present(fps)) then
            nTypes = size(fps(:))
            if (nTypes /= self%num_of_Types) then
               write(0,*) "Error: wrong size of array fps in `dataset_write_footer()'."
               stop
            end if
            has_setups = .true.
            write(self%file_iounit) has_setups
            do itype = 1, self%num_of_Types!ts%nTypes
               write(self%file_iounit) itype
               
               call fps(itype)%write(unit=self%file_iounit)
               !call save_Setup(stp(itype), unit=ts%unit)
            end do
         else
            has_setups = .false.
            write(self%file_iounit) has_setups
         end if

    end subroutine dataset_write_footer

    subroutine dataset_write_structure_info(self,filename, nAtoms, nTypes, energy)
        implicit none        
        class(Dataset)::self
        character(len=*), intent(in)    :: filename
        integer,          intent(in)    :: nAtoms
        integer,          intent(in)    :: nTypes
        real(real64), intent(in)    :: energy
        real(real64) :: E_atom

        if (self%current_file_position >= self%num_of_Strucs) then
            write(0,*) "Error: too many files for training set."
            stop
        else
            self%current_file_position = self%current_file_position + 1
        end if

        write(self%file_iounit) len_trim(filename)
        write(self%file_iounit) trim(filename)
        write(self%file_iounit) nAtoms, nTypes
        write(self%file_iounit) energy
        
        ! energy stats
        E_atom = energy/dble(nAtoms)
        if (self%current_file_position > 1) then
            self%E_min = min(self%E_min, E_atom)
            self%E_max = max(self%E_max, E_atom)
            self%E_av  = self%E_av + E_atom/dble(self%num_of_Strucs)
        else
            self%E_min = E_atom
            self%E_max = E_atom
            self%E_av  = E_atom/dble(self%num_of_Strucs)
        end if

        ! keep track of the atoms in the training set
        self%num_totalatoms = self%num_totalatoms + nAtoms

    end subroutine dataset_write_structure_info

    subroutine dataset_write_atom_info(self,itype, cooCart, forCart)
        implicit none        
        class(Dataset)::self
        integer,                        intent(in)    :: itype
        real(real64), dimension(3), intent(in)    :: cooCart
        real(real64), dimension(3), intent(in)    :: forCart

        write(self%file_iounit) itype
        write(self%file_iounit) cooCart(1:3)
        write(self%file_iounit) forCart(1:3)
        
    end subroutine dataset_write_atom_info

    subroutine dataset_write_fingerprint_values_info(self, nsf, sfval)
        implicit none
        class(Dataset)::self
        integer,                            intent(in)    :: nsf
        real(real64), dimension(nsf),   intent(in)    :: sfval
    
        write(self%file_iounit) nsf
        write(self%file_iounit) sfval(1:nsf)
    
      end subroutine dataset_write_fingerprint_values_info

      subroutine dataset_print_info(self)
        implicit none
        class(Dataset)::self
        integer::itype

        call print_centered_header("Training set info.")
        write(*,*)

        write(*,*) 'Training set file                   : ', trim(adjustl(self%filename))
        write(*,*) 'Number of structures in the data set: ', self%num_of_Strucs
    
        write(*,*) 'Atomic species in training set      : ', self%num_of_Types
        write(*,'(1x,"  Species :")', advance='no')
        do itype = 1, self%num_of_Types
           write(*,'(1x,A)', advance='no') trim(self%typeName(itype))
        end do
        write(*,*)
        write(*,*)

        if (self%normalized .or. (self%current_file_position == self%num_of_Strucs)) then
            write(*,*) 'Average energy (eV/atom) : ', self%E_av
            write(*,*) 'Minimum energy (eV/atom) : ', self%E_min
            write(*,*) 'Maximum energy (eV/atom) : ', self%E_max
            write(*,*)
         end if
     
        if (self%normalized) then
            write(*,*) 'The input and output values have been normalized to [-1.0, 1.0].'
            write(*,*) 'Structures outside of this interval will not be used for training.'
            write(*,*) '  Energy scaling factor: ', self%scale
            write(*,*) '  Atomic energy shift  : ', self%shift
        else
            write(*,*) 'The input and output values have not yet been normalized.'
        end if
        write(*,*)

      end subroutine dataset_print_info
end module mod_dataset
