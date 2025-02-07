module structureinfo
    use iso_fortran_env
    use bfio,only:FILEPATHLEN
    implicit none
    real(real64), parameter, private:: PI = 3.141592653589793d0

    type,public::Structuredata
        logical,                                    public :: pbc
        real(real64),       dimension(3,3),              public :: latticeVec
        real(real64),       dimension(3,3),              public :: recLattVec
        real(real64),       dimension(3),                public :: origin

        integer,                                    public :: num_of_atoms
        integer,                                    public :: num_of_types
        logical,                                    public :: hasForces

        real(real64),       dimension(:,:), allocatable, public :: abc_ith_atomic_coordinates !cooLatt
        real(real64),       dimension(:,:), allocatable, public :: abc_ith_atomic_forces !forLatt
        character(len=2), dimension(:),   allocatable, public :: atomTypeName
        integer,          dimension(:),   allocatable, public :: atomType

        logical,                                    public :: hasEnergy
        real(real64),                                    public :: cohesiveEnergy
        real(real64),                                    public :: totalEnergy
        character(len=FILEPATHLEN),                 public :: structurefilename
        logical,                                     private ::isinit
        contains 
        procedure::load => load_Structuredata
        procedure::calc_Cartesian_coordinates
        procedure::calc_Cartesian_forces
    end type Structuredata

    interface Structuredata
        module procedure init_Structuredata
    end interface Structuredata

    contains 

    type(Structuredata) function init_Structuredata() result(xsfdata)
        xsfdata%isInit = .false.
    end function init_Structuredata

    subroutine load_Structuredata(self,file) 
        use xsflib,  only: xsf_init,                        &
                       xsf_final,                       &
                       latticeVec_in   => latticeVec,   &
                       nAtoms_in       => nAtoms,       &
                       nTypes_in       => nTypes,       &
                       atomType_in     => atomType,     &
                       hasForces_in    => hasForces,    &
                       atomTypeName_in => atomTypeName, &
                       forCart_in      => forCart,      &
                       cooCart_in      => cooCart,      &
                       pbc_in          => pbc,          &
                       hasOutput,                       &
                       E_coh, E_tot  
        use fromaenet,only:geo_get_bounds,geo_recip_lattice
        implicit none
        class(Structuredata)::self
        character(len=*), intent(in) :: file
        logical                      :: fexists

        inquire(file=file, exist=fexists)
        if (.not. fexists) then
            write(0,*) "Error: file not found in `Structuredata': ", trim(file)
            stop
        end if
        call xsf_init(file)

        self%pbc = pbc_in

        if (self%pbc) then
            self%latticeVec(:,:) = latticeVec_in(:,:)
            self%origin = (/ 0.0d0, 0.0d0, 0.0d0 /)
        else
            call geo_get_bounds(cooCart_in, self%latticeVec, self%origin)
        end if

        self%recLattVec(:,:) = geo_recip_lattice(self%latticeVec)

        self%num_of_atoms    = nAtoms_in
        self%num_of_types    = nTypes_in
        self%hasForces = hasForces_in


        if ((nAtoms_in == 0) .or. (nTypes_in == 0)) then
            write(0,*) "Error: invalid input file: ", trim(file)
            write(0,*) "       nAtoms = ", nAtoms_in
            write(0,*) "       nTypes = ", nTypes_in
            call xsf_final()
            stop
        end if

        if (allocated(self%abc_ith_atomic_coordinates)) then
            deallocate(self%abc_ith_atomic_coordinates)
        end if
        if (allocated(self%atomType)) then
            deallocate(self%atomType)
        end if
        if (allocated(self%atomTypeName)) then
            deallocate(self%atomTypeName)
        end if
        if (allocated(self%abc_ith_atomic_forces)) then
            deallocate(self%abc_ith_atomic_forces)
        end if

        allocate(self%abc_ith_atomic_coordinates(3,self%num_of_atoms),&
                self%atomType(self%num_of_atoms), self%atomTypeName(self%num_of_types))

        if (self%hasForces) allocate(self%abc_ith_atomic_forces(3,self%num_of_atoms))

        self%atomType(:)     = atomType_in(:)
        self%atomTypeName(:) = atomTypeName_in(:)

        ! convert cartesian coordinates to lattice coordinates:
        self%abc_ith_atomic_coordinates(1:3,1:self%num_of_atoms) = matmul(self%recLattVec, &
                                    cooCart_in)/(2.0d0*PI)

        ! if forces are available, convert them too:
        if (self%hasForces) then
            self%abc_ith_atomic_forces(1:3,1:self%num_of_atoms) = matmul(self%recLattVec, forCart_in)/(2.0d0*PI)
        end if

        ! store cohesive energy, if available
        if (hasOutput) then
            self%cohesiveEnergy = E_coh
            self%totalEnergy    = E_tot
            self%hasEnergy      = .true.
        else
            self%cohesiveEnergy = 0.0d0
            self%totalEnergy    = 0.0d0
            self%hasEnergy      = .false.
        end if

        self%isInit = .true.

        call xsf_final

        self%structurefilename = trim(file)


    end subroutine load_Structuredata



  !--------------------------------------------------------------------!
  !                  cartesian coordinates and forces                  !
  !--------------------------------------------------------------------!

    function calc_Cartesian_coordinates(self,iatom) result(coo)
        implicit none
        class(Structuredata)::self
        integer,            intent(in) :: iatom
        double precision, dimension(3) :: coo

        if (.not. self%isInit) then
            coo(:) = 0.0d0
            return
        end if

        coo(:) = matmul(self%latticeVec, self%abc_ith_atomic_coordinates(1:3,iatom))

    end function calc_Cartesian_coordinates

!--------------------------------------------------------------------!

    function  calc_Cartesian_forces(self,iatom) result(force)
        implicit none
        class(Structuredata)::self

        integer,            intent(in) :: iatom
        double precision, dimension(3) :: force

        if (.not. self%isInit) then
            force(:) = 0.0d0
            return
        end if

        force(:) = matmul(self%latticeVec, self%abc_ith_atomic_forces(1:3,iatom))

    end function calc_Cartesian_forces


end module structureinfo