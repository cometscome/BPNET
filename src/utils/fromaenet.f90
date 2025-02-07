module fromaenet
    use iso_fortran_env
    implicit none
    real(real64), parameter, public :: PI = 3.141592653589793d0


    public::geo_get_bounds,geo_recip_lattice,geo_type_conv
    contains


!--------------------------------------------------------------------!
!        deduce bounds from coordinates of isolated structure        !
! Returns input coordinates shifted by 'shift'.                      !
!--------------------------------------------------------------------!

    subroutine geo_get_bounds(cooCart, scal, shift)

        implicit none

        double precision, dimension(:,:), intent(inout) :: cooCart
        double precision, dimension(3,3), intent(out)   :: scal
        double precision, dimension(3),   intent(out)   :: shift

        double precision :: x_min, x_max
        double precision :: y_min, y_max
        double precision :: z_min, z_max

        integer :: iat

        ! To avoid numerical problems we make the bounding box
        ! slightly larger than necessary and also consider this in
        ! the shift of the coordinates. This should guarantee that all
        ! scaled coordinates are in [0,1).
        double precision, parameter :: EPS = 2.0d-6

        x_min = minval(cooCart(1,:))
        x_max = maxval(cooCart(1,:))
        y_min = minval(cooCart(2,:))
        y_max = maxval(cooCart(2,:))
        z_min = minval(cooCart(3,:))
        z_max = maxval(cooCart(3,:))

        ! orthogonal bounding box
        scal(:,1) = (/ x_max - x_min + EPS, 0.0d0, 0.0d0  /)
        scal(:,2) = (/ 0.0d0, y_max - y_min + EPS, 0.0d0  /)
        scal(:,3) = (/ 0.0d0, 0.0d0, z_max - z_min + EPS  /)

        ! origin of the bounding box
        shift(1) = x_min - 0.5d0*EPS
        shift(2) = y_min - 0.5d0*EPS
        shift(3) = z_min - 0.5d0*EPS

        ! shift coordinates to bounding box
        do iat = 1, size(cooCart(1,:))
        cooCart(1:3,iat) = cooCart(1:3,iat) - shift(1:3)
        end do

    end subroutine geo_get_bounds

  !--------------------------------------------------------------------!
  !           calculation of the reciprocal lattice vectors            !
  !                                                                    !
  ! if (cryst == .true.) the crystallographic reciprocal lattice will  !
  ! be returned, i.e., the vectors are not scaled by 2*PI              !
  !--------------------------------------------------------------------!

    function geo_recip_lattice(avec, cryst) result(bvec)

        implicit none
    
        double precision, dimension(3,3), intent(in) :: avec
        double precision, dimension(3,3)             :: bvec
        logical, optional,                intent(in) :: cryst
    
        double precision :: V
    
        bvec(1,1:3) =  vproduct(avec(1:3,2), avec(1:3,3))
        bvec(2,1:3) = -vproduct(avec(1:3,1), avec(1:3,3))
        bvec(3,1:3) =  vproduct(avec(1:3,1), avec(1:3,2))
    
        V = geo_cell_volume(avec)
        bvec(:,:) = bvec(:,:)/V
    
        if (present(cryst)) then
           if (.not. cryst) bvec(:,:) = bvec(:,:)*2.0d0*PI
        else
           bvec(:,:) = bvec(:,:)*2.0d0*PI
        end if
    
    end function geo_recip_lattice

  !--------------------------------------------------------------------!
  !                            cell volume                             !
  !--------------------------------------------------------------------!

    function geo_cell_volume(avec) result(V)

        implicit none
    
        double precision, dimension(3,3), intent(in) :: avec
        double precision                             :: V
    
        V = avec(1,1)*avec(2,2)*avec(3,3) &
          + avec(2,1)*avec(3,2)*avec(1,3) &
          + avec(3,1)*avec(1,2)*avec(2,3) &
          - avec(3,1)*avec(2,2)*avec(1,3) &
          - avec(1,1)*avec(3,2)*avec(2,3) &
          - avec(2,1)*avec(1,2)*avec(3,3)
    
    end function geo_cell_volume    

  !------------------------------------------------------------------!
  !                       vector/cross product                       !
  !------------------------------------------------------------------!

    function vproduct(a,b) result(c)

        implicit none
    
        double precision, dimension(3), intent(in) :: a, b
        double precision, dimension(3)             :: c
    
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    
    end function vproduct

!--------------------------------------------------------------------!
!                        convert type number                         !
!                                                                    !
! For the case of two different atom type indices:                   !
! (1) the atom types of the parametrization                          !
! (2) the atom types of the currently loaded structure               !
! This routine allows to convert between the two.                    !
! 0 will be returned, if the species is not found.                   !
!--------------------------------------------------------------------!

    function geo_type_conv(it1, nT1, name1, nT2, name2) result(it2)

        implicit none
    
        integer,                          intent(in) :: it1
        integer,                          intent(in) :: nT1
        character(len=*), dimension(nT1), intent(in) :: name1
        integer,                          intent(in) :: nT2
        character(len=*), dimension(nT2), intent(in) :: name2
        integer                                      :: it2
    
        integer :: itype
    
        it2 = 0
        search : do itype = 1, nT2
           if (trim(name2(itype)) == trim(name1(it1))) then
              it2 = itype
              exit search
           end if
        end do search
    
    end function geo_type_conv    

end module fromaenet