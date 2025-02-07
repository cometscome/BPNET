module neighborlist
    use iso_fortran_env
    implicit none

    public::max_packed_spheres


    contains


    ! 最大球数を計算する関数
    function max_packed_spheres(rmin, rmax) result(num_spheres)
        implicit none
        real(kind=8), intent(in) :: rmin  ! 小球の半径
        real(kind=8), intent(in) :: rmax  ! 容器球の半径
        integer :: num_spheres            ! 最大球数

        real(kind=8) :: volume_ratio

        ! エラーチェック
        if (rmin <= 0.0 .or. rmax <= 0.0) then
            print *, "Error: Radii must be positive values."
            num_spheres = 0
            return
        end if

        if (rmin > rmax) then
            print *, "Error: rmin must be smaller than rmax."
            num_spheres = 0
            return
        end if

        ! 球の数を計算
        volume_ratio = (rmax / rmin)**3
        num_spheres = int(0.74 * volume_ratio)
    end function max_packed_spheres

end module neighborlist