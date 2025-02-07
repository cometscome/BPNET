module bspline
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer,parameter::bspline_order=3 !3rd order

    contains

    subroutine make_knotsvector(knots,d,npoints,rmax)
        real(dp), intent(out) :: knots(:)
        integer,intent(in)::d,npoints
        integer::n,i
        real(dp),intent(in)::rmax
        real(dp)::dr
        n = npoints + 2*d

        
        dr = rmax/dble(npoints-1)
        do i=1,3
            knots(i) = 0d0
            knots(n-i+1) = rmax
        end do
    
        do i=1,npoints
            knots(i+d) =  dble(i-1)*dr
        end do
    end subroutine

    subroutine bspline_basis_functions(values,x,knots,d)
        implicit none
        real(8),intent(in)::x
        integer,intent(in)::d
        real(8), intent(in) :: knots(:)
        real(8),intent(out) :: values(:)
        integer::i,i0
        values = 0d0

        i0 = find_knot_position(knots,x,size(knots),d)
        do i=i0-3,i0+1
            if (i >= 1 .and. i <= size(values)) then
                values(i) = bspline_basis(d, knots, i, x)
            end if
        end do

        !do i=1,size(values)
        !    values(i) = bspline_basis(d, knots, i, x)
        !end do
    end subroutine

    subroutine bspline_basis_functions_deriv(values,dvalues,x,knots,d)
        implicit none
        real(8),intent(in)::x
        integer,intent(in)::d
        real(8), intent(in) :: knots(:)
        real(8),intent(out) :: values(:),dvalues(:)
        integer::i,i0,n
        real(kind=8), dimension(2) :: result_value 

        values = 0d0
        dvalues = 0d0

        n = size(knots)

        !call bspline_basis_and_derivative_vectorized(knots, x, values, dvalues)
        !call bspline_basis_and_derivative(knots, x, values, dvalues)
        !return

        i0 = find_knot_position(knots,x,n,d)

        !do i=i0-3,i0+1
        do concurrent(i=i0-3:i0+1)
            if (i >= 1 .and. i <= size(values)) then
                call bspline_basis_and_derivative_non_recursive(knots, i, x, values(i), dvalues(i))

                !call bspline_basis_and_derivative_non_recursive(knots, i, x, values(i), dvalues(i))
                !call bspline_basis_with_derivative(d, knots, i, x, values(i), dvalues(i))
                !call bspline_basis_and_derivative_sub(d, knots, i, x, size(knots) - 1  ,values(i),  dvalues(i))
                !values(i) = bspline_basis(d, knots, i, x)
                !dvalues(i) = bspline_derivative(d, knots, i, x)
            end if
        end do

        !do i=1,size(values)
        !    values(i) = bspline_basis(d, knots, i, x)
        !    dvalues(i) = bspline_derivative(d, knots, i, x)
        !end do
    end subroutine

    recursive subroutine bspline_basis_with_derivative(d, knots, i, x, basis_value, derivative_value)
        integer, intent(in) :: d, i
        real(dp), intent(in) :: knots(:), x
        real(dp), intent(out) :: basis_value, derivative_value
        real(dp) :: left, right, left_deriv, right_deriv

        if (d == 0) then
            ! 0次基底関数
            if (knots(i) <= x .and. x < knots(i+1)) then
                basis_value = 1.0_dp
                derivative_value = 0.0_dp  ! 0次基底関数の微分は0
            else
                basis_value = 0.0_dp
                derivative_value = 0.0_dp
            end if
        else
            ! 再帰計算
            if (knots(i+d) - knots(i) /= 0.0_dp) then
                call bspline_basis_with_derivative(d-1, knots, i, x, left, left_deriv)
                left = (x - knots(i)) / (knots(i+d) - knots(i)) * left
                left_deriv = left_deriv / (knots(i+d) - knots(i)) + &
                            1.0_dp / (knots(i+d) - knots(i)) * left
            else
                left = 0.0_dp
                left_deriv = 0.0_dp
            end if

            if (knots(i+d+1) - knots(i+1) /= 0.0_dp) then
                call bspline_basis_with_derivative(d-1, knots, i+1, x, right, right_deriv)
                right = (knots(i+d+1) - x) / (knots(i+d+1) - knots(i+1)) * right
                right_deriv = right_deriv / (knots(i+d+1) - knots(i+1)) - &
                            1.0_dp / (knots(i+d+1) - knots(i+1)) * right
            else
                right = 0.0_dp
                right_deriv = 0.0_dp
            end if

            basis_value = left + right
            derivative_value = left_deriv + right_deriv
        end if
    end subroutine bspline_basis_with_derivative    

    recursive function bspline_basis_re(d, knots, i, x) result(basis_value)
        integer, intent(in) :: d, i
        real(dp), intent(in) :: knots(:), x
        real(dp) :: basis_value
        real(dp) :: left, right, left_denom, right_denom

        if (d == 0) then
            ! 0-th term
            basis_value = merge(1.0_dp, 0.0_dp, knots(i) <= x .and. x < knots(i+1))
        else
            ! Denominator pre-calculation to avoid repeated division
            left_denom = knots(i+d) - knots(i)
            right_denom = knots(i+d+1) - knots(i+1)

            ! Recursive left calculation
            if (left_denom /= 0.0_dp) then
                left = (x - knots(i)) / left_denom * bspline_basis_re(d-1, knots, i, x)
            else
                left = 0.0_dp
            end if

            ! Recursive right calculation
            if (right_denom /= 0.0_dp) then
                right = (knots(i+d+1) - x) / right_denom * bspline_basis_re(d-1, knots, i+1, x)
            else
                right = 0.0_dp
            end if

            ! Combine left and right
            basis_value = left + right
        end if
    end function bspline_basis_re 

    ! B-spline
    pure recursive function bspline_basis(d, knots, i, x) result(basis_value)
      integer, intent(in) :: d, i
      real(dp), intent(in) :: knots(:), x
      real(dp) :: basis_value
      real(dp) :: left, right
  
      if (d == 0) then
          ! 0-th term
          if (knots(i) <= x .and. x < knots(i+1)) then
              basis_value = 1.0_dp
          else
              basis_value = 0.0_dp
          end if
      else
          ! recursive
          if (knots(i+d) - knots(i) /= 0.0_dp) then
              left = (x - knots(i)) / (knots(i+d) - knots(i)) * bspline_basis(d-1, knots, i, x)
          else
              left = 0.0_dp
          end if
          if (knots(i+d+1) - knots(i+1) /= 0.0_dp) then
              right = (knots(i+d+1) - x) / (knots(i+d+1) - knots(i+1)) * bspline_basis(d-1, knots, i+1, x)
          else
              right = 0.0_dp
          end if
          basis_value = left + right
      end if
    end function bspline_basis

    subroutine bspline_basis_and_derivative_sub(d, knots, i, x, n,basis_value, derivative_value)
        implicit none
        integer, intent(in) :: d, i                  ! d: 次数, i: 基底関数のインデックス
        real(dp), intent(in) :: knots(:), x         ! knots: ノットベクトル, x: 評価点
        real(dp), intent(out) :: basis_value        ! 基底関数値
        real(dp), intent(out) :: derivative_value   ! 微分値
        integer,intent(in)::n
        real(dp):: Nf(0:d, 0:n), dNf(0:d, 0:n)    ! N: 基底関数, dN: 微分基底関数

        integer :: j, k
    
        !n = size(knots) - 1                         ! ノットベクトルのサイズ
        !allocate(N(0:d, 0:n), dN(0:d, 0:n))         ! 0次からd次までの値を格納
    
        ! 初期化: 0次B-spline基底関数
        do j = 0, n-1
            if (knots(j) <= x .and. x < knots(j+1)) then
                Nf(0, j) = 1.0_dp
            else
                Nf(0, j) = 0.0_dp
            end if
            dNf(0, j) = 0.0_dp  ! 微分の初期値
        end do
    
        ! 再帰関係をループで計算
        do k = 1, d
            do j = 0, n-k-1
                ! B-spline基底関数の値
                if (knots(j+k) - knots(j) /= 0.0_dp) then
                    Nf(k, j) = (x - knots(j)) / (knots(j+k) - knots(j)) * Nf(k-1, j)
                else
                    Nf(k, j) = 0.0_dp
                end if
                if (knots(j+k+1) - knots(j+1) /= 0.0_dp) then
                    Nf(k, j) = Nf(k, j) + (knots(j+k+1) - x) / (knots(j+k+1) - knots(j+1)) * Nf(k-1, j+1)
                end if
    
                ! 微分値の計算
                if (knots(j+k) - knots(j) /= 0.0_dp) then
                    dNf(k, j) = Nf(k-1, j) / (knots(j+k) - knots(j))
                else
                    dNf(k, j) = 0.0_dp
                end if
                if (knots(j+k+1) - knots(j+1) /= 0.0_dp) then
                    dNf(k, j) = dNf(k, j) - Nf(k-1, j+1) / (knots(j+k+1) - knots(j+1))
                end if
            end do
        end do
    
        ! 結果を返す
        basis_value = Nf(d, i)    ! 基底関数値
        derivative_value = dNf(d, i)   ! 微分値
    
        !deallocate(N, dN)
    end subroutine bspline_basis_and_derivative_sub

    integer function find_knot_position(x,x00,nn,dd) result(low)
        integer,intent(in)::nn,dd
        real(8),intent(in)::x(nn)
        real(8),intent(in)::x00
        integer :: high, mid

        low = 1+dd
        high = nn-dd
        if(x(high) -x00 .eq. 0d0) then
            low = high
            return
        end if


        do while (high - low > 1)
            mid = (low + high) / 2
            if (x(mid) > x00) then
                high = mid
            else
                low = mid
            end if
        end do
    end function


  
    ! B-spline
    pure recursive function bspline_derivative(d, knots, i, x) result(deriv_value)
      integer, intent(in) :: d, i
      real(dp), intent(in) :: knots(:), x
      real(dp) :: deriv_value
      real(dp) :: left, right
  
      if (d == 0) then
          ! zero for zero-th term
          deriv_value = 0.0_dp
      else
          ! left hand side
          if (knots(i+d) - knots(i) /= 0.0_dp) then
              left = d / (knots(i+d) - knots(i)) * bspline_basis(d-1, knots, i, x)
          else
              left = 0.0_dp
          end if
  
          ! right hand side
          if (knots(i+d+1) - knots(i+1) /= 0.0_dp) then
              right = d / (knots(i+d+1) - knots(i+1)) * bspline_basis(d-1, knots, i+1, x)
          else
              right = 0.0_dp
          end if
  
          deriv_value = left - right
      end if
    end function bspline_derivative


    subroutine bspline_basis_and_derivative(knots, x, basis_values, derivatives)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        implicit none
        real(dp), intent(in) :: knots(:)                 ! ノットベクトル
        real(dp), intent(in) :: x                       ! 評価点
        real(dp), intent(out) :: basis_values(size(knots)-4) ! 基底関数の値 (d=3)
        real(dp), intent(out) :: derivatives(size(knots)-4)  ! 基底関数の微分 (d=3)
        real(dp) :: n0(size(knots)-1), n1(size(knots)-2), n2(size(knots)-3), n3(size(knots)-4)
        real(dp) :: dn1(size(knots)-2), dn2(size(knots)-3), dn3(size(knots)-4)
        integer :: i
    
        ! 初期化: 0次基底関数 (d=0)
        !n0 = 0.0_dp
        do concurrent (i = 1:size(knots)-1)
        !do i = 1, size(knots)-1
            if (knots(i) <= x .and. x < knots(i+1)) then
                n0(i) = 1.0_dp
            else
                n0(i) = 0.0_dp
            end if
        end do
    
        ! 1次基底関数とその微分 (d=1)
        !dn1 = 0.0_dp
        do concurrent (i = 1:size(knots)-2)
        !do i = 1, size(knots)-2
            if (knots(i+1) - knots(i) /= 0.0_dp) then
                n1(i) = (x - knots(i)) / (knots(i+1) - knots(i)) * n0(i)
                dn1(i) = 1.0_dp / (knots(i+1) - knots(i)) * n0(i)
            else
                n1(i) = 0.0_dp
                dn1(i) = 0.0_dp
            end if
            if (knots(i+2) - knots(i+1) /= 0.0_dp) then
                n1(i) = n1(i) + (knots(i+2) - x) / (knots(i+2) - knots(i+1)) * n0(i+1)
                dn1(i) = dn1(i) - 1.0_dp / (knots(i+2) - knots(i+1)) * n0(i+1)
            end if
        end do
    
        ! 2次基底関数とその微分 (d=2)
        !dn2 = 0.0_dp
        do concurrent (i = 1:size(knots)-3)
        !do i = 1, size(knots)-3
            if (knots(i+2) - knots(i) /= 0.0_dp) then
                n2(i) = (x - knots(i)) / (knots(i+2) - knots(i)) * n1(i)
                dn2(i) = (1.0_dp / (knots(i+2) - knots(i))) * n1(i) + &
                     ((x - knots(i)) / (knots(i+2) - knots(i))) * dn1(i)
            else
                n2(i) = 0.0_dp
                dn2(i) = 0.0_dp
            end if
            if (knots(i+3) - knots(i+1) /= 0.0_dp) then
                n2(i) = n2(i) + (knots(i+3) - x) / (knots(i+3) - knots(i+1)) * n1(i+1)
                dn2(i) = dn2(i) - (1.0_dp / (knots(i+3) - knots(i+1))) * &
                    n1(i+1) + ((knots(i+3) - x) / (knots(i+3) - knots(i+1))) * dn1(i+1)
            end if
        end do
    
        ! 3次基底関数とその微分 (d=3)
        !dn3 = 0.0_dp
        do concurrent (i = 1:size(knots)-4)
        !do i = 1, size(knots)-4
            if (knots(i+3) - knots(i) /= 0.0_dp) then
                n3(i) = (x - knots(i)) / (knots(i+3) - knots(i)) * n2(i)
                dn3(i) = (1.0_dp / (knots(i+3) - knots(i))) * n2(i) +&
                     ((x - knots(i)) / (knots(i+3) - knots(i))) * dn2(i)
            else
                n3(i) = 0.0_dp
                dn3(i) = 0.0_dp
            end if
            if (knots(i+4) - knots(i+1) /= 0.0_dp) then
                n3(i) = n3(i) + (knots(i+4) - x) / (knots(i+4) - knots(i+1)) * n2(i+1)
                dn3(i) = dn3(i) - (1.0_dp / (knots(i+4) - knots(i+1))) *&
                 n2(i+1) + ((knots(i+4) - x) / (knots(i+4) - knots(i+1))) * dn2(i+1)
            end if
        end do
    
        ! 結果を返す
        basis_values = n3
        derivatives = dn3
    end subroutine bspline_basis_and_derivative



    pure function bspline_basis_non_recursive(knots, i, x) result(basis_value)
        integer, intent(in) :: i
        real(dp), intent(in) :: knots(:), x
        real(dp) :: basis_value
        real(dp) :: N0(3), N1(2), N2(2)
        integer :: j

        ! Step 0: Initialize the 0-th degree basis functions
        N0(:) = 0.0_dp
        do j = 0, 3
            if (knots(i+j) <= x .and. x < knots(i+j+1)) then
                N0(j+1) = 1.0_dp
            end if
        end do

        ! Step 1: Compute 1st degree basis functions
        N1(:) = 0.0_dp
        do j = 0, 2
            if (knots(i+j+1) - knots(i+j) /= 0.0_dp) then
                N1(j+1) = (x - knots(i+j)) / (knots(i+j+1) - knots(i+j)) * N0(j+1)
            end if
            if (knots(i+j+2) - knots(i+j+1) /= 0.0_dp) then
                N1(j+1) = N1(j+1) + (knots(i+j+2) - x) / (knots(i+j+2) - knots(i+j+1)) * N0(j+2)
            end if
        end do

        ! Step 2: Compute 2nd degree basis functions
        N2(:) = 0.0_dp
        do j = 0, 1
            if (knots(i+j+2) - knots(i+j) /= 0.0_dp) then
                N2(j+1) = (x - knots(i+j)) / (knots(i+j+2) - knots(i+j)) * N1(j+1)
            end if
            if (knots(i+j+3) - knots(i+j+1) /= 0.0_dp) then
                N2(j+1) = N2(j+1) + (knots(i+j+3) - x) / (knots(i+j+3) - knots(i+j+1)) * N1(j+2)
            end if
        end do

        ! Step 3: Compute 3rd degree basis functions
        basis_value = 0.0_dp
        if (knots(i+3) - knots(i) /= 0.0_dp) then
            basis_value = (x - knots(i)) / (knots(i+3) - knots(i)) * N2(1)
        end if
        if (knots(i+4) - knots(i+1) /= 0.0_dp) then
            basis_value = basis_value + (knots(i+4) - x) / (knots(i+4) - knots(i+1)) * N2(2)
        end if
    end function bspline_basis_non_recursive   
    
    
    pure subroutine bspline_basis_and_derivative_non_recursive(knots, i, x, basis_value, derivative_value)
        implicit none
        integer, intent(in) :: i
        real(dp), intent(in) :: knots(:), x
        real(dp), intent(out) :: basis_value, derivative_value
        real(dp) :: N0(3), N1(3), N2(3)
        real(dp) :: dN0(3), dN1(3), dN2(3)
        integer :: j
    
        ! Initialize arrays
        N0(:) = 0.0_dp
        dN0(:) = 0.0_dp
    
        ! Step 0: Compute 0-th degree basis functions
        do j = 0, 3
            if (knots(i+j) <= x .and. x < knots(i+j+1)) then
                N0(j+1) = 1.0_dp
            end if
        end do
    
        ! Step 1: Compute 1st degree basis functions and their derivatives
        N1(:) = 0.0_dp
        dN1(:) = 0.0_dp
        do j = 0, 2
            if (knots(i+j+1) - knots(i+j) /= 0.0_dp) then
                N1(j+1) = (x - knots(i+j)) / (knots(i+j+1) - knots(i+j)) * N0(j+1)
                dN1(j+1) = 1.0_dp / (knots(i+j+1) - knots(i+j)) * N0(j+1)
            end if
            if (knots(i+j+2) - knots(i+j+1) /= 0.0_dp) then
                N1(j+1) = N1(j+1) + (knots(i+j+2) - x) / (knots(i+j+2) - knots(i+j+1)) * N0(j+2)
                dN1(j+1) = dN1(j+1) - 1.0_dp / (knots(i+j+2) - knots(i+j+1)) * N0(j+2)
            end if
        end do
    
        ! Step 2: Compute 2nd degree basis functions and their derivatives
        N2(:) = 0.0_dp
        dN2(:) = 0.0_dp
        do j = 0, 1
            if (knots(i+j+2) - knots(i+j) /= 0.0_dp) then
                N2(j+1) = (x - knots(i+j)) / (knots(i+j+2) - knots(i+j)) * N1(j+1)
                dN2(j+1) = 1.0_dp / (knots(i+j+2) - knots(i+j)) * N1(j+1) + &
                     (x - knots(i+j)) / (knots(i+j+2) - knots(i+j)) * dN1(j+1)
            end if
            if (knots(i+j+3) - knots(i+j+1) /= 0.0_dp) then
                N2(j+1) = N2(j+1) + (knots(i+j+3) - x) / (knots(i+j+3) - knots(i+j+1)) * N1(j+2)
                dN2(j+1) = dN2(j+1) - 1.0_dp / (knots(i+j+3) - knots(i+j+1)) * &
                     N1(j+2) + (knots(i+j+3) - x) / (knots(i+j+3) - knots(i+j+1)) * dN1(j+2)
            end if
        end do
    
        ! Step 3: Compute 3rd degree basis functions and their derivatives
        basis_value = 0.0_dp
        derivative_value = 0.0_dp
        if (knots(i+3) - knots(i) /= 0.0_dp) then
            basis_value = (x - knots(i)) / (knots(i+3) - knots(i)) * N2(1)
            derivative_value = 1.0_dp / (knots(i+3) - knots(i)) * N2(1) + &
                 (x - knots(i)) / (knots(i+3) - knots(i)) * dN2(1)
        end if
        if (knots(i+4) - knots(i+1) /= 0.0_dp) then
            basis_value = basis_value + (knots(i+4) - x) / (knots(i+4) - knots(i+1)) * N2(2)
            derivative_value = derivative_value - 1.0_dp / (knots(i+4) -&
             knots(i+1)) * N2(2) + (knots(i+4) - x) / (knots(i+4) - knots(i+1)) * dN2(2)
        end if
    end subroutine bspline_basis_and_derivative_non_recursive
 
end module