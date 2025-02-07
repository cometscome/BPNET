module math_utils
    implicit none
  contains
    function square(x) result(res)
      real, intent(in) :: x
      real :: res
      res = x * x
    end function square
  end module math_utils