program test_math_utils
    use math_utils
    implicit none
    
    real :: result
    
    result = square(3.0)
    if (abs(result - 9.0) < 1.0e-6) then
      print *, "Test Passed"
      stop 0
    else
      print *, "Test Failed"
      stop 1
    end if
  end program test_math_utils