module module_activation_functions
    use iso_fortran_env
    implicit none

    type, abstract,public :: activation_function
    contains
        procedure(evaluate_a), deferred :: evaluate
        procedure(evaluate_a_prime), deferred :: evaluate_prime
        procedure(display_a), deferred :: display
    end type activation_function

    abstract interface
       real(real64) pure function evaluate_a(self,x) result(val)
            use iso_fortran_env
            import :: activation_function
            class(activation_function), intent(in) :: self
            real(real64),intent(in)::x
       end function evaluate_a

       real(real64) pure function evaluate_a_prime(self,x) result(val)
            use iso_fortran_env
            import :: activation_function
            class(activation_function), intent(in) :: self
            real(real64),intent(in)::x
        end function evaluate_a_prime

        subroutine display_a(self) 
            import :: activation_function
            class(activation_function), intent(in) :: self
        end subroutine

    end interface

    type, extends(activation_function) :: linear
        contains 
        procedure :: evaluate       => evaluate_linear
        procedure :: evaluate_prime => evaluate_linear_prime
        procedure :: display => display_linear
    end type

    type, extends(activation_function) :: relu
        contains 
        procedure :: evaluate       => evaluate_relu
        procedure :: evaluate_prime => evaluate_relu_prime
        procedure :: display => display_relu
    end type

    type, extends(activation_function) :: tanhyper
        contains 
        procedure :: evaluate       => evaluate_tanh
        procedure :: evaluate_prime => evaluate_tanh_prime
        procedure :: display => display_tanh
    end type

    contains 

    subroutine display_linear(self) 
        class(linear), intent(in) :: self
        write(*,*) "activation: linear"
    end subroutine

    subroutine display_relu(self) 
        class(relu), intent(in) :: self
        write(*,*) "activation: relu"
    end subroutine

    subroutine display_tanh(self) 
        class(tanhyper), intent(in) :: self
        write(*,*) "activation: tanh"
    end subroutine

    real(real64) pure function evaluate_linear(self,x) result(val)
        use iso_fortran_env
        class(linear), intent(in) :: self
        real(real64),intent(in)::x
        val = x
    end function evaluate_linear

    real(real64) pure function evaluate_linear_prime(self,x) result(val)
        use iso_fortran_env
        class(linear), intent(in) :: self
        real(real64),intent(in)::x
        val = 1d0
    end function evaluate_linear_prime

    real(real64) pure function evaluate_relu(self,x) result(val)
        use iso_fortran_env
        class(relu), intent(in) :: self
        real(real64),intent(in)::x
        if (x < 0d0) then
            val = 0d0
        else
            val = x
        endif
    end function evaluate_relu

    real(real64) pure function evaluate_relu_prime(self,x) result(val)
        use iso_fortran_env
        class(relu), intent(in) :: self
        real(real64),intent(in)::x
        
        if (x < 0d0) then
            val = 0d0
        else
            val = 1d0
        endif
    end function evaluate_relu_prime


    real(real64) pure function evaluate_tanh(self,x) result(val)
        use iso_fortran_env
        class(tanhyper), intent(in) :: self
        real(real64),intent(in)::x

        val = tanh(x)
    end function evaluate_tanh

    real(real64) pure function evaluate_tanh_prime(self,x) result(val)
        use iso_fortran_env
        class(tanhyper), intent(in) :: self
        real(real64),intent(in)::x
        
        val = 1d0-tanh(x)*tanh(x)
    end function evaluate_tanh_prime    


    subroutine set_activation(activation_type,activation)
        implicit none
        character(len=*),intent(in)::activation_type
        class(activation_function),allocatable,intent(out)::activation

        if (trim(adjustl(activation_type)) == "linear") then
            allocate(activation,source=linear())
            !write(*,*) "linear activation"
        else if (trim(adjustl(activation_type)) == "tanh" .or. trim(adjustl(activation_type)) == "tanhyper") then
            allocate(activation,source=tanhyper())
        else if (trim(adjustl(activation_type)) == "relu") then
            allocate(activation,source=relu())
        else
            write(*,*) "activation function ",trim(adjustl(activation_type))," is not supported"
            stop 1
        endif
    end



end module module_activation_functions