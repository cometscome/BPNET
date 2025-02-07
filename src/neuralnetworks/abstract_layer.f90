module module_abstract_layer
    implicit none

    type, abstract,public :: abstract_layer
    contains
        procedure(forward_a), deferred :: forward
        procedure(backward_a), deferred :: backward
        procedure(display_layer), deferred :: display
    end type abstract_layer

    type :: LayerWrapper
        class(abstract_layer), allocatable :: layer
        contains
        procedure :: forward
        procedure :: backward
        procedure :: display
    end type

    abstract interface
        subroutine  forward_a(self,xin,xout) 
            use iso_fortran_env
            import :: abstract_layer
            class( abstract_layer), intent(inout) :: self
            real(real64),intent(in)::xin(:)
            real(real64),intent(out)::xout(:)
        end subroutine forward_a

        subroutine  backward_a(self,dLdyT,dLdxT) 
            use iso_fortran_env
            import :: abstract_layer
            class( abstract_layer), intent(inout) :: self
            real(real64),intent(in)::dLdyT(:,:)
            real(real64),intent(out)::dLdxT(:,:)
        end subroutine backward_a
        
        subroutine display_layer(self) 
            import :: abstract_layer
            class(abstract_layer), intent(in) :: self
        end subroutine display_layer
    end interface

    contains 

    subroutine forward(self,xin,xout)
        use iso_fortran_env
        class(LayerWrapper),intent(inout)::self
        real(real64),intent(in)::xin(:)
        real(real64),intent(out)::xout(:)

        call self%layer%forward(xin,xout)
    end

    subroutine backward(self,dLdyT,dLdxT)
        use iso_fortran_env
        class(LayerWrapper),intent(inout)::self
        real(real64),intent(in)::dLdyT(:,:)
        real(real64),intent(out)::dLdxT(:,:)

        call self%layer%backward(dLdyT,dLdxT)
    end

    subroutine display(self) 
        class(LayerWrapper),intent(in)::self
        call self%layer%display()
    end subroutine display


end module module_abstract_layer