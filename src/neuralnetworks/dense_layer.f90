module module_dense_layer
    use iso_fortran_env
    use module_abstract_layer,only:abstract_layer
    use module_activation_functions,only:activation_function,linear,tanhyper,set_activation
    implicit none

    type, extends(abstract_layer) :: dense_layer
        integer::indim,outdim
        real(real64),allocatable::W(:,:)
        real(real64),allocatable::b(:)
        real(real64),allocatable::y(:)
        real(real64),allocatable::x(:)
        

        !real(real64),allocatable::dLdx(:)
        real(real64),allocatable::dLdW(:,:)
        real(real64),allocatable::dLdb(:)
        logical::need_paramderiv
        class(activation_function),allocatable::activation



        contains 
        procedure :: forward       => forward_dense
        procedure :: backward => backward_dense
        procedure :: display => dense_layer_display
    end type

    interface dense_layer
        procedure::init_dense_layer
        procedure::set_dense_layer
    end interface 


    contains

    type(dense_layer) function init_dense_layer(indim,outdim,activation,need_paramderiv) result(layer)
        implicit none
        integer,intent(in)::indim,outdim
        character(len=*),intent(in)::activation
        logical,intent(in),optional::need_paramderiv

        allocate(layer%W(outdim,indim))
        allocate(layer%b(outdim))
        allocate(layer%y(outdim),source=0d0)
        allocate(layer%x(indim),source=0d0)
        allocate(layer%dLdW(outdim,indim),source=0d0)
        allocate(layer%dLdb(outdim),source=0d0)
        !allocate(layer%dLdx(outdim),source=0d0)
        !write(*,*) activation

        call set_activation(activation,layer%activation)
        

        call random_number(layer%W)
        call random_number(layer%b)
        layer%indim = indim
        layer%outdim = outdim
        layer%need_paramderiv = .false.

        if (present(need_paramderiv)) then
            layer%need_paramderiv = need_paramderiv
        endif
    end function

    type(dense_layer) function set_dense_layer(W,b,activation,need_paramderiv) result(layer)
        implicit none
        real(real64),intent(in)::W(:,:),b(:)
        character(len=*),intent(in)::activation
        logical,intent(in),optional::need_paramderiv
        integer::indim,outdim

        outdim = SIZE(W, 1)
        indim  = SIZE(W, 2)

        allocate(layer%W(outdim,indim),source=W)
        allocate(layer%b(outdim),source=b)
        allocate(layer%y(outdim),source=0d0)
        allocate(layer%x(indim),source=0d0)
        allocate(layer%dLdW(outdim,indim),source=0d0)
        allocate(layer%dLdb(outdim),source=0d0)
        !allocate(layer%dLdx(outdim),source=0d0)
        !write(*,*) activation

        call set_activation(activation,layer%activation)
        
        layer%indim = indim
        layer%outdim = outdim
        layer%need_paramderiv = .false.

        if (present(need_paramderiv)) then
            layer%need_paramderiv = need_paramderiv
        endif
    end function

    subroutine dense_layer_display(self) 
        class(dense_layer), intent(in) :: self
        write(*,*) "Dense layer."
        write(*,*) "Input and output dimensions: ",self%indim,self%outdim
        call self%activation%display()
    end subroutine dense_layer_display

    subroutine forward_dense(self,xin,xout)
        use iso_fortran_env
        class(dense_layer), intent(inout) :: self
        real(real64),intent(in)::xin(:)
        real(real64),intent(out)::xout(:)
        integer::i

        xout(1:self%outdim) = matmul(self%W,xin) + self%b(1:self%outdim)
        self%y(1:self%outdim) = xout(1:self%outdim)
        do i=1,self%outdim
            xout(i) = self%activation%evaluate(xout(i))
        end do
        

        if (self%need_paramderiv) then
            self%x(1:self%indim) = xin(1:self%indim)
        end if
    end subroutine 

    subroutine backward_dense(self,dLdyT,dLdxT)
        use iso_fortran_env
        class(dense_layer), intent(inout) :: self
        real(real64),intent(in)::dLdyT(:,:)
        real(real64),intent(out)::dLdxT(:,:)
        integer::i,j
        real(8)::sigma_prime
        real(real64),allocatable::dLdyT_temp(:,:)

        allocate(dLdyT_temp,mold=dLdyT)
        do i=1,self%outdim
            sigma_prime = self%activation%evaluate_prime(self%y(i))
            dLdyT_temp(1,i) = dLdyT(1,i)*sigma_prime 
        end do

        dLdxT = matmul(dLdyT_temp,self%W)

        !dLdx = matmul(transpose(self%W),dLdy)
        if (self%need_paramderiv) then
            do j=1,self%indim
                do i=1,self%outdim
                    self%dLdW(i,j) = dLdyT_temp(1,i)*self%x(j)
                end do
            end do
        endif
        if (self%need_paramderiv) then
            do i=1,self%outdim
                !sigma_prime = self%activation%evaluate_prime(self%y(i))
                !dLdxT(1,i) = sigma_prime*dLdxT(1,i)
                self%dLdb(i) = dLdyT_temp(1,i)
            end do
        end if

    end subroutine 


end module module_dense_layer