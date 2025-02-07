module module_multi_layers
    use iso_fortran_env
    use module_aenet_potential,only:aenet_potential,choose_activationfunction_fromaenet
    use module_abstract_layer,only:abstract_layer,LayerWrapper
    use module_dense_layer,only:dense_layer
    implicit none

    type, public :: Multi_layers
        type(LayerWrapper),allocatable::layers(:)
        integer::numlayers
        integer::nnodes_max
        real(real64),allocatable::x_temps(:,:,:)
        integer,          dimension(:),   allocatable :: nnodes


        contains 
        procedure :: display => multi_layers_display
        procedure :: forward => multi_layers_forward
        procedure :: backward => multi_layers_backward
    end type

    interface Multi_layers
        module procedure init_Multi_layers_fromaenet
    end interface Multi_layers

    contains 


    type(Multi_layers) function init_Multi_layers_fromaenet(aenetparams) result(potential)
        implicit none
        type(aenet_potential),intent(in)::aenetparams
        integer::ilayer
        integer::indim,outdim
        character(len=128)::activation
        integer::iw1,iw2
        integer               :: nnodes1, nnodes2
        integer, dimension(2) :: Wshape

        potential%numlayers = aenetparams%nlayers
        allocate(potential%layers(potential%numlayers-1))

        iw1 = 1
        nnodes2 = 1 
        nnodes1 = aenetparams%nnodes(1)     
        
        potential%nnodes_max = aenetparams%nnodes_max
        allocate(potential%x_temps(1,1:potential%nnodes_max,2),source=0d0)
        allocate(potential%nnodes,source=aenetparams%nnodes)
        

        do ilayer=1,potential%numlayers-1
            indim = aenetparams%nnodes(ilayer)
            outdim = aenetparams%nnodes(ilayer+1)
            activation = choose_activationfunction_fromaenet(aenetparams%f_a(ilayer+1))

            iw2 = aenetparams%iw(ilayer+1)
            
            nnodes2 = aenetparams%nnodes(ilayer+1)

            Wshape(1:2) = (/ nnodes2, nnodes1 /)

            
            allocate(potential%layers(ilayer)%layer,source=dense_layer(&
                reshape(aenetparams%W(iw1:iw2),Wshape),&
                aenetparams%W(iw2-nnodes2+1:iw2),&
                activation))

            iw1 = iw2 + 1
            nnodes1 = nnodes2
        end do


    end

    subroutine multi_layers_display(self)
        implicit none
        class(Multi_layers),intent(in)::self
        integer::ilayer
        write(*,*) "nodes: ",self%nnodes

        do ilayer=1,self%numlayers-1
            write(*,*) ilayer,"-th layer"
            call self%layers(ilayer)%display()
        end do
    end subroutine

    subroutine multi_layers_forward(self,xin,xout)
        use iso_fortran_env
        class(Multi_layers), intent(inout) :: self
        real(real64),intent(in)::xin(:)
        real(real64),intent(out)::xout(:)
        integer::ilayer
        integer::index1,index2,indextemp

        index1 = 1
        index2 = 2
        
        self%x_temps(1,1:self%nnodes(1),index1) = xin(1:self%nnodes(1))
        do ilayer=1,self%numlayers-1
            call self%layers(ilayer)%forward(self%x_temps(1,1:self%nnodes(ilayer),index1),&
                    self%x_temps(1,1:self%nnodes(ilayer+1),index2))
            indextemp = index2
            index2 = index1
            index1 = indextemp
        end do
        xout(1:self%nnodes(self%numlayers)) = self%x_temps(1,1:self%nnodes(self%numlayers),index1)


    end subroutine 

    subroutine multi_layers_backward(self,dLdy,dLdx)
        use iso_fortran_env
        class(Multi_layers), intent(inout) :: self
        real(real64),intent(in)::dLdy(:)
        real(real64),intent(out)::dLdx(:)
        integer::ilayer
        integer::index1,index2,indextemp

        index1 = 1
        index2 = 2
        
        self%x_temps(1,1:self%nnodes(self%numlayers),index1) = dLdy(1:self%nnodes(self%numlayers))
        do ilayer=self%numlayers-1,1,-1
            !write(*,*) ilayer,self%x_temps(1:1,1:self%nnodes(ilayer+1),index1)
            call self%layers(ilayer)%backward(self%x_temps(1:1,1:self%nnodes(ilayer+1),index1),&
                    self%x_temps(1:1,1:self%nnodes(ilayer),index2))
            indextemp = index2
            index2 = index1
            index1 = indextemp
        end do
        dLdx(1:self%nnodes(1)) = self%x_temps(1,1:self%nnodes(1),index1)
    end subroutine 



end module module_multi_layers