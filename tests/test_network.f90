program test_network
    use module_activation_functions,only:linear,relu,tanhyper
    use module_dense_layer,only:dense_layer
    use module_abstract_layer,only:abstract_layer
    implicit none
    class(abstract_layer),allocatable::layer
    !type(dense_layer)::dense
    integer::indim,outdim
    real(8),allocatable::xin(:)
    real(8),allocatable::xout(:)

    indim = 44
    outdim = 10
    allocate(xin(indim))
    allocate(xout(outdim))
    !allocate(layer,source=dense_layer(indim,outdim,"linear"))
    allocate(layer,source=dense_layer(indim,outdim,"tanh"))

    call random_number(xin)
    
    call layer%forward(xin,xout)
    write(*,*) xout

    

    stop 1
end program