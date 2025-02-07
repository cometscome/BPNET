program test_network
    use bpnet,only:initialize_potentials
    implicit none
    integer::num_of_Types
    character(len=2), dimension(2)         :: typeName
    character(len=1024), dimension(2) :: network_filenames
    num_of_Types = 2
    typeName(1) = "Ti"
    network_filenames(1) = "Ti.fingerprint.stp"
    typeName(2) = "O"
    network_filenames(2) = "O.fingerprint.stp"
    write(*,*) typeName,network_filenames

    call initialize_potentials(num_of_Types, typeName, network_filenames)


    stop 1
end program