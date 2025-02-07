module bpnet
    use shared_potentials,only:initialize_potentials,&
            bpnet_atomic_energy_and_forces_novirial,&
            set_Rcmax_and_Rcmin,&
            deallocate_potentials,reload_potentials
    use sub_generate,only:generate_subroutine
    implicit none
end module