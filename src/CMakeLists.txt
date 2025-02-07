set(bpnet_PATH "${CMAKE_CURRENT_SOURCE_DIR}/")
set(bpnet_dir ${bpnet_PATH}/)

set(SOURCES 
${bpnet_PATH}/utils/bfio.f90
${bpnet_PATH}/utils/ioutils.f90
${bpnet_PATH}/utils/neighborlist.f90
${bpnet_PATH}/utils/io.f90
${bpnet_PATH}/utils/xsflib.f90
${bpnet_PATH}/utils/fromaenet.f90
${bpnet_PATH}/utils/sortlib.f90
${bpnet_PATH}/utils/lclist.f90
${bpnet_PATH}/utils/sub_trainbin2ASCII.f90
${bpnet_PATH}/utils/sub_nnASCII2bin.f90
${bpnet_PATH}/utils/cutoff.f90
${bpnet_PATH}/neuralnetworks/activation_functions.f90
${bpnet_PATH}/neuralnetworks/abstract_layer.f90
${bpnet_PATH}/neuralnetworks/dense_layer.f90
${bpnet_PATH}/neuralnetworks/multi_layers.f90
${bpnet_PATH}/neuralnetworks/aenet_potential.f90
${bpnet_PATH}/chebyshev/chebyshev.f90
${bpnet_PATH}/BPpotential.f90
${bpnet_PATH}/shared_potentials.f90
${bpnet_PATH}/sub_generate.f90
${bpnet_PATH}/sub_predict.f90
${bpnet_PATH}/inputdata.f90
${bpnet_PATH}/fingerprints.f90
${bpnet_PATH}/fingerprint_basis.f90
${bpnet_PATH}/spline/bspline.f90
${bpnet_PATH}/chebyshev/fingerprints_sub_chebyshevbasis.f90
${bpnet_PATH}/chebyshev/fingerprint_basis_chebyshev.f90
${bpnet_PATH}/spline/fingerprints_sub_splinebasis.f90
${bpnet_PATH}/spline/fingerprint_basis_spline.f90
${bpnet_PATH}/splinetensor/fingerprints_sub_splinetensorbasis.f90
${bpnet_PATH}/LJ/fingerprints_sub_LJ.f90
${bpnet_PATH}/LJ/fingerprint_basis_LJ.f90
${bpnet_PATH}/multibasis/fingerprints_sub_multibasis.f90
${bpnet_PATH}/multibasis/fingerprint_basis_multi.f90
${bpnet_PATH}/mod_dataset.f90
${bpnet_PATH}/structureinfo.f90
${bpnet_PATH}/bpnet.f90
)



add_library(bpnet STATIC)
target_sources(bpnet
PRIVATE
${SOURCES}
)

install(TARGETS bpnet EXPORT bpnetConfig DESTINATION lib)
install(FILES ${SOURCES} DESTINATION include)

export(EXPORT bpnetConfig FILE bpnetConfig.cmake)

#set(CMAKE_Fortran_FLAGS "-g -fcheck=all -Wall -fbacktrace")# -finit-real=nan -finit-integer=-999")
set(CMAKE_Fortran_FLAGS "-O3 -pedantic")# -finit-real=nan -finit-integer=-999")


add_executable(bpnet_generate.x ${bpnet_dir}/bpnet_generate.f90)
target_link_libraries(bpnet_generate.x PRIVATE bpnet)

add_executable(bpnet_predict.x ${bpnet_dir}/bpnet_predict.f90)
target_link_libraries(bpnet_predict.x PRIVATE bpnet)


add_executable(nnASCII2bin.x ${bpnet_dir}/nnASCII2bin.f90)
target_link_libraries(nnASCII2bin.x PRIVATE bpnet)

add_executable(trainbin2ASCII.x ${bpnet_dir}/trainbin2ASCII.f90)
target_link_libraries(trainbin2ASCII.x PRIVATE bpnet)

# math_utils module
add_library(math_utils math_utils.f90)

