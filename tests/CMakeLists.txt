# math_utilsライブラリの使用を明示
#add_executable(test_math_utils test_math_utils.f90)
#target_link_libraries(test_math_utils PRIVATE math_utils)

# math_utilsモジュールファイルのパスを指定
#target_include_directories(test_math_utils PRIVATE ${CMAKE_BINARY_DIR}/mod)

# テストを追加
#add_test(NAME MathUtilsTest COMMAND test_math_utils)

# math_utilsライブラリの使用を明示
add_executable(test_network test_network.f90)
target_link_libraries(test_network PRIVATE bpnet)

add_executable(test_bpnet test_bpnet.f90)
target_link_libraries(test_bpnet PRIVATE bpnet)

# math_utilsモジュールファイルのパスを指定
#target_include_directories(test_math_utils PRIVATE ${CMAKE_BINARY_DIR}/mod)

# テストを追加
add_test(NAME NetworkTest 
        COMMAND test_network 
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
#add_test(NAME BPnetTest COMMAND test_bpnet 
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests)
#configure_file(${CMAKE_SOURCE_DIR}/tests/O.fingerprint.stp ${CMAKE_BINARY_DIR}/tests/O.fingerprint.stp  COPYONLY)
#configure_file(${CMAKE_SOURCE_DIR}/tests/Ti.fingerprint.stp ${CMAKE_BINARY_DIR}/tests/Ti.fingerprint.stp  COPYONLY)