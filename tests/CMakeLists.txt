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

# math_utilsモジュールファイルのパスを指定
#target_include_directories(test_math_utils PRIVATE ${CMAKE_BINARY_DIR}/mod)

# テストを追加
add_test(NAME NetworkTest COMMAND test_network)