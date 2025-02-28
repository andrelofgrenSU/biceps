include_directories(
	${CMAKE_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_SOURCE_DIR}/include
	${Boost_INCLUDE_DIRS}
)

add_executable(mesh_test_suite ${CMAKE_CURRENT_SOURCE_DIR}/mesh_test_suite.cpp)
target_link_libraries(mesh_test_suite Eigen3::Eigen biceps ${Boost_LIBRARIES} ${UMFPACK_LIBRARIES})

add_executable(fem_2d_test_suite ${CMAKE_CURRENT_SOURCE_DIR}/fem_2d_test_suite.cpp)
target_link_libraries(fem_2d_test_suite Eigen3::Eigen biceps ${Boost_LIBRARIES} ${UMFPACK_LIBRARIES})

add_executable(poisson_test_suite ${CMAKE_CURRENT_SOURCE_DIR}/poisson_test_suite.cpp)
target_link_libraries(poisson_test_suite Eigen3::Eigen biceps ${Boost_LIBRARIES} ${UMFPACK_LIBRARIES})

add_executable(stokes_test_suite ${CMAKE_CURRENT_SOURCE_DIR}/stokes_test_suite.cpp)
target_link_libraries(stokes_test_suite Eigen3::Eigen biceps ${Boost_LIBRARIES} ${UMFPACK_LIBRARIES})

add_custom_target(
	test_mesh DEPENDS biceps COMMAND 
	./mesh_test_suite 
)

add_custom_target(
	test_fem_2d DEPENDS biceps COMMAND 
	./fem_2d_test_suite
)

add_custom_target(
	test_poisson DEPENDS biceps COMMAND 
	./poisson_test_suite 
)
 
add_custom_target(
	test_stokes DEPENDS biceps COMMAND 
	./stokes_test_suite
)
