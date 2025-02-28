add_library(
	biceps
	${CMAKE_CURRENT_SOURCE_DIR}/logger.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/interval_mesh.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/structured_mesh.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/fem_1d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/fem_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/fem_function_1d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/fem_function_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/poisson_fem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pstokes_fem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/free_surface_fem.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/time_integrator.cpp
	# ${CMAKE_CURRENT_SOURCE_DIR}/io_handler.cpp
)

# Generate different float type header based on whether USE_LONG_DOUBLE is set.
set(COMMAND_ARGS -E -P ${CMAKE_SOURCE_DIR}/include/float_type_config.hpp -o ${CMAKE_BINARY_DIR}/float_type.hpp)
if(${USE_LONG_DOUBLE})
    list(APPEND COMMAND_ARGS -DUSE_LONG_DOUBLE)
endif()
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/float_type.hpp
    COMMAND c++ ${COMMAND_ARGS}
    DEPENDS ${CMAKE_SOURCE_DIR}/include/float_type_config.hpp
)
add_custom_target(preprocess_float_type_header DEPENDS ${CMAKE_BINARY_DIR}/float_type.hpp)
add_dependencies(biceps preprocess_float_type_header)

include_directories(${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/include ${Boost_INCLUDE_DIRS})
target_link_libraries(
	biceps Eigen3::Eigen ${Boost_LIBRARIES} ${XDMF_LIBRARY} ${XDMF_CORE_LIBRARY}
)
install(TARGETS biceps DESTINATION lib)
