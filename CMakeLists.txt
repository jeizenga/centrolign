cmake_minimum_required(VERSION 3.10)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_SYSTEM: ${CMAKE_SYSTEM}")

set(CMAKE_CXX_STANDARD 11)


if (${DEV}) 
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -ggdb")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -ggdb")
	add_compile_options(-fsanitize=address -fsanitize=undefined)
	add_link_options(-fsanitize=address -fsanitize=undefined)
else()
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(TEST_SRC_DIR ${SRC_DIR}/test)
set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

include_directories(inc)
set(SOURCES
	${SRC_DIR}/gesa.cpp
	${SRC_DIR}/modify_graph.cpp
	${SRC_DIR}/graph.cpp
	${SRC_DIR}/main.cpp
	${SRC_DIR}/path_graph.cpp
	${SRC_DIR}/utility.cpp
	${SRC_DIR}/determinize.cpp
	${SRC_DIR}/anchorer.cpp
	${SRC_DIR}/chain_merge.cpp
	${SRC_DIR}/alignment.cpp
	${SRC_DIR}/stitcher.cpp
	${SRC_DIR}/logging.cpp
)

project(Centrolign)
add_library(centrolign SHARED ${SOURCES})
target_include_directories(centrolign PUBLIC ${INC_DIR})

add_executable(centrolign_exec ${SRC_DIR}/main.cpp)
target_link_libraries(centrolign_exec centrolign)
set_property(TARGET centrolign_exec PROPERTY OUTPUT_NAME centrolign)

install(TARGETS centrolign LIBRARY 
	DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        PRIVATE_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")


install(TARGETS centrolign_exec RUNTIME 
	DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/)

# tests are added after the install directions
set(TESTS
	test_rmq
	test_graph_algs
	test_utility_algs
	test_path_graph
	test_gesa
	test_chain_merge
	test_anchoring
	test_max_search_tree
	test_alignment
	test_stitcher
)

if (${DEV}) 
	foreach(TEST ${TESTS})
		add_executable(${TEST} ${TEST_SRC_DIR}/${TEST}.cpp)
		target_include_directories(${TEST} PRIVATE ${INC_DIR})
		target_link_libraries(${TEST} centrolign)
	endforeach()
endif()






