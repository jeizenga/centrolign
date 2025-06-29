cmake_minimum_required(VERSION 3.10)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_SYSTEM: ${CMAKE_SYSTEM}")

find_package(Git)

set(CMAKE_CXX_STANDARD 11)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(MACOSX TRUE)
else()
	set(MACOSX FALSE)
endif()


add_compile_options(-ffast-math)

if (${DEV}) 
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -ggdb")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -ggdb")
	add_compile_options(-fsanitize=address -fsanitize=undefined)
	add_link_options(-fsanitize=address -fsanitize=undefined)
else()
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
	if (${PROFILE})
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ggdb")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
		
		if (NOT ${MACOSX})
			find_package(PkgConfig REQUIRED)
			pkg_check_modules (JEMALLOC jemalloc)
			
			pkg_search_module(JEMALLOC REQUIRED jemalloc)
			include_directories(${JEMALLOC_INCLUDE_DIRS})
		endif()
		
	endif()
endif()

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(TEST_SRC_DIR ${SRC_DIR}/test)
set(SCRIPT_SRC_DIR ${SRC_DIR}/scripts)
set(INC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

# adapted from https://stackoverflow.com/questions/1435953/how-can-i-pass-git-sha1-to-compiler-as-definition-using-cmake

# the commit's SHA1, and whether the building workspace was dirty or not
execute_process(COMMAND "${GIT_EXECUTABLE}" describe --match=NeVeRmAtCh --always --abbrev=40 --dirty
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_HASH
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

# the date of the commit
execute_process(COMMAND "${GIT_EXECUTABLE}" log -1 --format=%ad --date=local
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_DATE
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

# the subject of the commit
execute_process(COMMAND "${GIT_EXECUTABLE}" log -1 --format=%s
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE GIT_COMMIT_SUBJECT
  ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)

configure_file("${SRC_DIR}/version.cpp.pre" "${CMAKE_CURRENT_BINARY_DIR}/version.cpp" @ONLY)

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
	${SRC_DIR}/tree.cpp
	${SRC_DIR}/core.cpp
	${SRC_DIR}/trie.cpp
	${SRC_DIR}/structure_tree.cpp
	${SRC_DIR}/superbubble_distance_oracle.cpp
	${SRC_DIR}/simplifier.cpp
	${SRC_DIR}/test_util.cpp
	${SRC_DIR}/compacted_graph.cpp
	${SRC_DIR}/parameters.cpp
	${SRC_DIR}/gfa.cpp
	${SRC_DIR}/path_esa.cpp
	${SRC_DIR}/esa.cpp
	${SRC_DIR}/tree_distance_oracle.cpp
	${SRC_DIR}/bonder.cpp
	${SRC_DIR}/union_find.cpp
	${SRC_DIR}/inconsistency_identifier.cpp
	${SRC_DIR}/execution.cpp
	${SRC_DIR}/match_finder.cpp
	${SRC_DIR}/induced_match_finder.cpp
	${SRC_DIR}/packed_vector.cpp
	${SRC_DIR}/thread_pool.cpp
	${CMAKE_CURRENT_BINARY_DIR}/version.cpp
)

project(Centrolign)
add_library(centrolign SHARED ${SOURCES})
target_include_directories(centrolign PUBLIC ${INC_DIR})

add_executable(centrolign_exec ${SRC_DIR}/main.cpp)
target_link_libraries(centrolign_exec centrolign)
set_property(TARGET centrolign_exec PROPERTY OUTPUT_NAME centrolign)

if ("${PROFILE}" AND NOT "${DEV}" AND NOT "${MACOSX}")
	target_link_libraries(centrolign ${JEMALLOC_LIBRARIES})
	target_link_libraries(centrolign_exec ${JEMALLOC_LIBRARIES})
endif()

install(TARGETS centrolign LIBRARY 
	DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        PRIVATE_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")


install(TARGETS centrolign_exec RUNTIME 
	DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/)

# tests and scripts are added after the install directions
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
	test_tree
	test_gfa
	test_ruq
	test_superbubbles
	test_trie
	test_simplifier
	test_compacted_graph
	test_orthogonal_max_search_tree
	test_parameters
	test_path_esa
	test_superbubble_distance_oracle
	test_partition_algorithms
	test_match_finding
	test_bonder
	test_cactus
	test_inconsistency_identifier
	test_bit_packed_vector
)

if (${DEV}) 
	foreach(TEST ${TESTS})
		add_executable(${TEST} ${TEST_SRC_DIR}/${TEST}.cpp)
		target_include_directories(${TEST} PRIVATE ${INC_DIR})
		target_link_libraries(${TEST} centrolign)
	endforeach()
endif()

set(SCRIPT_FILES
	make_var_mat
	tree_pair_dist
	sim_centromere
	compare_truth_aln
	tree_compare
	find_universal_bridges
	remove_samples
)

if (${SCRIPTS}) 
	foreach(SCRIPT ${SCRIPT_FILES})
		add_executable(${SCRIPT} ${SCRIPT_SRC_DIR}/${SCRIPT}.cpp)
		target_include_directories(${SCRIPT} PRIVATE ${INC_DIR})
		target_link_libraries(${SCRIPT} centrolign)
	endforeach()
endif() 




