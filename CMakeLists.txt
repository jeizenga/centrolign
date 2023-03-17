cmake_minimum_required(VERSION 3.10)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_SYSTEM: ${CMAKE_SYSTEM}")

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

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
)

project(Centrolign)
add_library(centrolign SHARED ${SOURCES})
target_include_directories(centrolign PUBLIC ${INC_DIR})

add_executable(centrolign_exec src/main.cpp)
set_property(TARGET centrolign_exec PROPERTY OUTPUT_NAME centrolign)

install(TARGETS centrolign LIBRARY 
	DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        PRIVATE_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")


install(TARGETS centrolign_exec RUNTIME 
	DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/)

# tests are added after the install directions
#add_executable(test_rmq ${TEST_SRC_DIR}/test_rmq.cpp)




