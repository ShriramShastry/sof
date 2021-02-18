# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required(VERSION 3.10)

project(SOF_MATH_BENCH C)

include(../../scripts/cmake/misc.cmake)

add_executable(math_bench
	main.c
)

sof_append_relative_path_definitions(math_bench)

set(sof_source_directory "${PROJECT_SOURCE_DIR}/../..")
set(sof_install_directory "${PROJECT_BINARY_DIR}/sof_ep/install")
set(sof_binary_directory "${PROJECT_BINARY_DIR}/sof_ep/build")

set(config_h ${sof_binary_directory}/library_autoconfig.h)

target_compile_options(math_bench PRIVATE -g -O3 -Wall -Werror -Wmissing-prototypes
  -Wimplicit-fallthrough -DCONFIG_LIBRARY -imacros${config_h})

target_link_libraries(math_bench PRIVATE -ldl -lm)

install(TARGETS math_bench DESTINATION bin)


include(ExternalProject)

ExternalProject_Add(sof_ep
	DOWNLOAD_COMMAND ""
	SOURCE_DIR "${sof_source_directory}"
	PREFIX "${PROJECT_BINARY_DIR}/sof_ep"
	BINARY_DIR "${sof_binary_directory}"
	CMAKE_ARGS -DCONFIG_LIBRARY=ON
		-DCONFIG_LIBRARY_STATIC=ON
		-DCMAKE_INSTALL_PREFIX=${sof_install_directory}
		-DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
		-DCONFIG_H_PATH=${config_h}
	BUILD_ALWAYS 1
	BUILD_BYPRODUCTS "${sof_install_directory}/lib/libsof.so"
)

ExternalProject_Add_Step(
	sof_ep defconfig
	COMMAND ${CMAKE_COMMAND} --build . --target library_defconfig
	DEPENDEES configure
	DEPENDERS build
	WORKING_DIRECTORY "${sof_binary_directory}"
)

add_library(sof_library STATIC IMPORTED)
set_target_properties(sof_library PROPERTIES IMPORTED_LOCATION "${sof_install_directory}/lib/libsof.a")
add_dependencies(sof_library sof_ep)

target_link_libraries(math_bench PRIVATE sof_library)
target_include_directories(math_bench PRIVATE ${sof_install_directory}/include)

set_target_properties(math_bench
	PROPERTIES
	INSTALL_RPATH "${sof_install_directory}/lib"
	INSTALL_RPATH_USE_LINK_PATH TRUE
)
