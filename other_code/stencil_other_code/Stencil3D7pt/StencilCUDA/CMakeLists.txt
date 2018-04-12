################################################################################
#                                                                              #
# Copyright (c) 2011-2014, University of Delaware                              # 
# All rights reserved.                                                         #
#                                                                              #
# Redistribution and use in source and binary forms, with or without           # 
# modification, are permitted provided that the following conditions           # 
# are met:                                                                     #
#                                                                              #
# 1. Redistributions of source code must retain the above copyright            # 
# notice, this list of conditions and the following disclaimer.                # 
#                                                                              #
# 2. Redistributions in binary form must reproduce the above copyright         # 
# notice, this list of conditions and the following disclaimer in the          # 
# documentation and/or other materials provided with the distribution.         # 
#                                                                              #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS          # 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT            # 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS            # 
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE               # 
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,         # 
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,         # 
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;             # 
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER             # 
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT           # 
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN            # 
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE              # 
# POSSIBILITY OF SUCH DAMAGE.                                                  # 
#                                                                              #
################################################################################

cmake_minimum_required( VERSION 2.6 )

ENABLE_TESTING()

#set(CMAKE_C_COMPILER "/usr/bin/gcc-4.8")
#set(CMAKE_CXX_COMPILER "/usr/bin/g++-4.8")

project ( StencilCUDA )



#include_directories(/usr/include/c++/4.4.4/include)

#include_directories(/usr/include/c++/4.8.2/
#					include)

include_directories(include)

include_directories(/usr/local/cude/include)

if( "${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local" )
message( "HOME" ) 
set (CMAKE_INSTALL_PREFIX $ENV{HOME} )
endif( "${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local" )

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")



find_package(CUDA QUIET REQUIRED)

message(¨Special Release¨)
set(CMAKE_BUILD_TYPE  ¨RELEASE¨)
option(DEBUG "DEBUG" OFF)
if (DEBUG EQUAL ON)
	message("DEBUG VERSION!!!")
	set (CMAKE_BUILD_TYPE "DEBUG")
else()
	message("Release Version")
	set (CMAKE_BUILD_TYPE "RELEASE")
endif()


##### Reset the CXX FLAGS so it is empty
###if (DEBUG EQUAL ON)
###	#SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};--ptxas-options=-v -O0 -G -g -ccbin /usr/bin/cc)
###
###	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};--ptxas-options=-v -O0 -G -g -ccbin /usr/bin/g++-4.8)
###	SET( CMAKE_C_FLAGS  "-Wall -Wextra  -O0 -g -std=c11 -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=200809L " )
###	SET( CMAKE_CXX_FLAGS  "-Wall -Wextra  -O0 -g" )
###	message("DEBUG MODE")
###else()
###	#SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};--ptxas-options=-v -O3 -ccbin /usr/bin/cc)
###	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};--ptxas-options=-v -O3 -ccbin /usr/bin/g++-4.8)
###	SET( CMAKE_C_FLAGS  "-Wall -Wextra -O3 -std=c11 -D_XOPEN_SOURCE -D_POSIX_C_SOURCE=200809L " )
###	SET( CMAKE_CXX_FLAGS  "-Wall -Wextra -O3" )
###endif()

set(CUDA_PROPAGATE_HOST_FLAGS OFF)
set(CUDA_HOST_COMPILER /usr/bin/g++-4.8)

if (DEBUG EQUAL ON)	
#	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ;--ptxas-options=-v -O0 -G -g -std=c++11)
	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ;--ptxas-options=-v  -O0 -G -g -std=c++11;--default-stream legacy )
else()
#	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -Xcompiler -fPIC; --ptxas-options=-v; -std=c++11;--verbose  )
	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -Xcompiler -fPIC; --ptxas-options=-v;  -std=c++11)
#	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -Xcompiler -fPIC; --ptxas-options=-v; -std=c++11;--default-stream legacy  )
#	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -Xcompiler -fPIC; --ptxas-options=-v; -std=c++11;--default-stream per-thread  )
    #list(APPEND CUDA_NVCC_FLAGS --compiler-bindir /usr/bin/g++-4.8; -std=c++11  )


endif()


#SET(LIB_TYPE STATIC) #set the lib type
CUDA_ADD_LIBRARY(CUDA_Stencil_lib ${LIB_TYPE} 
		StencilCuda.cu)


#cuda_add_executable(StencilCUDA
#		main.c
#		timer.c
#		)

add_executable(StencilCUDA
		main.c
		timer.c
		)

target_link_libraries(StencilCUDA  -L/usr/local/cuda/lib64 ${CUDA_LIBRARIES} CUDA_Stencil_lib -lpthread  -lm -lrt -lcuda -lcudart  -lstdc++)
