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

find_package(CUDA QUIET REQUIRED)

#set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE ON)
#set (BUILD_SHARED_LIBS OFF)
set(CUDA_PROPAGATE_HOST_FLAGS OFF)

set(CUDA_HOST_COMPILER /usr/bin/g++-4.8)

if (DEBUG EQUAL ON)	
#	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ;--ptxas-options=-v -O0 -G -g -std=c++11)
	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ;--ptxas-options=-v  -O0 -G -g -std=c++11;--default-stream legacy )
else()
#	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -Xcompiler -fPIC; --ptxas-options=-v; -std=c++11;--verbose  )
	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -Xcompiler -fPIC; --ptxas-options=-v; -std=c++11  )
#	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -Xcompiler -fPIC; --ptxas-options=-v; -std=c++11;--default-stream legacy  )
#	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -Xcompiler -fPIC; --ptxas-options=-v; -std=c++11;--default-stream per-thread  )
    #list(APPEND CUDA_NVCC_FLAGS --compiler-bindir /usr/bin/g++-4.8; -std=c++11  )
endif()

include_directories(/usr/local/cude/include)

#set(CUDA_SEPARABLE_COMPILATION ON)
cuda_add_library(StencilCudalib2.3.1.3 ${LIB_TYPE}
				 StencilCuda.cu 
		)

add_executable(StencilCuda2.3.1.3
			   StencilMain.cpp
               StencilTP.cpp
			   StencilCPUKernel.cpp
			   )
#target_link_libraries(Hello  darts -L/usr/local/cuda/lib64 ${CUDA_LIBRARIES}  mylib  -lcuda -lcudart -lm -lrt --verbose)
target_link_libraries(StencilCuda2.3.1.3  darts -L/usr/local/cuda/lib64 ${CUDA_LIBRARIES}  StencilCudalib2.3.1.3  -lcuda -lcudart -lm -lrt ) 
