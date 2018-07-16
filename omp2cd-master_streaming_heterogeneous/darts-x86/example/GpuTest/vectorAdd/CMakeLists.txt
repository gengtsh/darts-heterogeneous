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

#SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ; -m64 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52 -gencode arch=compute_60,code=sm_60 -gencode arch=compute_60,code=compute_60)

if (DEBUG EQUAL ON)
	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ;--ptxas-options=-v -O0 -G -g -std=c++11)
else()
	SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};  --ptxas-options=-v; -std=c++11  )
endif()

include_directories(/usr/local/cude/include)
include_directories(./)

#set(CUDA_SEPARABLE_COMPILATION ON)

cuda_add_executable(vectorAdd
		    vectorAdd.cu)
#target_link_libraries(Hello darts  ${CUDA_LIBRARIES} ) 
#target_link_libraries(Hello  darts -L/usr/local/cuda/lib64 ${CUDA_LIBRARIES}  mylib  -lcuda -lcudart -lm -lrt --verbose)
target_link_libraries(vectorAdd  -L/usr/local/cuda/lib64 ${CUDA_LIBRARIES} -lcuda -lcudart  ) 
