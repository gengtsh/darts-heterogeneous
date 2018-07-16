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

add_executable	(Luleshv1
				luleshMain.cpp
				LagrangeLeapFrogTP.cpp
				LagrangeNodalTP.cpp
		
				CalcForceForNodesTP.cpp
				CalcAccelerationForNodesTP.cpp
				ApplyAccelerationBoundaryConditionsForNodesTP.cpp
				CalcVelocityForNodesTP.cpp
				CalcPositionForNodesTP.cpp
				CalcVolumeForceForElemsTP.cpp
				InitStressTermsForElemsTP.cpp
				IntegrateStressForElemsTP.cpp
				IntegrateStressForElemsP2TP.cpp
				CalcHourglassControlForElemsTP.cpp
				CalcFBHourglassForceForElemsTP.cpp
		
				LagrangeElementsTP.cpp
				CalcLagrangeElementsTP.cpp
				CalcLagrangeElementsP2TP.cpp
				CalcQForElemsTP.cpp
				CalcMonotonicQRegionForElemsTP.cpp
				ApplyMaterialPropertiesForElemsTP.cpp
				UpdateVolumesForElemsTP.cpp

				lulesh.cc
				lulesh-comm.cc
				lulesh-init.cc
				lulesh-util.cc
				lulesh-viz.cc
				luleshFunction.cpp
				)
target_link_libraries(Luleshv1 darts)
set (LULESH_MPI "-DUSE_MPI=0")
#set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LULESH_MPI} " )
target_compile_options(Luleshv1 BEFORE PRIVATE ${LULESH_MPI})

#aux_source_directory(. DIR_SIRS)
#add_executable (Lulesh ${DIR_SIRS})
#target_link_libraries(Lulesh darts)
