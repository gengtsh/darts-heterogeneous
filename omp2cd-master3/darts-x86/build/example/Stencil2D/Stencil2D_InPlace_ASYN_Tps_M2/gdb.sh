
b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2D_main_tmp.cpp:main


b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DPartition.h:56
b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DPartition.cpp:Stencil2DPartitionChunks::fire
b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DPartition.cpp:Stencil2DPartitionCheckTP::fire
b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DPartition.cpp:Stencil2DPartitionSync::fire
#

b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DRowDecomposition.h:95
#b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DRowDecomposition.cpp:Stencil2DRowCheck::fire
#b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DRowDecomposition.cpp:Stencil2DRowLoopCopyUp::fire
#b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DRowDecomposition.cpp:Stencil2DRowLoopCopyDown::fire
#b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DRowDecomposition.cpp:Stencil2DRowLoop::fire
b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DRowDecomposition.cpp:Stencil2DRowSync::fire
#
#b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DKernel.cpp:stencil2d_seq
#b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DKernel.cpp:computeInner_stencil2d
#b /home/tgeng/omp2cd-master/darts-x86/example/Stencil2D/Stencil2D_InPlace_ASYN_Tps_M2/Stencil2DKernel.cpp:copyLine_stencil2d
##
#b /home/tgeng/DARTS/src/runtime/Runtime.cpp:157
#b /home/tgeng/DARTS/src/runtime/Runtime.cpp:280
#b /home/tgeng/DARTS/src/runtime/Runtime.cpp:283
#b /home/tgeng/DARTS/src/runtime/Runtime.cpp:105
#
#b /home/tgeng/DARTS/src/scheduler/MSchedPolicy.cpp:132
#b /home/tgeng/DARTS/src/scheduler/MSchedPolicy.cpp:153
#b /home/tgeng/DARTS/src/scheduler/MSchedPolicy.cpp:45
#b /home/tgeng/DARTS/src/scheduler/MSchedPolicy.cpp:66
#
#
#b /home/tgeng/DARTS/src/scheduler/TPSchedPolicy.cpp:203
#b /home/tgeng/DARTS/src/scheduler/TPSchedPolicy.cpp:214
#b /home/tgeng/DARTS/src/scheduler/TPSchedPolicy.cpp:233
#b /home/tgeng/DARTS/src/scheduler/TPSchedPolicy.cpp:78
#b /home/tgeng/DARTS/src/scheduler/TPSchedPolicy.cpp:119





set print thread-events on
show print thread-events
#set scheduler-locking on
#run 

#info variables  // to list all global and static variable names
#info locals     // to list all local variables of current stack frame (name and values), including static variables in that function
#info args       // to list arguments of the current stack frame (name and value)

