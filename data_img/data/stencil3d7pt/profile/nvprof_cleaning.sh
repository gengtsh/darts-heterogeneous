#!/bin/bash

declare -a Machines=( "ccsl" )

declare -a Traces=("cuDeviceGetAttribute" "cudaSetupArgument" "cudaMallocHost" "cudaGetDeviceCount" "cudaGetDeviceProperties" "cudaMemGetInfo" "cudaMalloc" "cudaMemcpy" "CUDA memcpy HtoD" "CUDA memcpy DtoH" "cudaFree" "cudaLaunch (gpu_stencil37_hack1_cp_slices" "gpu_stencil37_hack1_cp_slice" "cudaLaunch (gpu_stencil37_hack1_cp_rows" "gpu_stencil37_hack1_cp_rows" "cudaLaunch (gpu_stencil37_hack1_cp_cols" "gpu_stencil37_hack1_cp_cols" "cudaLaunch (gpu_stencil37_hack2" "gpu_stencil37_hack2")

for machine in "${Machines[@]}"; do
    
    if [[ "${machine}" -eq "ccsl" ]]; then
           declare -a ZXY=( "200_200_50" "200_200_100" "200_200_200" "400_400_800" "800_800_200" "800_800_400" "800_800_800" )
        coreThreadMachine="_7_1"
    fi

    if [[ "${machine}" -eq "f3" ]]; then
        declare -a ZXY=("200_200_50" "200_200_100" "200_200_200" "400_400_800" "800_800_200" "800_800_400" "800_800_800"  "800_1000_1000" "1000_1000_1000")
        coreThreadMachine="_31_1"

    fi
                   
   echo "machine, size, metric, invocations, min, max, mean" > tempFileMetric     
   for sizesZXY in "${ZXY[@]}"; do
       echo $sizesZXY
       declare -a arrayTempMetrics
       counter=0
   	for metric in "${Traces[@]}"; do
           tempTraceValues=`grep  "${metric}" ${machine}_nvprof_traces/StencilCudaGpu37_${sizesZXY}${coreThreadMachine}.txt | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(length(d),  min(d), max(d), mean(d), sep="\n")'`
           tempTraceValues=`echo ${tempTraceValues} | sed -e 's/ /,/g'` 
           echo "${machine}, ${sizesZXY}, ${metric}, ${tempTraceValues}" >> tempFileMetric
        done
    done
    
    mv tempFileMetric ${machine}_nvprof_traces/nvprofTraces.csv     
done

