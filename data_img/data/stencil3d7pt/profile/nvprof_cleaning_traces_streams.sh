#!/bin/bash

stream_version="streams"

declare -a Machines=( "ccsl" "f3" "supermicro" "debian" )

declare -a Traces=("cuDeviceGetAttribute" "cudaSetupArgument" "cudaMallocHost" "cudaGetDeviceCount" "cudaGetDeviceProperties" "cudaMemGetInfo" "cudaLaunch" "cudaMalloc" "cudaMemcpy" "CUDA memcpy HtoD" "CUDA memcpy DtoH" "cudaFree" "cudaLaunch (gpu_stencil37_hack2_cp_slices" "\"gpu_stencil37_hack2_cp_slices" "cudaLaunch (gpu_stencil37_hack2_cp_rows" "\"gpu_stencil37_hack2_cp_rows" "cudaLaunch (gpu_stencil37_hack2_cp_cols" "\"gpu_stencil37_hack2_cp_cols" "cudaLaunch (gpu_stencil37_hack2" "\"gpu_stencil37_hack2")

declare -a total_times_traces=("total_CUDA memcpy HtoD" "total_CUDA memcpy DtoH" "total_kernel_execution_time" "total_CUDA_functions" "total_running_time")

declare -a ZXY=("50_100_100" "50_200_200" "100_100_100" "100_200_200" "200_200_200" "50_400_400"  "100_400_400" "200_400_400" "400_400_400" "800_400_400" "50_800_800" "100_800_800" "200_800_800" "400_800_800" "800_800_800" "50_1000_1000" "100_1000_1000" "200_1000_1000" "400_1000_1000" "800_1000_1000" "1000_1000_1000")
#declare -a ZXY=("50_100_100" "50_200_200" )

for machine in "${Machines[@]}"; do
    
    if [[ "${machine}" == "ccsl" ]]; then
        coreThreadMachine="7_1_5"
    fi

    if [[ "${machine}" == "f3" ]]; then
        coreThreadMachine="31_1_5_"
    fi

    if [[ "${machine}" == "supermicro" ]]; then
        coreThreadMachine="39_1_5"
    fi

    if [[ "${machine}" == "debian" ]]; then
        coreThreadMachine="11_1_5"
    fi
                   
    echo "machine,size,iteration,metric,invocations,min,max,mean,sum" > tempFileMetric     
    
    for sizesZXY in "${ZXY[@]}"; do
        IFS='_' read -ra ADDR <<< "$sizesZXY"
        sizesYXZ="${ADDR[2]}_${ADDR[1]}_${ADDR[0]}"            
        echo "${sizesYXZ}"
        for iter in `seq 1 9`; do
       	    for metric in "${Traces[@]}"; do

               tempTraceValues=`grep  "${metric}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(length(d),  min(d), max(d), mean(d), sum(d), sep=",")'`
               tempTraceValues=`echo ${tempTraceValues} | sed -e 's/ /,/g'` 

               echo "${machine},${sizesYXZ},${iter},${metric},${tempTraceValues}" >> tempFileMetric
            done

       	    for metric_total in "${total_times_traces[@]}"; do
                if [[ "${metric_total}" == "total_CUDA memcpy HtoD" ]]; then
                        tempTraceValues=`grep  "CUDA memcpy HtoD" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(length(d),  min(d), max(d), mean(d), sum(d), sep=",")'`
                fi

                    if [[ "${metric_total}" == "total_CUDA memcpy DtoH" ]]; then
                        tempTraceValues=`grep  "CUDA memcpy DtoH" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(length(d),  min(d), max(d), mean(d), sum(d), sep=",")'`
                    fi  

                    if [[ "${metric_total}" == "total_kernel_execution_time" ]]; then
                        K1=`grep  "${Traces[13]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                        K2=`grep  "${Traces[15]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                        K3=`grep  "${Traces[17]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                        K4=`grep  "${Traces[19]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                        total_sum=`echo $K1 $K2 $K3 $K4 | awk '{print $1 + $2 + $3 + $4}'`
                        tempTraceValues="0, 0, 0, 0, $total_sum"
					fi
                        
					if [[ "${metric_total}" == "total_CUDA_functions" ]]; then
                            K0=`grep  "${Traces[0]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                            K1=`grep  "${Traces[1]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                            K2=`grep  "${Traces[2]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                            K3=`grep  "${Traces[3]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                            K4=`grep  "${Traces[4]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                           K5=`grep  "${Traces[5]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                            K6=`grep  "${Traces[6]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                            K7=`grep  "${Traces[7]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                            K8=`grep  "${Traces[11]}" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | sed -e 's/,/ /g' |  awk '{print $2}' | Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'cat(sum(d), sep=",")'`
                            total_sum=`echo $K0 $K1 $K2 $K3 $K4 $K5 $K6 $K7 | awk '{print $1 + $2 + $3 + $4 $5 + $6 + $7 + $8 + $9}'`
                            tempTraceValues="0, 0, 0, 0, $total_sum"
						fi
						if [[ "${metric_total}" ==  "total_running_time" ]]; then
                            KMalloc1=`grep "cudaMalloc3D" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | head -1 | awk -F, '{print $1}'`
                            KMalloc2=`grep "cudaMalloc3D" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | head -1 | awk -F, '{print $2}'`                  
                            KFree1=`grep "cudaFree" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | tail -1 | awk -F, '{print $1}'`
                            KFree2=`grep "cudaFree" "${machine}/${stream_version}/nvprof_trace_StencilCudaGpu37_${sizesYXZ}_${coreThreadMachine}_${iter}.txt" | tail -1 | awk -F, '{print $2}'`
                            
                            total_sum=`echo $KMalloc1 $KMalloc2 $KFree1 $KFree2 | awk '{print $3 - $1 + $2 + $4}'`
                            tempTraceValues="0, 0, 0, 0, $total_sum"                            
						fi        
            
                echo "${machine},${sizesYXZ},${iter},${metric_total},${tempTraceValues}" >> tempFileMetric

            done
        done
    done   
    mv tempFileMetric ${machine}/nvprofTraces_${stream_version}.csv     
done

