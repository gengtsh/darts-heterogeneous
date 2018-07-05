#!/bin/bash

function join { local IFS="$1"; shift; echo "$*"; }

declare -a machines=( "ccsl" )

declare -a zxy=( "200_200_50" "200_200_100" "200_200_200" "400_400_800" "800_800_200" "800_800_400" "800_800_800" )

for machine in "${machines[@]}"; do

	declare -a metrics=("CPU_CLK_UNHALTED" "INST_RETIRED" "LLC_MISSES" "LLC_REFS" "l2_rqsts" "cpl_cycles" "icache" "br_inst_exec" "br_misp_exec" "resource_stalls" "fp_assist" "l2_trans" "mem_uops_retired") 
	header=`join , ${metrics[@]}` 
    echo "${machine}, size, ${header}" > tempFileMetric 
    
    for sizesZXY in "${zxy[@]}"; do
		#counter=0
        grep "Counted" ${machine}_opreport_cpu/operf_StencilCudaCpu37_${sizesZXY}_7_1.txt > varAll

        tempSizeMetric="${machine}, ${sizesZXY},"
		for metric in "${metrics[@]}"; do
		        var=`grep -hnr "${metric}" varAll`

        #counter=$((counter + 1))

        	if [ -n "${var}" ]; then       
         
                IFS=':' read -r -a lineMetric <<< $var
        		countMetric=`echo ${var} | awk '{print $NF}'`
        
                #lineMetric=`grep Counted ${machine}_opreport_cpu/operf_StencilCudaCpu37_${zxy}_7_1.txt 
        
        		#temp=`tail -n 1 ${machine}_opreport_cpu/operf_StencilCudaCpu37_${zxy}_7_1.txt`       
        
                tempMetric=`cat ${machine}_opreport_cpu/operf_StencilCudaCpu37_${zxy}_7_1.txt | awk '{print $0}' | tail -1`
        
                tempValueMetric=`echo $tempMetric | awk -v temp=$(($lineMetric*4-2)) ' {print $((temp))*2}'` 
        
            else
                tempValueMetric=0
            fi
        
            tempSizeMetric="$tempSizeMetric $tempValueMetric,"
    
        done
        
        echo $tempSizeMetric >> tempFileMetric      
    done
    #rm -rf var temp
done

#"${apps[@]}";
