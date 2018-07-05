#!/bin/bash

function join { local IFS="$1"; shift; echo "$*"; }

declare -a Machines=( "f3" )

for machine in "${Machines[@]}"; do
    
    if [[ "${machine}" -eq "ccsl" ]]; then

	    declare -a Metrics=("CPU_CLK_UNHALTED" "INST_RETIRED" "LLC_MISSES" "LLC_REFS" "l2_rqsts" "cpl_cycles" "icache" "br_inst_exec" "br_misp_exec" "resource_stalls" "fp_assist" "l2_trans" "mem_uops_retired") 
        declare -a ZXY=( "200_200_50" "200_200_100" "200_200_200" "400_400_800" "800_800_200" "800_800_400" "800_800_800" )
        coreThreadMachine="_7_1"
    fi

    if [[ "${machine}" -eq "f3" ]]; then
        declare -a Metrics=("CPU_CLK_UNHALTED" "INST_RETIRED" "LLC_MISSES" "LLC_REFS" "l2_rqsts" "l1d" "cpl_cycles" "icache" "br_inst_exec" "br_misp_exec" "resource_stalls" "fp_assist" "hw_interrupts" "l2_trans" "mem_uops_retired")
        declare -a ZXY=("200_200_50" "200_200_100" "200_200_200" "400_400_800" "800_800_200" "800_800_400" "800_800_800"  "800_1000_1000" "1000_1000_1000")
        coreThreadMachine="_31_1"

    fi
                    
	header=`join , ${Metrics[@]}` 
    echo "machine, size, ${header}" > tempFileMetric 
    
    for sizesZXY in "${ZXY[@]}"; do
		#counter=0
        grep "Counted" ${machine}_opreport_cpu/operf_StencilCudaCpu37_${sizesZXY}${coreThreadMachine}.txt > varAll

        tempSizeMetric="${machine}, ${sizesZXY},"
        tempMetric=`cat ${machine}_opreport_cpu/operf_StencilCudaCpu37_${sizesZXY}${coreThreadMachine}.txt | awk  '{print $0}' | tail -1`
		for metric in "${Metrics[@]}"; do
		        var=`grep -hnr "${metric}" varAll`

        	if [ -n "${var}" ]; then       
         
                IFS=':' read -r -a lineMetric <<< $var
        		countMetric=`echo ${var} | awk '{print $NF}'`
                tempValueMetric=`echo $tempMetric | awk -v temp=$(($lineMetric*4-2)) ' {print $((temp))}'`
                tempValueMetric=`echo $((${tempValueMetric} * ${countMetric}))`
                var=""
        
            else
                tempValueMetric=0
            fi
        
            tempSizeMetric="$tempSizeMetric $tempValueMetric,"
        done
        
        echo $tempSizeMetric >> tempFileMetric      
    done
    mv tempFileMetric ${machine}_opreport_cpu/operfMetric.csv
    rm -rf varAll
done

#"${apps[@]}";
