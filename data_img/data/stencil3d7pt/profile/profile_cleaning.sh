#!/bin/bash

declare -a arr=( "ccsl" )


declare -a zxy=( "200_200_50" "200_200_100" "200_200_200" "400_400_800" "800_800_200" "800_800_400" "800_800_800" )

for machine in "${arr[@]}"; do
     

    grep -c "Counted" ${machine}_opreport_cpu/operf_StencilCudaCpu37_200_200_50_7_1.txt  > var

    tail -n +$(($( cat var)+ 5)) < ${machine}_opreport_cpu/operf_StencilCudaCpu37_200_200_100_7_1.txt > temp

    cat temp | awk '{print $1}' | tail -1 

    rm -rf var temp
done

#"${apps[@]}";
