#!/bin/bash



grep -c "Counted" ccsl_opreport_cpu/operf_StencilCudaCpu37_200_200_50_7_1.txt  > var

tail -n +$(($( cat var)+ 5)) < ccsl_opreport_cpu/operf_StencilCudaCpu37_200_200_100_7_1.txt > temp

cat temp | awk '{print $0}'

rm var temp
