#!/bin/sh

printf " == LOADING THE MODULES ==\n"
printf "\tmodule load gcc/8.1\n"
printf "\tmodule load intel_libomp_oss\n"
module load gcc/8.1
module load intel_libomp_oss

make


