#!/usr/bin/perl

use warnings;
use strict;

my $nThreads = 32;
my $nThreads1 = $nThreads-1;
my $tb = 3;


my @threads = ();
for my $j($tb..$nThreads1){
	push @threads,[$j,1];
}

#for my $j(0..15){
#	push @threads,[$j,2];
#}
##my @threads = (
##	[0,1],
##	[1,1],
##	[3,1],
##	[7,1],
##	[15,1],
##	[31,1],
##	[1,2],
##	[3,2],
##	[7,2],
##	[15,2]
##
##);
##
my $n_iter = 30;
my $n_reps = 20;
my $its = 5;

#my $filename = shift @ARGV;
#$filename = "output" unless ($filename);

#   Possible kernels:
#   OldLaplaceSolver NewLaplaceSolver SeqLaplaceSolver CompleteLaplaceSolver 
#   NaivePtrLaplaceSolver OMPLaplaceSolver ReallocLaplaceSolver 
#   LaplaceSolverCDLoop NaiveLaplaceSolver

#my @kernels = qw/Seq OMP AllocOnce CDLoop FineGrain InPlace Naive NaiveTPs NaiveTPsPtr/; # Seq/; 

#my @kernels = qw/Seq OMP AllocOnce CDLoop FineGrain InPlace Naive NaiveTPsPtr/; # Seq/; 
#my @kernels = qw/Naive NaiveTPsPtr FineGrain InPlace InPlaceTPs/; # Seq/; 

#my @kernels = qw/Seq OMP Naive NaiveTPsPtr FineGrain InPlace InPlaceTPs/; # Seq/; 

#my @kernels = qw/Seq OMP OMPFG OMPFGV2 FineGrain FineGrainM2V2 FineGrainM2 InPlace Naive NaiveTPsPtr InPlaceTPs FineGrainTPs/; # Seq/; 
#my @kernels = qw/Seq OMP OMPFG OMPFGV2 FineGrainM2V2 FineGrainM2 InPlace Naive NaiveTPsPtr InPlaceTPs FineGrainTPs/; # Seq/; 
#my @kernels = qw/FineGrainM2 /; # Seq/; 
#my @kernels = qw/StencilCudaHybrid3 StencilCudaGpu StencilCudaCpu /; # Seq/; 
#my @kernels = qw/Seq OMP/; # Seq/; 

#my @kernels = qw/SeqLaplaceSolver     OMPLaplaceSolver 
#                 ReallocLaplaceSolver LaplaceSolverCDLoop 
#                 NaiveLaplaceSolver/;

#my @kernels = qw/StencilCudaHybrid3 StencilCudaGpu StencilCudaCpu/; # Seq/; 
#my @kernels = qw/StencilCudaCpu StencilCudaHybrid3  /; # Seq/; 
#my @outputs = ();
#
#for my $i (1000,3000,5000,25000) {
##for my $i (1000,2000,3000,4000,5000,6000,7000,8000,9000,10000) {
#    for my $thd (@threads) {
#        my ($n_cu,$n_su) = @$thd;
#        $ENV{'DARTS_NUM_CU'}    = $n_cu;
#        $ENV{'DARTS_NUM_SU'}    = $n_su;
#        $ENV{'OMP_PROC_BIND'}   = 'true';
#        $ENV{'OMP_NUM_THREADS'} = ($n_cu+1)*$n_su;
#        for my $kernel (@kernels) {
#            my $ker = './' . $kernel;
#            if (-e $ker) {
#                open my $fh, '>>', "${kernel}_${i}_${n_cu}_${n_su}.txt" or die $!;
#                for (my $idx = 0; $idx < $its; ++$idx) {
#                    print $fh  `numactl --interleave=0-1 $ker $i $i $n_iter $n_reps`;
#                }
#            } else {
#                print "Couldn't find kernel named $ker\n";
#            }
#        }
#    }
#}
#

my @kernels = qw/OMP_GCC8.1 OMP_GCC8.1_IniLoop OMP_GCC8.1_IniLoop_SIMD OMP_INTEL OMP_INTEL_IniLoop OMP_INTEL_IniLoop_SIMD /; 

do('/software/apps/lmod/lmod/init/perl');
module('load gcc/8.1');
module('load intel_libomp_oss');
#system('gcc -v');
#system('echo $LIBIOMP5_PATH');

my @tppolicy = qw/compact scatter/;

#compact

for my $tpp (@tppolicy){

    for my $i (1000,3000,5000) {
    
        for my $thd (@threads) {
            my ($n_cu,$n_su) = @$thd;
    		my $n_total = ($n_cu+1)*$n_su;
            $ENV{'OMP_NUM_THREADS'} = $n_total;
            for my $kernel (@kernels) {
    	        my $ker = './' . $kernel;
    	        if (-e $ker) {
                    if ($kernel =~/GCC/){
                        if($tpp eq "compact" ){
                            $ENV{'GOMP_CPU_AFFINITY'}="0 16 1 17 2 18 3 19 4 20 5 21 6 22 7 23 8 24 9 25 10 26 11 27 12 28 13 29 14 30 15 31";
                        }elsif($tpp eq "scatter"){
                            $ENV{'GOMP_CPU_AFFINITY'}="0 8 1 9 2 10 3 11 4 12 5 13 6 14 7 15 16 24 17 25 18 26 19 27 20 28 21 29 22 30 23 31";
                        }
                        
                        $ENV{'OMP_PROC_BIND'}   = 'true';
                        delete $ENV{'KMP_AFFINITY'};
                        #print "GCC:",$ENV{'GOMP_CPU_AFFINITY'},",",$ENV{KMP_AFFINITY},",",$ENV{'OMP_PROC_BINE'},"\n";

                    }elsif($kernel =~/INTEL/){
                        $ENV{'KMP_AFFINITY'}= $tpp;
                        delete $ENV{'GOMP_CPU_AFFINITY'};
                        delete $ENV{'OMP_PROC_BIND'};
                        #print "INTEL:",$ENV{'GOMP_CPU_AFFINITY'},",",$ENV{KMP_AFFINITY},",",$ENV{'OMP_PROC_BINE'},"\n";
                    }
    	            open my $fh, '>>', "${tpp}_${kernel}_${i}_${n_cu}_${n_su}.txt" or die $!;
    	            for (my $idx = 0; $idx < $its; ++$idx) {
    	                print $fh  `numactl --interleave=0-1 $ker $i $i $n_iter $n_reps`;
    	            }
    	        } else {
    	            print "Couldn't find kernel named $ker\n";
    	        }
    	    }
        }
    }
}


