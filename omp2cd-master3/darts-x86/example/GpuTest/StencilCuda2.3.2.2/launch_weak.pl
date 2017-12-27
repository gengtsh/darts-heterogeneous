#!/usr/bin/perl

use warnings;
use strict;

my @threads = (
    [31,1],
#	[15,2]#,
#   [11,4],
#   [05,8]
);

my $n_iter = 30;
my $n_reps = 20;
my $its = 5;

#my $filename = shift @ARGV;
#$filename = "output" unless ($filename);

#   Possible kernels:
#   OldLaplaceSolver NewLaplaceSolver SeqLaplaceSolver CompleteLaplaceSolver 
#   NaivePtrLaplaceSolver OMPLaplaceSolver ReallocLaplaceSolver 
#   LaplaceSolverCDLoop NaiveLaplaceSolver

#my @kernels = qw/Seq OMP OMPFG OMPFGV2 FineGrain FineGrainM2V2 FineGrainM2 InPlace Naive NaiveTPsPtr InPlaceTPs FineGrainTPs/; # Seq/; 
my @kernels = qw/StencilCudaHybrid3 StencilCudaGpu StencilCudaCpu/; # Seq/; 
#my @kernels = qw/StencilCudaHybrid3 /; # Seq/; 
my @outputs = ();
my $sz_start = 1000;
my $sz_end = 51000;
my $sz_step = 2000;
for (my $i= $sz_start;$i<$sz_end;$i=$i+$sz_step){
	for my $thd (@threads) {
        my ($n_cu,$n_su) = @$thd;
        $ENV{'DARTS_NUM_CU'}    = $n_cu;
        $ENV{'DARTS_NUM_SU'}    = $n_su;
        $ENV{'OMP_PROC_BIND'}   = 'true';
        $ENV{'OMP_NUM_THREADS'} = ($n_cu+1)*$n_su;
		my $n_total = ($n_cu+1)*$n_su;
        if ($n_total== 16){
            #$ENV{'GOMP_CPU_AFFINITY'}="0-7 16-23";
			$ENV{'DARTS_AFFINITY'}="0-7 16-23";
        }else{
            #$ENV{'GOMP_CPU_AFFINITY'}="0-31";
			$ENV{'DARTS_AFFINITY'}="0-31";
        }
		
        for my $kernel (@kernels) {
            my $ker = './' . $kernel;
            if (-e $ker) {
                open my $fh, '>>', "${kernel}_${i}_${n_cu}_${n_su}.txt" or die $!;
                for (my $idx = 0; $idx < $its; ++$idx) {
                    print $fh  `numactl --interleave=0-1 $ker $i $i $n_iter $n_reps `;
                }
            } else {
                print "Couldn't find kernel named $ker\n";
            }
        }
    }
}



