#!/usr/bin/perl

use warnings;
use strict;

my @threads = (
#    [31,1],
	[31,1],
#	[15,2]#,
#   [11,4],
#   [05,8]
);

my $n_iter = 30;
my $n_reps = 20;
my $its = 10;

#my $filename = shift @ARGV;
#$filename = "output" unless ($filename);

#   Possible kernels:
#   OldLaplaceSolver NewLaplaceSolver SeqLaplaceSolver CompleteLaplaceSolver 
#   NaivePtrLaplaceSolver OMPLaplaceSolver ReallocLaplaceSolver 
#   LaplaceSolverCDLoop NaiveLaplaceSolver

#my @kernels = qw/Seq OMP OMPFG OMPFGV2 FineGrain FineGrainM2V2 FineGrainM2 InPlace Naive NaiveTPsPtr InPlaceTPs FineGrainTPs/; # Seq/; 
#my @kernels = qw/StencilCUDA/; # Seq/; 
#my @kernels = qw/OMP_GCC8.1 OMP_GCC8.1_IniLoop OMP_GCC8.1_IniLoop_SIMD/; # Seq/; 
my @kernels = qw/OMP_INTEL OMP_INTEL_IniLoop OMP_INTEL_IniLoop_SIMD/ ; # Seq/; 
my @outputs = ();
my $sz_start = 1000;
my $sz_end = 11000;
my $sz_step = 1000;
for (my $i= $sz_start;$i<$sz_end;$i=$i+$sz_step){
	for my $thd (@threads) {
        my ($n_cu,$n_su) = @$thd;
        $ENV{'DARTS_NUM_CU'}    = $n_cu;
        $ENV{'DARTS_NUM_SU'}    = $n_su;
        $ENV{'OMP_PROC_BIND'}   = 'true';
        $ENV{'OMP_NUM_THREADS'} = ($n_cu+1)*$n_su;
        #$ENV{'KMP_AFFIINITY'} = 'compact';
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



