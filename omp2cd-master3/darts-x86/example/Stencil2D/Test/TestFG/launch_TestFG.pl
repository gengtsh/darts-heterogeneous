#!/usr/bin/perl

use warnings;
use strict;

##my @threads = (
##    [31,1],
##	  [15,2]#,
##);

my @threads = ();
for my $j(1..31){
	push @threads,[$j,1];
}

my $n_iter = 30;
my $n_reps = 20;

#my $filename = shift @ARGV;
#$filename = "output" unless ($filename);

#   Possible kernels:
#   OldLaplaceSolver NewLaplaceSolver SeqLaplaceSolver CompleteLaplaceSolver 
#   NaivePtrLaplaceSolver OMPLaplaceSolver ReallocLaplaceSolver 
#   LaplaceSolverCDLoop NaiveLaplaceSolver

my @kernels = qw/TestFG/; # Seq/; 
my @outputs = ();

for (my $i =10;$i<110;$i=$i+5) {
	for my $thd (@threads) {
        my ($n_cu,$n_su) = @$thd;
        $ENV{'DARTS_NUM_CU'}    = $n_cu;
        $ENV{'DARTS_NUM_SU'}    = $n_su;
        $ENV{'OMP_PROC_BIND'}   = 'true';
        $ENV{'OMP_NUM_THREADS'} = ($n_cu+1)*$n_su;
        for my $kernel (@kernels) {
            my $ker = './' . $kernel;
            if (-e $ker) {
                open my $fh, '>>', "${kernel}_${i}_${n_cu}_${n_su}.txt" or die $!;
                for (my $idx = 0; $idx < 10; ++$idx) {
					#print $fh  `numactl --interleave=0-1 $ker $i $i $n_iter $n_reps`;
					print $fh  `numactl --interleave=0-1 $ker $i $n_reps`;
                }
            } else {
                print "Couldn't find kernel named $ker\n";
            }
        }
    }
}



