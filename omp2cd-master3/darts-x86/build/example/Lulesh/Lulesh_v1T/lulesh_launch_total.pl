#!/usr/bin/perl

use warnings;
use strict;

my @threads = (
#	[1,1],
#	[3,1],
#	[7,1],
#	[15,1],
	[31,1],
	[15,1],
#	[15,2]#,
#   [11,4],
#   [05,8]
);


my $n_s = '-s';
my $n_i = '-i';
my $n_t = '-t';
my $n_iter = 20;
my $n_rep = 10;
my @n_trees = (2,4,8,16,32);
#my $filename = shift @ARGV;
#$filename = "output" unless ($filename);


#my @kernels = qw/lulesh2.0 Luleshv1/;  
my @kernels = qw/lulesh2.0/;  
#my @kernels = qw/Luleshv1T/; 
my @outputs = ();
my $szInit=200;
my $szEnd=410;
my $szStep=20;
##omp
my $ompKernel = 'lulesh2.0';
for (my $s=$szInit;$s<$szEnd; $s=$s+$szStep) {
    for my $thd (@threads) {
        my ($n_cu,$n_su) = @$thd;
        $ENV{'DARTS_NUM_CU'}    = $n_cu;
        $ENV{'DARTS_NUM_SU'}    = $n_su;
        $ENV{'OMP_PROC_BIND'}   = 'true';
        $ENV{'OMP_NUM_THREADS'} = ($n_cu+1)*$n_su;
        my $n_total = ($n_cu+1)*$n_su;
		if ($n_total== 16){
			$ENV{'GOMP_CPU_AFFINITY'}="0-7 16-23"; 
		}else{
			$ENV{'GOMP_CPU_AFFINITY'}="0-31"; 
		}


		#for my $kernel (@kernels) {
            my $ker = './' . $ompKernel;
            if (-e $ker) {
                open my $fh, '>>', "${ompKernel}_${n_cu}_${n_su}_s${s}_i${n_iter}.txt" or die $!;
			
                for (my $idx = 0; $idx < $n_rep; ++$idx) {
					print $fh  `numactl --interleave=0-1 $ker $n_s $s $n_i $n_iter`;
				}
            } else {
                print "Couldn't find kernel named $ker\n";
            }
        #}
    }
}
#luleshv1T
my $dartsKernel='Luleshv1T';
for my $ntr (@n_trees){
	for (my $s=$szInit;$s<$szEnd; $s=$s+$szStep) {
	    for my $thd (@threads) {
	        my ($n_cu,$n_su) = @$thd;
	        $ENV{'DARTS_NUM_CU'}    = $n_cu;
	        $ENV{'DARTS_NUM_SU'}    = $n_su;
	
			#for my $kernel (@kernels) {
	            my $ker = './' . $dartsKernel;
	            if (-e $ker) {
	                open my $fh, '>>', "${dartsKernel}_${n_cu}_${n_su}_s${s}_i${n_iter}_t${ntr}.txt" or die $!;
				
	                for (my $idx = 0; $idx < $n_rep; ++$idx) {
						print $fh  `numactl --interleave=0-1 $ker $n_s $s $n_i $n_iter $n_t $ntr`;
					}
	            } else {
	                print "Couldn't find kernel named $ker\n";
	            }
	        #}
	    }
	}
}
#my @threads = (
#    [31,1],
#	 [15,2]#,
##   [11,4],
##   [05,8]
#);
#
#my $n_iter = 30;
#my $n_reps = 20;
#
##my $filename = shift @ARGV;
##$filename = "output" unless ($filename);
#
##   Possible kernels:
##   OldLaplaceSolver NewLaplaceSolver SeqLaplaceSolver CompleteLaplaceSolver 
##   NaivePtrLaplaceSolver OMPLaplaceSolver ReallocLaplaceSolver 
##   LaplaceSolverCDLoop NaiveLaplaceSolver
#
##my @kernels = qw/Seq OMP FineGrain InPlace Naive NaiveTPsPtr InPlaceTPs/; # Seq/; 
#my @kernels = qw/InPlaceTPs/; # Seq/; 
#my @outputs = ();
#
#for my $i (1000,2000,3000,4000,5000,6000,7000,8000,9000,10000) {
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
#                for (my $idx = 0; $idx < 10; ++$idx) {
#                    print $fh  `numactl --interleave=0-1 $ker $i $i $n_iter $n_reps`;
#                }
#            } else {
#                print "Couldn't find kernel named $ker\n";
#            }
#        }
#    }
#}



