#!/usr/bin/perl

use warnings;
use strict;


my @threads = ();
for my $j(1..31){
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
	# my @kernels = qw/stencil.darts.exe stencil.omp.exe/; # Seq/; 
	my @kernels = qw/stencil.darts.exe/; # Seq/; 
#my @kernels = qw/Seq OMP/; # Seq/; 

#my @kernels = qw/SeqLaplaceSolver     OMPLaplaceSolver 
#                 ReallocLaplaceSolver LaplaceSolverCDLoop 
#                 NaiveLaplaceSolver/;

my @outputs = ();

#for my $i (1000) {
for my $i (1000,3000,5000) {
#for my $i (1000,2000,3000,4000,5000,6000,7000,8000,9000,10000) {
    # for my $thd (@threads) {
    for my $thd (4, 8, 10, 20) {

        # my ($n_cu,$n_su) = @$thd;
		# my $tps = 2;
		my $tps = $thd;
		$ENV{'DARTS_TPPOLICY'}="4";
		$ENV{'DARTS_MCPOLICY'}="2";
		$ENV{'DARTS_TPSLEEP'}="16";
		$ENV{'DARTS_MCSLEEP'}="16384";
		$ENV{'DARTS_AFFINITY'}="0-19";
		$ENV{'DARTS_NUMTPS'}=$tps;
		$ENV{'OMP_NUM_THREADS'}=$thd;
		$ENV{'DARTS_TPLOOP_TH0'}="4";
		$ENV{'DARTS_TPLOOP_TH1'}="16384";
		$ENV{'OMP_SCHEDULE'}=qw/static/;
        #$ENV{'DARTS_NUM_CU'}    = $n_cu;
        #$ENV{'DARTS_NUM_SU'}    = $n_su;
        #$ENV{'OMP_PROC_BIND'}   = 'true';
        #$ENV{'OMP_NUM_THREADS'} = ($n_cu+1)*$n_su;
        for my $kernel (@kernels) {
            my $ker = './' . $kernel;
            if (-e $ker) {
                # open my $fh, '>>', "${kernel}_${i}_${n_cu}_${n_su}.txt" or die $!;
                open my $fh, '>>', "${kernel}_i_${i}_threads_${thd}_tps_${tps}.txt" or die $!;
                for (my $idx = 0; $idx < 10; ++$idx) {
					print "i = $i, threads = $thd\n";
					print "numactl --interleave=0-1 $ker $i $i $n_iter $n_reps\n";
                    print $fh  `numactl --interleave=0-1 $ker $i $i $n_iter $n_reps`;
                }
            } else {
                print "Couldn't find kernel named $ker\n";
            }
        }
    }
}



