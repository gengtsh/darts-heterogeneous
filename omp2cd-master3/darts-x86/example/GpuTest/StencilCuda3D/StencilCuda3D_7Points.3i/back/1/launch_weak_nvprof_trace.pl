#!/usr/bin/perl

use warnings;
use strict;
use Cwd;

my $currpath = Cwd::abs_path();


my @threads = (
    [31,1],
#	[15,2]#,
#   [11,4],
#   [05,8]
);


my $n_iter = 10;
my $n_reps = 5; 
my $its = 1;


my @sz_zxy = (
    [50,100,100],
    [50,200,200],
    [100,100,100],
    [100,200,200],
    [200,200,200],
    [50 ,400,400],
    [100,400,400],
    [200,400,400],
    [400,400,400],
    [800,400,400],
    [50 ,800,800],
    [100,800,800],
    [200,800,800],
    [400,800,800],
    [800,800,800],
    [50 ,1000,1000],
    [100,1000,1000],
    [200,1000,1000],
    [400,1000,1000],
    [800,1000,1000],
    [1000,1000,1000]


);
#my $filename = shift @ARGV;
#$filename = "output" unless ($filename);

#   Possible kernels:
#   OldLaplaceSolver NewLaplaceSolver SeqLaplaceSolver CompleteLaplaceSolver 
#   NaivePtrLaplaceSolver OMPLaplaceSolver ReallocLaplaceSolver 
#   LaplaceSolverCDLoop NaiveLaplaceSolver

#my @kernels = qw/Seq OMP OMPFG OMPFGV2 FineGrain FineGrainM2V2 FineGrainM2 InPlace Naive NaiveTPsPtr InPlaceTPs FineGrainTPs/; # Seq/; 
my @kernels = qw/StencilCudaGpu37 /; # Seq/; 
#my @kernels = qw/StencilCudaHybrid3 /; # Seq/; 
my @outputs = ();
my $sz_start = 1000;
my $sz_end = 51000;
my $sz_step = 2000;

my $nvprof = "nvprof_trace";
#for (my $i= $sz_start;$i<$sz_end;$i=$i+$sz_step){
for my $sz(@sz_zxy){
    my ($k,$i,$j) = @$sz;
    for my $thd (@threads) {
        my ($n_cu,$n_su) = @$thd;
        $ENV{'DARTS_NUM_CU'}    = $n_cu;
        $ENV{'DARTS_NUM_SU'}    = $n_su;
        $ENV{'OMP_PROC_BIND'}   = 'true';
        $ENV{'OMP_NUM_THREADS'} = ($n_cu+1)*$n_su;
		my $n_total = ($n_cu+1)*$n_su;
        #if ($n_total <= 16){
        #    #$ENV{'GOMP_CPU_AFFINITY'}="0-7 16-23";
		#	$ENV{'DARTS_AFFINITY'}="0-7 16-23";
        #}else{
        #    #$ENV{'GOMP_CPU_AFFINITY'}="0-31";
		#	$ENV{'DARTS_AFFINITY'}="0-31";
        #}
		$ENV{'DARTS_AFFINITY'}="0-${n_cu}";
        for my $kernel (@kernels) {
            my $ker = './' . $kernel;
            if (-e $ker) {
				for(my $m=1;$m<$n_iter;++$m){
					open my $fh, '>>', "${kernel}_${i}_${j}_${k}_${n_cu}_${n_su}_${n_reps}_${m}.txt" or die $!;
					for (my $idx = 0; $idx < $its; ++$idx) {
						print $fh  `numactl --interleave=0-1 $ker $i $j $k $m $n_reps `;
						system `nvprof --print-gpu-trace --print-api-trace --csv -u ns $ker $i $j $k $m $n_reps 2> ${nvprof}_${kernel}_${i}_${j}_${k}_${n_cu}_${n_su}_${n_reps}_${m}.txt`; 
	
					}
				}
            } else {
                print "Couldn't find kernel named $ker\n";
            }
        }
    }
}



