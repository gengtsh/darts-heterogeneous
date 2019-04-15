#!/usr/bin/perl

use warnings;
use strict;
use Cwd;
use List::Util qw[min max];
use POSIX;

my $currpath = Cwd::abs_path();


my $sv = $ARGV[0];
if (not defined $sv){
    die "need provide server: f3 or ccsl or debian or supermicro\n";
}

if (defined $sv){
    print "the server name: '$sv'\n";
}

my $max_x;
my $max_y;
my $max_z;

if($sv eq "f3"){
    $max_x=1010;
    $max_y=1010;
    $max_z=1010;
}elsif($sv eq "ccsl"){

    $max_x=810;
    $max_y=810;
    $max_z=810;

}elsif($sv eq "debian"){
    $max_x=1010;
    $max_y=1010;
    $max_z=1010;
}elsif($sv eq "supermicro"){
    $max_x=1010;
    $max_y=1010;
    $max_z=1010;
}    



my @sz_zxy = ();
push @sz_zxy, [50,50,50];
my $sp = 100;
for (my $z=100;$z<$max_z;$z=$z+$sp){
    for (my $x=100;$x<$max_x;$x=$x+$sp){
        #for (my $y=50;$y<$max_y;$y=$y+$sp){
            my $y = $x;
            push @sz_zxy,[$z,$x,$y];
        #}
    }
}

#map{print $_->[0],",",$_->[1],",",$_->[2],"\n";}@sz_zxy;

my @threads = (
    [2,1],
#	[15,2]#,
#   [11,4],
#   [05,8]
);


my $n_iter = 2;
my $n_reps = 1; 
my $its = 1;
my $ntMax= 30;
my $ntVal;
my $basicBlock=16;

#my @sz_zxy = (
#    [50  ,50  ,50  ],
#    [50  ,100 ,100 ],
#    [50  ,200 ,200 ],
#    [50  ,400 ,400 ],
#    [50  ,800 ,800 ],
#    [50  ,1000,1000],
#    [100 ,100 ,100 ],
#    [100 ,200 ,200 ],
#    [100 ,400 ,400 ],
#    [100 ,800 ,800 ],
#    [100 ,1000,1000],
#    [200 ,100 ,100 ],
#    [200 ,200 ,200 ],
#    [200 ,400 ,400 ],
#    [200 ,800 ,800 ],
#    [200 ,1000,1000],	
#    [400 ,100 ,100 ],
#    [400 ,200 ,200 ],
#    [400 ,400 ,400 ],
#    [400 ,800 ,800 ],
#    [400 ,1000,1000],
#    [800 ,100 ,100 ],
#    [800 ,200 ,200 ],
#    [800 ,400 ,400 ],
#    [800 ,800 ,800 ],
#    [800 ,1000,1000],	
#    [1000,100 ,100 ],
#    [1000,200 ,200 ],
#    [1000,400 ,400 ],
#    [1000,800 ,800 ],
#    [1000,1000,1000]
#
#);
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
                for(my $IsStream=0;$IsStream<2;++$IsStream){
				    for(my $m=1;$m<$n_iter;++$m){
                        if($IsStream==0){
                            my $nt = 0;
                            #open my $fh, '>>', "${kernel}_${i}_${j}_${k}_${n_cu}_${n_su}_${m}_${n_reps}_${IsStream}_${nt}.txt" or die $!;
				    	    open my $fh, '>>', "${kernel}_${i}_${j}_${k}_${m}_${n_reps}_${IsStream}_${nt}.txt" or die $!;
				    	    for (my $idx = 0; $idx < $its; ++$idx) {
				    		    print $fh  `numactl --interleave=0-1 $ker $i $j $k $m $n_reps $IsStream $nt `;
                                #system `nvprof --print-gpu-trace --print-api-trace --csv -u ns $ker $i $j $k $m $n_reps $IsStream $nt &> ${nvprof}_${kernel}_${i}_${j}_${k}_${n_cu}_${n_su}_${m}_${n_reps}_${IsStream}_$nt.txt`; 
				    		    system `nvprof --print-gpu-trace --print-api-trace --csv -u ns $ker $i $j $k $m $n_reps $IsStream $nt &> ${nvprof}_${kernel}_${i}_${j}_${k}_${m}_${n_reps}_${IsStream}_${nt}.txt`; 
	                        }
                        }elsif($IsStream==1){
                            $ntVal = min(ceil(max($i,$j,$k)/$basicBlock),$ntMax);
                            for(my $nt=1;$nt<$ntVal+1;++$nt){
                                #open my $fh, '>>', "${kernel}_${i}_${j}_${k}_${n_cu}_${n_su}_${m}_${n_reps}_${IsStream}_${nt}.txt" or die $!;
				    	        open my $fh, '>>', "${kernel}_${i}_${j}_${k}_${m}_${n_reps}_${IsStream}_${nt}.txt" or die $!;
				    	        for (my $idx = 0; $idx < $its; ++$idx) {
				    		        print $fh  `numactl --interleave=0-1 $ker $i $j $k $m $n_reps `;
                                    #system `nvprof --print-gpu-trace --print-api-trace --csv -u ns $ker $i $j $k $m $n_reps $IsStream $nt &> ${nvprof}_${kernel}_${i}_${j}_${k}_${n_cu}_${n_su}_${m}_${n_reps}_${IsStream}_$nt.txt`; 
				    		        system `nvprof --print-gpu-trace --print-api-trace --csv -u ns $ker $i $j $k $m $n_reps $IsStream $nt &> ${nvprof}_${kernel}_${i}_${j}_${k}_${m}_${n_reps}_${IsStream}_${nt}.txt`; 
	                            }
                            }
				    	}
				    }
                }
            } else {
                print "Couldn't find kernel named $ker\n";
            }
        }
    }
}



