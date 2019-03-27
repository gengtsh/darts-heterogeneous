#!/usr/bin/perl

use warnings;
use strict;

my $nThreads = 3;
#my $nThreads = 32;
my $nThreads1 = $nThreads-1;
my $tb = 1;

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
my $its = 1;


#my @kernels = qw/OMP_GCC8.1 OMP_GCC8.1_IniLoop OMP_GCC8.1_IniLoop_SIMD OMP_INTEL OMP_INTEL_IniLoop OMP_INTEL_IniLoop_SIMD /; 
#my @kernels = qw/OMP  OMP_GCC8.1_InitLoop_SIMD OMPFGV2 Naive NaiveBk NaiveTPsPtr NaiveTPsPtrBk FineGrainM2V2 FineGrainM2V2Bk/; 
#my @kernels = qw/OMP  OMP_GCC8.1_InitLoop_SIMD OMPFGV2 Naive NaiveTPsPtr FineGrainM2V2/; 
#my @kernels = qw/OMP  OMP_GCC8.1_InitLoop_SIMD/; 

#my @kernels = qw/Naivelik NaiveTPsPtrlik FineGrainM2V2lik/; 
my @kernels = qw/NaiveTPsPtrlik/; 

##do('/slowbro/software/apps/lmod/7.8/lmod/lmod/init/perl');
#do('/software/apps/lmod/lmod/init/perl');
#module('load gcc/8.1');
##module('load intel_libomp_oss');
##system('gcc -v');
##system('echo $LIBIOMP5_PATH');

$ENV{'LD_LIBRARY_PATH'} = '$LD_LIBRARY_PATH:/software/apps/likwid/4.3.2/lib:/software/apps/gcc/8.1/lib64';
#print $ENV{'LD_LIBRARY_PATH'};
#system('echo $LD_LIBRARY_PATH');

my @tppolicy = qw/compact scatter/;

#compact
my $likwid_bin = '/software/apps/likwid/4.3.2/bin/likwid-perfctr';

my @thdcmp = (0, 16, 1, 17, 2, 18, 3, 19, 4, 20, 5, 21, 6, 22, 7, 23, 8, 24, 9, 25, 10, 26, 11 ,27, 12, 28, 13, 29, 14, 30, 15, 31 );
my @thdscc = (0, 8, 1, 9, 2, 10, 3, 11, 4 ,12, 5 ,13, 6, 14, 7, 15, 16, 24, 17, 25, 18, 26, 19, 27, 20, 28 ,21, 29, 22, 30, 23, 31);

my $lk = 'likwid';

my @likwid_group = qw/CACHES DATA L2 L2CACHE L3 L3CACHE MEM NUMA FLOPS_DP FLOPS_SP TLB_DATA UOPS /;


for my $tpp (@tppolicy){
    my @thdpp;

    #for my $i (1000,3000,5000) {
    for my $i (1000) {
    
        for my $thd (@threads) {
            my ($n_cu,$n_su) = @$thd;
    		my $n_total = ($n_cu+1)*$n_su;
            $ENV{'OMP_NUM_THREADS'} = $n_total;

			$ENV{'DARTS_NUM_CU'}    = $n_cu;
			$ENV{'DARTS_NUM_SU'}    = $n_su;
			for my $kernel (@kernels) {
    	        my $ker = './' . $kernel;
    	        if (-e $ker) {

                        if($tpp eq "compact" ){
                           
                            @thdpp = @thdcmp;

                            #print $ENV{'GOMP_CPU_AFFINITY'},"\n";
                            #$ENV{'GOMP_CPU_AFFINITY'}="0 16 1 17 2 18 3 19 4 20 5 21 6 22 7 23 8 24 9 25 10 26 11 27 12 28 13 29 14 30 15 31";
                            #$ENV{'DARTS_AFFINITY'}="0 16 1 17 2 18 3 19 4 20 5 21 6 22 7 23 8 24 9 25 10 26 11 27 12 28 13 29 14 30 15 31";
                        }elsif($tpp eq "scatter"){
					
                            @thdpp = @thdscc;

                            #$ENV{'GOMP_CPU_AFFINITY'}="0 8 1 9 2 10 3 11 4 12 5 13 6 14 7 15 16 24 17 25 18 26 19 27 20 28 21 29 22 30 23 31  ";
                            #$ENV{'DARTS_AFFINITY'}="0 8 1 9 2 10 3 11 4 12 5 13 6 14 7 15 16 24 17 25 18 26 19 27 20 28 21 29 22 30 23 31";
                        }
                        
                        $ENV{'GOMP_CPU_AFFINITY'}= join(" ",@thdpp);
                        $ENV{'DARTS_AFFINITY'}= join(" ",@thdpp);

                        $ENV{'OMP_PROC_BIND'}   = 'true';
                        #print "n_cu: ", $n_cu,"\n";

                        my $likthd = join(",",@thdpp[0..$n_cu]);
                        #open my $fh, '>>', "${lk}_${tpp}_${kernel}_${i}_${n_cu}_${n_su}.txt" or die $!;
                        for (my $idx = 0; $idx < $its; ++$idx) {
                            #print $fh  `numactl --interleave=0-1 $ker $i $i $n_iter $n_reps`;
                            for my $group (@likwid_group){
                                my $fh =  "${lk}_${tpp}_${kernel}_${i}_${n_cu}_${n_su}_${group}.csv";
                                #print "group: ",$group,"\n";
                                #print "likwid_bin:", $likwid_bin,"\n";
                                #print "numactl --interleave=0-1 $likwid_bin -m -M 0 -C $likthd -c $likthd  -g $group -o $fh $ker $i $i $n_iter $n_reps ", "\n";
                                
                                #`numactl --interleave=0-1 $likwid_bin -m -M 0 -C $likthd -c $likthd  -g $group -o $fh $ker $i $i $n_iter $n_reps `
                                `$likwid_bin -m -M 0 -C $likthd -c $likthd  -g $group -o $fh $ker $i $i $n_iter $n_reps `
                            }

    	                }
    	        } else {
    	            print "Couldn't find kernel named $ker\n";
    	        }
    	    }
        }
    }
}


