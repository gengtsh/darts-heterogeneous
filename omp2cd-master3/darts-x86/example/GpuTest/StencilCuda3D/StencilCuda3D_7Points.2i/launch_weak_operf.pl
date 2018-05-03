#!/usr/bin/perl

use warnings;
use strict;
use Cwd;

my $currpath = Cwd::abs_path();
#print "$currpath\n";
my $opreport = 'opreport';
my $results = 'results';
my $opreportpath ="$currpath/$opreport";
my $resultspath ="$currpath/$results";
#print "$opreportpath\n";
unless(-d $opreportpath){
    mkdir $opreportpath, 0755;
}
unless(-d $resultspath){
    mkdir $resultspath, 0755;
}
my $errors=0;

while ($_ = glob("$opreportpath/*.*")){
    next if -d $_;
#    print "delete $_\n";
    unlink($_) or ++$errors, warn("can't remove $_: $!");
}

while ($_ = glob("$resultspath/*.*")){
    next if -d $_;
#    print "delete $_\n";
    unlink($_) or ++$errors, warn("can't remove $_: $!");
}
exit(1) if $errors;

my @threads = (
    [31,1],
#	[15,2]#,
#   [11,4],
#   [05,8]
);

my $n_iter = 30;
my $n_reps = 10;
my $its = 1;

my @sz_zxy = (
    [50,200,200],
    [100,200,200],
    [200,200,200],
    [800,400,400],
    [200,800,800],
    [400,800,800],
    [800,800,800],
    [800,1000,1000],
    [1000,1000,1000]


);
#my $filename = shift @ARGV;
#$filename = "output" unless ($filename);
#my @kernels = qw/StencilCudaCpu37 StencilCudaGpu37 StencilCudaHybrid37/;  
my @kernels = qw/StencilCudaCpu37/;  
my @outputs = ();
my $sz_start = 1000;
my $sz_end = 51000;
my $sz_step = 2000;

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
		
		$ENV{'DARTS_AFFINITY'}="0-31";
        for my $kernel (@kernels) {
            my $ker = './' . $kernel;
            if (-e $ker) {
                my $fh_operf= "operf_${kernel}_${i}_${j}_${k}_${n_cu}_${n_su}.txt";
                open my $fh, '>>', "$resultspath/${kernel}_${i}_${j}_${k}_${n_cu}_${n_su}.txt" or die $!;
                for (my $idx = 0; $idx < $its; ++$idx) {
                    print $fh  `numactl --interleave=0-1 $ker $i $j $k $n_iter $n_reps `;
                    system `operf -e CPU_CLK_UNHALTED:6000 -e INST_RETIRED:6000 -e LLC_MISSES:6000 -e LLC_REFS:6000  -e l2_rqsts:200000 -e l1d:2000000 -e cpl_cycles:2000000 -e icache:200000 -e br_inst_exec:200000:0xff -e br_misp_exec:200000:0xff -e resource_stalls:2000000 -e fp_assist:100000 -e hw_interrupts:100000 -e l2_trans:200000 -e mem_uops_retired:2000000  -t $ker $i $j $k $n_iter $n_reps`;
                    #system `opreport -m tid -g -l -s sample > $opreportpath/$fh_operf`;
                    system `opreport -m tid -x -a -g -l -s sample > $opreportpath/$fh_operf`;
                
                }
            } else {
                print "Couldn't find kernel named $ker\n";
            }
        }
    }
}



