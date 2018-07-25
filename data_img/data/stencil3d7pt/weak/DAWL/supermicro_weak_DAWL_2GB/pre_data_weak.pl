#!/usr/bin/perl

use warnings;
use strict;

#my @vals = (1 .. 10);
##my @vals = (1,3,5);

#my @kernels = qw/StencilCudaCpu StencilCudaGpu StencilCudaHybrid3/;
my @kernels = qw/StencilSeq37 StencilCUDA37 StencilCudaCpu37 StencilCudaGpu37 StencilCudaHybrid37/;




#my @stats = ([@kernels]);

my @threads=(
	[39,1],
#	[15,2]
	);
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
	
	
my $tc=40;
my $tc1=$tc-1;
my $n_rept=10;
my $startPos=8;
my $endPos = $n_rept+$startPos-1;
my $sz_start = 1000;
my $sz_end = 50000;	
my $sz_step = 2000; 


	
for my $thread (@threads) {
	my @t_exe_time;
	my @t_speedup;
	push @t_exe_time, ["#size",@kernels];
	push @t_speedup, ["#size",@kernels];
	
	my ($n_cu,$n_su) = @$thread;
	#for (my $v=$sz_start; $v<$sz_end;$v=$v+$sz_step){
	for my $sz (@sz_zxy){
	    my ($z,$x,$y) = @$sz;
		my $ratio = 0;
		my $value = $z.'*'.$x.'*'.$y;
	    my @line  = ($value);
		my @line_s  = ($value);
	    for my $k (@kernels) {
	        my $filename = $k. '_'.$x.'_'.$y. '_'. $z. '_' . $n_cu . '_' . $n_su . '.txt';
	        unless (-e $filename) {
	        	print "$filename does not exist.\n"; 
	        	next;
	        }
	        open my $fh, '<', $filename or die "Cannot open $filename: $!";
	        print "Processing $filename...\n";
			my @data_t;
			while (<$fh>) {
				chomp;
				my @data = (split /\s*,\s*/,$_)[$startPos..$endPos];
				my @data_s = sort{$a <=> $b}@data;
				shift @data_s;
				pop @data_s;
				@data_t=(@data_t,@data_s);
			}
			my $avg=0;
			$avg +=$_ for @data_t;
			$avg /=@data_t;
			push @line, sprintf ("%.2f",$avg);
			push @line_s, $line[1]/sprintf ("%.2f",$avg);
			
	    }
	    push @t_exe_time, [@line];
		push @t_speedup, [@line_s];
		
	}
	
	my $file_out = $n_cu.'_'.$n_su.'_'."weak_execution.dat";
	open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
	print $fh join(',',@$_) . "\n" for (@t_exe_time);
	
	my $file_out_s = $n_cu.'_'.$n_su.'_'.'weak_speedup'.'.dat';
	open my $fh_s, '>', $file_out_s or die "couldn't open $file_out_s: $!";
	print $fh_s join(',',@$_) . "\n" for (@t_speedup);
	
}
