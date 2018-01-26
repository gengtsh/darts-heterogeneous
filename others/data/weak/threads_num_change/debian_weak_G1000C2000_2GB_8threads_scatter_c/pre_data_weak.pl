#!/usr/bin/perl

use warnings;
use strict;

#my @vals = (1 .. 10);
##my @vals = (1,3,5);

#my @kernels = qw/StencilCudaCpu StencilCudaGpu StencilCudaHybrid3/;
my @kernels = qw/StencilSEQ StencilCUDA StencilCudaCpu StencilCudaGpu StencilCudaHybrid3/;




#my @stats = ([@kernels]);

my @threads=(
	[11,1],
#	[15,2]
	);

	
	
my $tc=12;
my $tc1=$tc-1;
my $n_rept=20;
my $startPos=8;
my $endPos = $n_rept+$startPos-1;
my $sz_start = 1000;
my $sz_end = 37000;	
my $sz_step = 2000; 


	
for my $thread (@threads) {
	my @t_exe_time;
	my @t_speedup;
	push @t_exe_time, ["#size",@kernels];
	push @t_speedup, ["#size",@kernels];
	
	my ($n_cu,$n_su) = @$thread;
	for (my $v=$sz_start; $v<$sz_end;$v=$v+$sz_step){
	    my $value = $v;
		my $ratio = 0;
	    my @line  = ($value);
		my @line_s  = ($value);
	    for my $k (@kernels) {
			if ($k eq 'StencilCudaHybrid3'){
				$n_cu = 7;
			}else{
				$n_cu = 11;
			}
	        my $filename = $k . '_' . $value . '_' . $n_cu . '_' . $n_su . '.txt';
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
	$n_cu=11;
	my $file_out = $n_cu.'_'.$n_su.'_'."weak_execution.dat";
	open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
	print $fh join(',',@$_) . "\n" for (@t_exe_time);
	
	my $file_out_s = $n_cu.'_'.$n_su.'_'.'weak_speedup'.'.dat';
	open my $fh_s, '>', $file_out_s or die "couldn't open $file_out_s: $!";
	print $fh_s join(',',@$_) . "\n" for (@t_speedup);
	
}

#for (my $v=$sz_start;$v<$sz_end;$v=$v+$sz_step){
#	push @vals, $v;
#}
#
#my $st = @kernels;
#
#for my $thread (@threads) {
#	my ($n_cu,$n_su) = @$thread;
#	my $Multi_file=$n_cu.'_'.$n_su.'_'.'weak_execution.dat';
#	unless (-e $Multi_file) {
#		print "$Multi_file does not exist.\n"; 
#		next;
#	}
#	open my $fh_multi, '<', $Multi_file or die "Cannot open $Multi_file: $!";
#	print "Processing $Multi_file...\n";
#	
#	my @tmp_multi  = map { (split /\s*,\s*/, $_) } <$fh_multi>;
#	
#	my @stats;
#	push @stats,[@size];
#	for (my $idx=1;$idx<@vals+1;++$idx){
#		my @lines;
#		push @lines, $tmp_multi[$idx*$st];
#		push @lines, 1;
#		my $sp;
#		for (my $i=2;$i<$st;++$i){
#			$sp=$tmp_multi[$st*$idx+1]/$tmp_multi[$st*$idx+$i];
#			push @lines,sprintf ("%.2f",$sp);
#		}
#		push @stats, [@lines];
#	}
#	my $file_out = $n_cu.'_'.$n_su.'_'.'weak_speedup'.'.dat';
#	open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
#	print $fh join(',',@$_) . "\n" for (@stats);
#
#
#}


