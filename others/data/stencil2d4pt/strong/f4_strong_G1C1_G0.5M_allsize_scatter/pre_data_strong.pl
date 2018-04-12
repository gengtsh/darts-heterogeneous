#!/usr/bin/perl

use warnings;
use strict;
use Cwd;


my @vals = ();
my $v_bg = 1000;
my $v_ed = 43000;
my $v_sp = 2000;
for (my $i= $v_bg;$i<$v_ed;$i = $i+$v_sp){
	push @vals, $i;
}
#@vals = (500,1000,3000,5000,25000);
my $vs = @vals;
my $vs1 = @vals+1;

#my @kernels = qw/StencilCudaCpu StencilCudaGpu StencilCudaHybrid3/;
#my @kernels = qw/StencilSEQ StencilCUDA StencilCudaCpu StencilCudaGpu StencilCudaHybrid3/;
my @kernels = qw/StencilSEQ StencilCUDA StencilCudaCpu StencilCudaHybrid3/;
my @kernels_cp = @kernels;
my @threads=();
my @threads_c=();
my $tb=3;
my $tc=32;
my $tc1=$tc-1;
my $n_rept=20;
my $startPos=8;
my $endPos = $n_rept+$startPos-1;


for my $j($tb..$tc1){
	push @threads, [$j,1];
	push @threads_c, $j;
}
#for my $j(0..15){
#	push @threads, [$j,2];
#}
#my @SU=(1,2);

#my $dir = getcwd;
#opendir (my $dh, $dir) or die "Could not open '$dir' for reading '$!'\n";;


#seq

my $k_seq = shift @kernels;
#my @files = grep(/^$k_seq/,readdir ($dh));
#my @time_seq;
my %time_seq;
for my $value (@vals){
	my @files = glob("$k_seq\_$value*");
	#print @files,"\n";
	if(@files==0){
		print "$k_seq\_$value\_does not exist.\n";
		next;
	}	
	
	for my $filename (@files){
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
		#push @v_seq, sprintf ("%.2f",$avg);
		$time_seq{$value} = sprintf ("%.2f",$avg);
		
	}

}
#print $time_seq{"1000"};

#cuda

my $k_cuda = shift @kernels;
#my @time_cuda;
my %time_cuda;
for my $value (@vals){
	my @files = glob("$k_cuda\_$value*");
	#print @files,"\n";
	for my $filename (@files){
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
		#push @v_cuda, sprintf ("%.2f",$avg);
		$time_cuda{$value}=sprintf ("%.2f",$avg);
		
	}

}


##others
#
my @SU=(1);
my $n_cu;

my %total;
for my $n_su(@SU){
	for my $value (@vals){
		my %value_kernels;
		for my $k (@kernels) {
			my @k_value;
			for my $thread (@threads) {	
				$n_cu = $thread->[0];
				my $total_cores=($n_cu+1)*$n_su;
				if($total_cores >$tc){
					next;
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
				push @k_value, sprintf ("%.2f",$avg);
				
			}
			$value_kernels{$k}=[@k_value]
		}
		#push @total,[@value_kernels];
		$total{$value}= \%value_kernels;
	}
}

my %total_exe_time;
my %total_speedup;
#execution time
for my $n_su(@SU){
	
	for my $value (@vals){
		my @v_exe_time;
		my @v_speedup;
		push @v_exe_time, ["#cores",@kernels_cp];
		push @v_speedup, ["#cores",@kernels_cp];
		#for my $thread (@threads) {
		for (my $i =0;$i < @threads; ++$i){
			my $thread = $threads[$i];
			$n_cu = $thread->[0];		
			my $total_cores=($n_cu+1)*$n_su;
			if($total_cores >$tc){
				next;
			}
			my @line;
			my @line_s;
			push @line, $total_cores,$time_seq{$value},$time_cuda{$value};
			push @line_s, $total_cores,$time_seq{$value}/$time_seq{$value},$time_seq{$value}/$time_cuda{$value};
			for my $k (@kernels) {
				push @line,$total{$value}->{$k}[$i];
				push @line_s,$time_seq{$value}/$total{$value}->{$k}[$i];
			}
			push @v_exe_time,[@line];
			push @v_speedup,[@line_s];
			
		
		}
		$total_exe_time{$value}= [@v_exe_time];
		$total_speedup{$value}= [@v_speedup];
		
		my $file_out = $value.'_'.$n_su.'_'.'strong_execution'.'.dat';
		open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
		print $fh join(',',@$_) . "\n" for (@v_exe_time);
		my $file_out_s = $value.'_'.$n_su.'_'.'strong_speedup'.'.dat';
		open my $fh_s, '>', $file_out_s or die "couldn't open $file_out_s: $!";
		print $fh_s join(',',@$_) . "\n" for (@v_speedup);			
	}
}

############################################################################################################
my %total_k;
for my $n_su(@SU){
	for my $k (@kernels) {
		my %k_value;
		for my $value (@vals){
			my @value_k;
		
			for my $thread (@threads) {	
				$n_cu = $thread->[0];
				my $total_cores=($n_cu+1)*$n_su;
				if($total_cores >$tc){
					next;
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
				push @value_k, sprintf ("%.2f",$avg);
				
			}
			$k_value{$value}=[@value_k]
		}
		#push @total,[@value_kernels];
		$total_k{$k}= \%k_value;
	}
}
$total_k{$k_seq}= \%time_seq;
$total_k{$k_cuda}= \%time_cuda;

my %total_exe_time_k;
my %total_speedup_k;
#execution time
for my $n_su(@SU){

	for my $k (@kernels_cp) {
		my @v_exe_time;
		my @v_speedup;
		push @v_exe_time, ["#cores",@vals];
		push @v_speedup, ["#cores",@vals];
		
		#for my $thread (@threads) {
		for (my $i =0;$i < @threads; ++$i){
			my $thread = $threads[$i];
			$n_cu = $thread->[0];		
			my $total_cores=($n_cu+1)*$n_su;
			if($total_cores >$tc){
				next;
			}
			my @line;
			my @line_s;
			push @line, $total_cores;
			push @line_s, $total_cores;
			
			for my $value (@vals){
				if(($k eq $k_seq)||($k eq $k_cuda)){
					push @line,$total_k{$k}->{$value};
					push @line_s,$time_seq{$value}/$total_k{$k}->{$value};				
				}else{
					push @line,$total_k{$k}->{$value}[$i];
					push @line_s,$time_seq{$value}/$total_k{$k}->{$value}[$i];				
				}

			}
			push @v_exe_time,[@line];
			push @v_speedup,[@line_s];
			
		
		}
		$total_exe_time{$k}= [@v_exe_time];
		$total_speedup{$k}= [@v_speedup];
		
		my $file_out = $k.'_'.$n_su.'_'.'strong_execution'.'.dat';
		open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
		print $fh join(',',@$_) . "\n" for (@v_exe_time);
		my $file_out_s = $k.'_'.$n_su.'_'.'strong_speedup'.'.dat';
		open my $fh_s, '>', $file_out_s or die "couldn't open $file_out_s: $!";
		print $fh_s join(',',@$_) . "\n" for (@v_speedup);			
	
	}
}

