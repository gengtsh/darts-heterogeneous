#!/usr/bin/perl

use warnings;
use strict;

use Cwd;

#use Cwd qw(abs_path);
my $currpath = Cwd::abs_path();
print "$currpath\n";

#my @servers = qw/f4 supermicro debian hive ccsl/;
my @servers = qw/f4/;
my $platform = 'darts';
my $sw = 'weak';
my $GC= 'G1C1';
my $type = 'hybrid';
my $sche = 'scatter';
my @gcvalue = ('G0.5M','G2000');
my @gpumem = ('2GB','4GB');



#for my $server (@servers){
#	my $fdname = $server.'_'.$sw.'_'.$GC;
#	my $path = "$currpath/$fdname";
#	chdir($path) or die "Can't chdir to $path $!";
#	print "$path\n";
#	chdir($currpath) or die "Can't chdir to $currpath $!";
#}


#my @kernels = qw/StencilCudaCpu StencilCudaGpu StencilCudaHybrid3/;
#my @kernels = qw/StencilSEQ StencilCUDA StencilCudaCpu StencilCudaGpu StencilCudaHybrid3/;
my @kernels = qw/StencilSEQ StencilCudaHybrid3/;


my @threads=(
	[31,1],
	[39,1],
	[11,1],
	[11,1],
	[8,1],
#	[15,2]
	);

	
	

my $sz_start = 1000;
my $sz_end = 50000;	
my $sz_step = 2000; 


my %exe_time;
my %speedup;	


for my $server (@servers){
	my $thread;
	if ($server eq "f4"){
		$thread = [31,1];
		$sz_end = 50000;
	}elsif($server eq "supermicro"){
		$thread = [39,1];
		$sz_end = 50000;
	}elsif($server eq "debian"){
		$thread = [11,1];
		$sz_end = 37000;
	}elsif($server eq "hive"){
		$thread = [11,1];
		$sz_end = 37000;
	}elsif($server eq "ccsl"){
		$thread = [7,1];
		$sz_end = 26000;
	}


	for my $gmem (@gpumem){
		for my $gcv (@gcvalue){
			my $fdname = $server.'_'.$sw.'_'.$platform.'_'.$type.'_'.$GC.'_'.$gcv.'_'.$gmem.'_'.$sche;
			my $path = "$currpath/$fdname";
			chdir($path) or die "Can't chdir to $path $!";	
			print "$path\n";
			
			my %t_exe_time;
			my %t_speedup;
			my %seq_exe_time;
			my @line;
			my @line_s;
			my ($n_cu,$n_su) = @$thread;
			my $tc=$n_cu+$n_su;
			my $tc1=$tc-1;
			my $n_rept=20;
			my $startPos=8;
			my $endPos = $n_rept+$startPos-1;		
			for (my $v=$sz_start; $v<$sz_end;$v=$v+$sz_step){
				my $value = $v;
				for my $k (@kernels) {
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
					if($k eq 'StencilSEQ'){
						$seq_exe_time{$value} = sprintf ("%.2f",$avg);
					}else{
						$t_exe_time{$value}=sprintf ("%.2f",$avg);
						$t_speedup{$value} = $seq_exe_time{$value}/sprintf ("%.2f",$avg);
					
					}			
				}
	
				
			}
			
			$exe_time{$gcv.'_'.$gmem} = \%t_exe_time;
			$speedup{$gcv.'_'.$gmem} = \%t_speedup;
			$exe_time{'StencilSEQ'} =\%seq_exe_time;
			
		}
		

	}

	
}
	
	
	
#output


for my $server (@servers){
	my $thread;
	if ($server eq "f4"){
		$thread = [31,1];
		$sz_end = 50000;
	}elsif($server eq "supermicro"){
		$thread = [39,1];
		$sz_end = 50000;
	}elsif($server eq "debian"){
		$thread = [11,1];
		$sz_end = 37000;
	}elsif($server eq "hive"){
		$thread = [11,1];
		$sz_end = 37000;
	}elsif($server eq "ccsl"){
		$thread = [7,1];
		$sz_end = 26000;
	}
	my ($n_cu,$n_su) = @$thread;
	my @t_exe_time;
	my @t_speedup;

	my @tmp;
	for my $gmem (@gpumem){
		for my $gcv (@gcvalue){
			push @tmp,$gcv.'_'.$gmem;
		}
	}
	push @t_exe_time, ["#size","Seq",@tmp];
	push @t_speedup, ["#size","Seq",@tmp];		
	for (my $v=$sz_start; $v<$sz_end;$v=$v+$sz_step){
		my $value = $v;
		my @line  = ($value,$exe_time{'StencilSEQ'}->{$value});
		my @line_s  = ($value,1);
		
		for my $gmem (@gpumem){
			for my $gcv (@gcvalue){
				push @line,$exe_time{$gcv.'_'.$gmem}->{$value};
				push @line_s,$speedup{$gcv.'_'.$gmem}->{$value};
			}
		}
		push @t_exe_time, [@line];
		push @t_speedup, [@line_s];		
	}

		
	
	chdir($currpath) or die "Can't chdir to $currpath $!";	
	my $file_out = $server.'_'.$n_cu.'_'.$n_su.'_'."weak_execution.dat";
	open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
	print $fh join(',',@$_) . "\n" for (@t_exe_time);
	
	my $file_out_s = $server.'_'.$n_cu.'_'.$n_su.'_'.'weak_speedup'.'.dat';
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


