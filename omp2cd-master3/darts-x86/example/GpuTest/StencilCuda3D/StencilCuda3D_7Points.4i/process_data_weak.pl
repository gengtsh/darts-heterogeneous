#!/usr/bin/perl

use warnings;
use strict;

use Cwd;

#use List::MoreUtils qw(first_index);
use List::Util qw[min max];
use POSIX;

my $currpath = Cwd::abs_path();
print "$currpath\n";

#my @servers = qw/debian f3 supermicro  ccsl/;
my @servers = qw/f3/;
#my @gputest = qw/gpu_streams gpu_streamNot/;


#my @vals = (1 .. 10);
##my @vals = (1,3,5);

#my @kernels = qw/StencilCudaCpu StencilCudaGpu StencilCudaHybrid3/;
#my @kernels = qw/StencilSeq37 StencilCUDA37 StencilCudaCpu37 StencilCudaGpu37 StencilCudaHybrid37/;
#my @kernels = qw/StencilCudaCpu37 StencilCudaGpu37/;
my @kernels = qw/StencilCudaGpu37/;



#my @stats = ([@kernels]);

my @threads=(
	[2,1],
#	[15,2]
	);

	
sub func_find_sz_max_xyz{
	my $sv = $_[0];
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
	return ($max_x,$max_y,$max_z);
}
	
sub func_calc_sz_zxy{
	my ($max_x,$max_y,$max_z) = @_;
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
	return @sz_zxy;
}


#my @sz_zxy = (
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
	
	
my $tc=31;
my $tc1=$tc-1;
my $n_rept=1;
my $startPos=8;
my $endPos = $n_rept+$startPos-1;


my $n_cu ;
my $n_su ;
my $n_reps = 1; 
my $n_iter = 2; 

my $ntMax= 30;
my $ntVal;
my $basicBlock=16;

my %hash_table;


### remove repeat name
#my %TPnamesHash;
#my @TPnamesDiff = grep{!$TPnamesHash{$_}++} @TPnames;
#my $stream = "gpu_streams";

for my $sv(@servers){
	#my $svpath = "$currpath/$sv/$stream";
	my $svpath = "$currpath/$sv";
	chdir($svpath) or die "Can't chdir to $svpath $!";
	print "$svpath\n";
	my %hash_server;
	
	#if($sv eq "f3"){
	#	$n_cu = 31;
	#	$n_su = 1;
	#
	#}elsif($sv eq "ccsl"){
	#	$n_cu = 7;
	#	$n_su = 1;
	#}elsif($sv eq "debian"){
	#	$n_cu = 11;
	#	$n_su = 1;	
	#}elsif($sv eq "supermicro"){
	#	$n_cu = 39;
	#	$n_su = 1;	
	#} 
	my ($max_x,$max_y,$max_z) = func_find_sz_max_xyz($sv);
	my @sz_zxy = func_calc_sz_zxy($max_x,$max_y,$max_z);

	for my $kernel (@kernels) {
		my %hash_kernel;
		my $prefix = ${kernel};
		for my $sz(@sz_zxy){
			my %hash_sz;
			
			
			my ($k,$i,$j) = @$sz;
			my $xyz=$i.'_'.$j.'_'.$k;
			for(my $IsStream=0;$IsStream<2;++$IsStream){
				my %hash_IsStream;
				for(my $m=1;$m<$n_iter;++$m){
					my %hash_iter;
					
					if($IsStream==0){
                        my $nt = 0;
						my %hash_nt;
						
						my $filename = $prefix.'_'.${i}.'_'.${j}.'_'.${k}.'_'.$m.'_'.${n_reps}.'_'.${IsStream}.'_'.${nt}.'.txt';
						
						unless (-e $filename) {
							print "$filename does not exist.\n"; 
							next;
						}
						open my $fh, '<', $filename or die "Cannot open $filename: $!";
						print "Processing $filename...\n";
						while (<$fh>) {
							chomp;
							my @data = (split /\s*,\s*/,$_)[$startPos..$endPos];
							my $data1 = $data[0];
							
							#$hash_nt0{$nt}=$data1;
							$hash_iter{$nt}=$data1;

						}

						
						
					}elsif($IsStream==1){
						$ntVal = min(ceil(min($i,$j,$k)/$basicBlock),$ntMax);
						#my %hash_nt;
                        for(my $nt=1;$nt<$ntVal+1;++$nt){
							
							my $filename = $prefix.'_'.${i}.'_'.${j}.'_'.${k}.'_'.$m.'_'.${n_reps}.'_'.${IsStream}.'_'.${nt}.'.txt';
							
							unless (-e $filename) {
							print "$filename does not exist.\n"; 
							next;
							}
							open my $fh, '<', $filename or die "Cannot open $filename: $!";
							print "Processing $filename...\n";	
							
							while (<$fh>) {
								chomp;
								my @data = (split /\s*,\s*/,$_)[$startPos..$endPos];
								my $data1 = $data[0];
							
							
								$hash_iter{$nt}=$data1;
							}		
						}
						#$hash_iter{$m}= \%hash_nt;
					}
					if(%hash_iter){	
						$hash_IsStream{$m}=\%hash_iter;
						#if($IsStream==0){
						#	print $hash_IsStream{$m}->{0},"&&&&&&&&&&&&&&&&&&&&&",", m=",$m,",IsStream:",$IsStream,",nt:",0,"\n";
						#}
					}
				}
				$hash_sz{$IsStream}=\%hash_IsStream;
			}
				#print $hash_sz{0}->{1}->{0},"++++++++++++\n";
				#print $hash_sz{1}->{1}->{1},"++++++++++++\n";			
			if(%hash_sz){
				$hash_kernel{$xyz} = \%hash_sz;
				#print $hash_kernel{$sz}->{0}->{1}->{0},"++++++++++++\n";
				#print $hash_kernel{$sz}->{1}->{1}->{1},"++++++++++++\n";
				
			}		
		}
		$hash_server{$kernel}=\%hash_kernel;
	}
	$hash_table{$sv}=\%hash_server;
}


#map{print "features: ", $_,"\n"} keys %hash_features ;

###################print output##############################
my $postpose = "postpose_weak_execution";
my @output_server;

my @out_features;
push @out_features,"machine","kernel","size(x)","size(y)","size(z)","iteration","IsStream","nTile","totalExeTime";

push @output_server,[@out_features];
#map{print "output features: ", $_,"\n"} @out_features ;

for my $sv(@servers){
	my ($max_x,$max_y,$max_z) = func_find_sz_max_xyz($sv);
	my @sz_zxy = func_calc_sz_zxy($max_x,$max_y,$max_z);
	for my $kernel (@kernels) {
		
		for my $sz(@sz_zxy){

			
			my ($k,$i,$j) = @$sz;
			my $xyz=$i.'_'.$j.'_'.$k;

			next if not exists $hash_table{$sv}->{$kernel}->{$xyz};
			
		
			for(my $IsStream=0;$IsStream<2;++$IsStream){
				for(my $m=1;$m<$n_iter;++$m){
					if($IsStream==0){
						my $nt = 0;
						
						my @out_val;
						my $val;

						if (exists $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{0}){
							$val = $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{0};
							push @out_val,$sv,$kernel,$i,$j,$k,$m,$IsStream,$nt,$val;
							push @output_server,[@out_val];
						}
						
						
					
					}elsif($IsStream==1){
						$ntVal = min(ceil(max($i,$j,$k)/$basicBlock),$ntMax);
                        for(my $nt=1;$nt<$ntVal+1;++$nt){
							my @out_val;
							my $val;
							if (exists $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{$nt}){
								$val = $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{$nt};
								push @out_val,$sv,$kernel,$i,$j,$k,$m,$IsStream,$nt,$val;
								push @output_server,[@out_val];
								
							}
												
						}
					}
				}
			}
		
		}
		
	}
}

chdir($currpath) or die "Can't chdir to $currpath $!";
print "$currpath\n";

my $file_out = $postpose.'.csv';
open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
print $fh join(',',@$_) . "\n" for (@output_server);


	
#for my $thread (@threads) {
#	my @t_exe_time;
#	my @t_speedup;
#	push @t_exe_time, ["#size",@kernels];
#	push @t_speedup, ["#size",@kernels];
#	
#	my ($n_cu,$n_su) = @$thread;
#	#for (my $v=$sz_start; $v<$sz_end;$v=$v+$sz_step){
#	for my $sz (@sz_zxy){
#	    my ($z,$x,$y) = @$sz;
#		my $ratio = 0;
#		my $value = $z.'*'.$x.'*'.$y;
#		my $value1 = $x.'_'.$y.'_'.$z;
#	    my @line  = ($value);
#		my @line_s  = ($value);
#		for my $k (@kernels) {
#			my @files = glob("$k\_$value1\_*");
#			#print @files,"\n";
#			if(@files==0){
#				print "$k\_$value1\_does not exist.\n";
#				next;
#			}	
#			for my $filename (@files){
#				unless (-e $filename) {
#					print "$filename does not exist.\n"; 
#					next;				
#				}
#
#				open my $fh, '<', $filename or die "Cannot open $filename: $!";
#				print "Processing $filename...\n";
#				my @data_t;
#				while (<$fh>) {
#					chomp;
#					my @data = (split /\s*,\s*/,$_)[$startPos..$endPos];
#					my @data_s = sort{$a <=> $b}@data;
#					shift @data_s;
#					pop @data_s;
#					@data_t=(@data_t,@data_s);
#				}
#				my $avg=0;
#				$avg +=$_ for @data_t;
#				$avg /=@data_t;
#				push @line, sprintf ("%.2f",$avg);
#				push @line_s, $line[1]/sprintf ("%.2f",$avg);
#			}
#			
#	    }
#	    push @t_exe_time, [@line];
#		push @t_speedup, [@line_s];
#		
#	}
#	
#	my $file_out = "weak_execution.dat";
#	open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
#	print $fh join(',',@$_) . "\n" for (@t_exe_time);
#	
#	my $file_out_s = 'weak_speedup'.'.dat';
#	open my $fh_s, '>', $file_out_s or die "couldn't open $file_out_s: $!";
#	print $fh_s join(',',@$_) . "\n" for (@t_speedup);
#	
#}
