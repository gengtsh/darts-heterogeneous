#!/usr/bin/perl

use warnings;
use strict;
use Cwd;
#use List::MoreUtils qw(first_index);
use List::Util qw[min max first];
use POSIX;

my $currpath = Cwd::abs_path();
print "$currpath\n";

#my @servers = qw/debian f3 supermicro  ccsl/;
my @servers = qw/f3/;
#my @gputest = qw/gpu_streams gpu_streamNot/;

my @types = qw/nvprof execution/;
my %hash_types=($types[0]=>"metric",$types[1]=>"totalExeTime");

my $ppose   = "postpose";
my @kernels = qw/StencilCudaGpu37/;

my %hash_table;
my %hash_metrics;

my @feature_right;
my @feature_left;


sub func_process_file2{
	my $fh = $_[0];
	my $ty = $_[1];
	
	

	my %hash_fh;
	my @aFeatureLeft;
	my @aFeatureRight;
	my $metricIndex;
	my $metricIndexMinus1;
	my $metricIndexPlus1;
	my $cnt=0;
	my @LastDataLeft=();
	my $FeatureMid;
	
	my $PorE = $hash_types{$ty};
	#print "PorE:", $PorE,"\n";
	
	while (<$fh>) {
		chomp;
		#print "line number:", $.,"\n";
		#next if ($_ =~/Profiling/ || $_ =~ /profiling/);
		my @data = (split /\s*,\s*/,$_);
		
		#print "cnt: ",$cnt,"," if ($PorE =~/metric/);
		if ($_ =~/machine/){
			
			$metricIndex =  first{$data[$_] =~ /$PorE/}0..$#data;
			$metricIndexMinus1 = $metricIndex -1;
			$metricIndexPlus1  = $metricIndex +1;
			
			@aFeatureLeft = @data[0...$metricIndexMinus1];
			@aFeatureRight = @data[$metricIndexPlus1..$#data];
			$FeatureMid = $data[$metricIndex];
			
			$hash_fh{"FeatureLeft"}=\@aFeatureLeft;
			$hash_fh{"FeatureRight"}=\@aFeatureRight;
			
			
		}else{
			my @dataLeft  = @data[0..$metricIndexMinus1];
			$metricIndexPlus1 = ($PorE =~/metric/)? ($#data-$#aFeatureRight): $metricIndexPlus1;
			my $name = join("",@data[$metricIndex..($metricIndexPlus1-1)]);
			#print "name:", $name,"\n";
			
			my $Index = ($PorE =~/metric/)?$metricIndexPlus1:$metricIndex;
			my @dataRight = @data[$Index..$#data];

			next if (@data == 0);
			if (!(@dataLeft ~~ @LastDataLeft)){
				
				#print "dataLeft:\n";
				#map {print $_,"," }@dataLeft; 
				#print "\n";
				#print "LastDataLeft:\n";
				#map {print $_,"," }@LastDataLeft;
				#print "\n";
				@LastDataLeft=@dataLeft;
				
				++$cnt;
				$hash_fh{$cnt}->{"dataLeft"}=\@dataLeft;
			}
			my $mts = ($PorE =~/metric/)? $name:$FeatureMid; 
				
			++$hash_metrics{$mts} if ($PorE =~/metric/);
			$hash_fh{$cnt}->{$mts}=\@dataRight;
			
		}		
	}
	@feature_right = @aFeatureRight if ($PorE =~/metric/);
	@feature_left = @aFeatureLeft if ($PorE =~/metric/);
	#print "cnt: ", $cnt,"\n";
	return \%hash_fh;
}



### remove repeat name
#my %TPnamesHash;
#my @TPnamesDiff = grep{!$TPnamesHash{$_}++} @TPnames;
#my $stream = "gpu_streams";

#opendir (my $currDir, $currpath) or die "can't open $currpath";
#my @file = grep{/.csv/} readdir($currDir);
#close $currDir;

my @file = glob("*.csv");
map{print "directory: ", $currpath, ", files: ", $_, "\n"} @file;
 
{
	#my $svpath = "$currpath/$sv";
	#chdir($svpath) or die "Can't chdir to $svpath $!";
	#print "$svpath\n";
	
	for my $ty (@types){
		
		my @ff = grep{/$ty/}@file;
		if (@ff == 0){
			print "$ty files does not exist.\n";
			next;
		}
		my %hash_tp;
		for my $filename (@ff){
			unless (-e $filename) {
				print "$filename does not exist.\n"; 
				next;				
			}

			open my $fh, '<', $filename or die "Cannot open $filename: $!";
			print "Processing $filename...\n";
			
			$hash_tp{$filename}= func_process_file2($fh,$ty);
		
		}
		$hash_table{$ty}=\%hash_tp;
	}
}

#map{print "metrics: ", $_,"\n"} keys %hash_metrics ;

###################print output##############################

my $postpose = "postpose_integrate";
my @out_features;
my %hash_out_data;
my @out_data;
my @aMetrics;
my @out_features_right;

for my $key (keys %hash_metrics){
	push @aMetrics,$key;
}

push @out_data,["metrics:",@aMetrics];
push @out_data,["feature:",@feature_right];

for my $am (@aMetrics){
	map {push @out_features_right,$am."_".$_} @feature_right;
}
push @out_features, @feature_left,@out_features_right,"TrueExe";

#map {print "out_features: ", $_,"\n" }@out_features;

push @out_data, [@out_features];


	
for my $sv(@servers){


	my @ff = grep{/$sv/}@file;
	if (@ff == 0){
		print "$sv server files does not exist.\n";
		next;
	}
	
	#my $ty0Idx =  first{$ff[$_] =~ /$types[0]/}0..$#ff;
	#my $ty1Idx =  first{$ff[$_] =~ /$types[1]/}0..$#ff;
	#my @fns;
	my @fn_hash;	
	for my $idx (0..$#types){
		my $fn_idx = first {$ff[$_] =~ /$types[$idx]/} 0..$#ff;
		$fn_hash[$idx]=$hash_table{$types[$idx]}->{$ff[$fn_idx]};
			
	}
	my $fn_hash0=$fn_hash[0];
	my $h_sz = keys %$fn_hash0;
	my $h_sz_data = $h_sz-2;
	#print "hash_size:", $h_sz;
	for my $sz (1..$h_sz_data){
		my @line;
		my $nv_left  = $fn_hash[0]->{$sz}->{"dataLeft"};
		my $exe_left = $fn_hash[1]->{$sz}->{"dataLeft"};
		push @line,@$nv_left;
		for my $m (@aMetrics){
			my $nv_right_t = $fn_hash[0]->{$sz}->{$m};
			push @line,@$nv_right_t;
		}
		if(@$nv_left ~~ @$exe_left){
			my $tmp = $fn_hash[1]->{$sz}->{$hash_types{$types[1]}};
			push @line,@$tmp;
		}
		push @out_data,[@line];
	}

}
#
##chdir($currpath) or die "Can't chdir to $currpath $!";
##print "$currpath\n";
##
my $file_out = $postpose.'.csv';
open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
print $fh join(',',@$_) . "\n" for (@out_data);


##
##
##for my $sv(@servers){
##	my ($max_x,$max_y,$max_z) = func_find_sz_max_xyz($sv);
##	my @sz_zxy = func_calc_sz_zxy($max_x,$max_y,$max_z);
##	for my $kernel (@kernels) {
##		
##		#my @out_kel;
##		
##	
##		#my $postfix = $postpose.'_'.$kernel;
##		for my $sz(@sz_zxy){
##		
##			my ($k,$i,$j) = @$sz;
##			my $xyz=$i.'_'.$j.'_'.$k;
##		
##			next if not exists $hash_table{$sv}->{$kernel}->{$xyz};
##		
##			for(my $IsStream=0;$IsStream<2;++$IsStream){		
##				for(my $m=1;$m<$n_iter;++$m){
##				
##					if($IsStream==0){
##						my $nt = 0;
##						for my $metric  (keys %hash_metrics){
##							my @out_val;
##							push @out_val,$sv,$kernel,$i,$j,$k,$m,$IsStream,$nt,$metric;
##							#print "metric:", $metric,"\n";
##							for my $fet (@features){
##								#print "features:", $fet,"\n";
##								my $val;
##								
##															
##								if (exists $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{$nt}->{$metric}->{$fet} ){
##									$val = $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{$nt}->{$metric}->{$fet} ;
##								}else{
##									$val = 0;
##								}
##								push @out_val, $val;
##							}
##							push @output_server,[@out_val];
##							
##						}	
##					
##					}elsif ($IsStream==1){
##						$ntVal = min(ceil(min($i,$j,$k)/$basicBlock),$ntMax);
##                        for(my $nt=1;$nt<$ntVal+1;++$nt){
##
##							for my $metric  (keys %hash_metrics){
##								my @out_val;
##								push @out_val,$sv,$kernel,$i,$j,$k,$m,$IsStream,$nt,$metric;
##								for my $fet (@features){
##									my $val;
##							
##									if (exists $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{$nt}->{$metric}->{$fet} ){
##										$val = $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{$nt}->{$metric}->{$fet};
##									}else{
##										$val = 0;
##									}
##									push @out_val, $val;
##								}
##								push @output_server,[@out_val];
##							}						
##						}
##					
##					}	
##				}
##			}
##		}
##		
##
##		
##	}
##}
##
##chdir($currpath) or die "Can't chdir to $currpath $!";
##print "$currpath\n";
##
##my $file_out = $postpose.'.csv';
##open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
##print $fh join(',',@$_) . "\n" for (@output_server);

