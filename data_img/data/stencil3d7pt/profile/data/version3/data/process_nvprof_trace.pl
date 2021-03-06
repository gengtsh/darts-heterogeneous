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


my $nvprof = "nvprof_trace";
my @kernels = qw/StencilCudaGpu37/;


#my @sz_zxy = (
#   [50  ,50  ,50  ],
#   [50  ,100 ,100 ],
#   [50  ,200 ,200 ],
#   [50  ,400 ,400 ],
#   [50  ,800 ,800 ],
#   [50  ,1000,1000],
#   [100 ,100 ,100 ],
#   [100 ,200 ,200 ],
#   [100 ,400 ,400 ],
#   [100 ,800 ,800 ],
#   [100 ,1000,1000],
#   [200 ,100 ,100 ],
#   [200 ,200 ,200 ],
#   [200 ,400 ,400 ],
#   [200 ,800 ,800 ],
#   [200 ,1000,1000],	
#   [400 ,100 ,100 ],
#   [400 ,200 ,200 ],
#   [400 ,400 ,400 ],
#   [400 ,800 ,800 ],
#   [400 ,1000,1000],
#   [800 ,100 ,100 ],
#   [800 ,200 ,200 ],
#   [800 ,400 ,400 ],
#   [800 ,800 ,800 ],
#   [800 ,1000,1000],	
#   [1000,100 ,100 ],
#   [1000,200 ,200 ],
#   [1000,400 ,400 ],
#   [1000,800 ,800 ],
#   [1000,1000,1000]
#
#);




my %hash_table;
my $Name = "Name";



#my $machine="ccsl";

my %hash_metrics;
my %hash_features;

my %hash_emptyfiles;

my $n_cu ;
my $n_su ;
my $n_reps = 1; 
my $n_iter = 2; 

my $ntMax= 30;
my $ntVal;
my $basicBlock=16;

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


sub func_process_file{
	my $fh = $_[0];
	
	my @properties;
	my %hash_fh;
	my $nameIndex;
	my $callNum="callNum";
	my $totalExe="totalExe";
	my $avgExe="avgExe";
	my $maxExe="maxExe";
	my $minExe="minExe";
	my $currExe; #="Duration";
	my $Start;
	my $size;
	
	my $totalThroughput = "totalThroughput";
	my $avgThroughput = "avgThroughput";
	my $maxThroughput = "maxThroughput";
	my $minThroughput = "minThroughput";
	my $maxGridx ="maxGridx";
	my $maxGridy ="maxGridy";
	my $maxGridz ="maxGridz";
	my $minGridx ="minGridx";
	my $minGridy ="minGridy";
	my $minGridz ="minGridz";
	
	my $gridx;
	my $gridy;
	my $gridz;
	my $Throughput;
	
	while (<$fh>) {
		chomp;
		#print "line number:", $.,"\n";
		next if ($_ =~/Profiling/ || $_ =~ /profiling/);
		
		next if ($_ =~/^$/ )|| ($_=~/Error|Warning|No kernels/);
		
		next if ($_ =~/^\d*\s+\*\s+\d*\s+\*\s+\d*/);
		
		next if ($_ =~/^\==/);
		
	
	
		if ($_ =~/Start/){
			#print "line number#:", $.,"\n";
			@properties = (split /\s*,\s*/,$_);
			#print $#properties,"\n";
			
			$currExe=$properties[1];
			
			$gridx	= $properties[2];
			$gridy	= $properties[3];
			$gridz	= $properties[4];
			$Throughput = $properties[12];
			#$nameIndex =  first_index{$_ =~ /Name/}@properties;
			$nameIndex =  first{$properties[$_] =~ /Name/}0..$#properties;


			
			#print "nameIndex: ", $nameIndex,"\n";
			#print "\n";
			
			#pop @properties if $properties[$#properties]=~/ID/;
			#map {print $_,","}@properties;
			
			#print "\n"; 
		
		}else{
			next if $_ =~/ns/;
			
			#print "line number##:", $.,"\n";
			#%hash_property = ();
			
			my @data = (split /\s*,\s*/,$_);
			#my $name = $data[$#data];
			#print "data: ",@data,"\n";
			#print $#data,"\n";
			#print $#properties,"\n";
			
			my %hash_property = map { $properties[$_] => $data[$_]}0..$nameIndex-1;
			map { $hash_property{$_} = 0 if $hash_property{$_} eq ""  } keys %hash_property;
			my $dataLastIndex = $#data - (($properties[$#properties]=~/ID/)? 1:0);
			my $name = join(",",@data[$nameIndex..$dataLastIndex]);
			
		
			#if($.==350){
			#	print $name,"\n";
			#	$name =~ s/\[([0-9]+)]//;
			#	print $name,"\n";
			#	#last;
			#}
			#$name = "" . $name;
			
			$name =~ s/\[([0-9]+)]//;
			$name =~ s/\s+//g;
			
			#$hash_property{$properties[$nameIndex]} = $name;
			#map {print $properties[$_],$data[$_]}0..$#data;
			#map {print $data[$_],","}0..$#data;
			
			$hash_property{ $callNum	}= 1;
			$hash_property{ $totalExe	}= $hash_property{$currExe};
			$hash_property{ $maxExe		}= $hash_property{$currExe};
			$hash_property{ $minExe		}= $hash_property{$currExe};
			$hash_property{ $avgExe		}= $hash_property{$currExe};
			
			$hash_property{ $maxGridx	}= $hash_property{$gridx};
			$hash_property{ $maxGridy	}= $hash_property{$gridy};
			$hash_property{ $maxGridz	}= $hash_property{$gridz};
			$hash_property{ $minGridx	}= $hash_property{$gridx};
			$hash_property{ $minGridy	}= $hash_property{$gridy};
			$hash_property{ $minGridz	}= $hash_property{$gridz};
			
			$hash_property{ $totalThroughput}= $hash_property{$Throughput};
			$hash_property{ $maxThroughput	}= $hash_property{$Throughput};
			$hash_property{ $minThroughput	}= $hash_property{$Throughput};
			$hash_property{ $avgThroughput	}= $hash_property{$Throughput};
			
			
			if(exists $hash_fh{$name}){
				
				$hash_property{$callNum}	= $hash_fh{$name}->{$callNum}+1;
				$hash_property{$totalExe}	= $hash_fh{$name}->{$totalExe}+$hash_property{$currExe};
				$hash_property{$avgExe}		= $hash_property{$totalExe}/$hash_property{$callNum};
				
				$hash_property{$maxExe}		= $hash_fh{$name}->{$maxExe} if $hash_fh{$name}->{$maxExe} > $hash_property{$maxExe};
				$hash_property{$minExe}		= $hash_fh{$name}->{$minExe} if $hash_fh{$name}->{$minExe} < $hash_property{$minExe};	
		
				
				$hash_property{$maxGridx}		= $hash_fh{$name}->{$maxGridx} if $hash_fh{$name}->{$maxGridx} > $hash_property{$maxGridx};
				$hash_property{$maxGridy}		= $hash_fh{$name}->{$maxGridx} if $hash_fh{$name}->{$maxGridy} > $hash_property{$maxGridy};
				$hash_property{$maxGridy}		= $hash_fh{$name}->{$maxGridy} if $hash_fh{$name}->{$maxGridz} > $hash_property{$maxGridz};
				
				$hash_property{$minGridx}		= $hash_fh{$name}->{$minGridx} if $hash_fh{$name}->{$minGridx} < $hash_property{$minGridx};	
				$hash_property{$minGridy}		= $hash_fh{$name}->{$minGridy} if $hash_fh{$name}->{$minGridy} < $hash_property{$minGridy};
				$hash_property{$minGridy}		= $hash_fh{$name}->{$minGridy} if $hash_fh{$name}->{$minGridy} < $hash_property{$minGridy};
		
				$hash_property{$maxThroughput}		= $hash_fh{$name}->{$maxThroughput} if $hash_fh{$name}->{$maxThroughput} > $hash_property{$maxThroughput};
				$hash_property{$minThroughput}		= $hash_fh{$name}->{$minThroughput} if $hash_fh{$name}->{$minThroughput} < $hash_property{$minThroughput};
				$hash_property{$totalThroughput}	= $hash_fh{$name}->{$totalThroughput}+$hash_property{$Throughput};
				$hash_property{$avgThroughput}		= $hash_property{$totalThroughput}/$hash_property{$callNum};
				
				
				#if ($.== 20||$.==21||$.==55){
				#	print $hash_property{$callNum},",";
				#}
			}
			
			$hash_fh{$name} = \%hash_property;
			
			#delete $hash_property{$delete_name};
			
			#if ($.== 20){
			#	#print $name;
			#	
			#	print map { "$_ => $hash_property{$_}\n" } keys %hash_property;
			#	#print $hash_iter{$name}->{$callNum},"\n";
			#	#print $hash_iter{$name},"\n";
			#	#print "\n";
			#}
			

			
			map { $hash_features{$_}++} keys %hash_property;
		
		}		
	}
	my @deleteIndex;
	my @deleteFeature = ("Start","Size", "Duration", "Grid X", "Grid Y", "Grid Z","Stream");
	foreach my $dp (@deleteFeature){
		#print "pb:", $dp,"\n";
		#push @deleteIndex, first_index{$_ =~ $dp}@properties;
		push @deleteIndex, first{$properties[$_] =~ $dp}0...$#properties;
	}
	
	#map {print "delete Index: ",$_,"\n"} @deleteIndex;
	map {delete $hash_features{$properties[$_]} } @deleteIndex;
		
	
	#map{print "sub_func: keys: ", $_,"\n"} keys %hash_fh ;
	return \%hash_fh;
					

}


### remove repeat name
#my %TPnamesHash;
#my @TPnamesDiff = grep{!$TPnamesHash{$_}++} @TPnames;
#my $stream = "gpu_streams";

my $max_x;
my $max_y;
my $max_z;

 
for my $sv(@servers){
	#my $svpath = "$currpath/$sv/$stream";
	my $svpath = "$currpath/$sv";
	chdir($svpath) or die "Can't chdir to $svpath $!";
	print "$svpath\n";
	my %hash_server;
	

#	if($sv eq "f3"){
#		$n_cu = 31;
#		$n_su = 1;
#	
#	}elsif($sv eq "ccsl"){
#		$n_cu = 7;
#		$n_su = 1;
#	}elsif($sv eq "debian"){
#		$n_cu = 11;
#		$n_su = 1;	
#	}elsif($sv eq "supermicro"){
#		$n_cu = 39;
#		$n_su = 1;	
#	}
	my ($max_x,$max_y,$max_z) = func_find_sz_max_xyz($sv);
	my @sz_zxy = func_calc_sz_zxy($max_x,$max_y,$max_z);
	
	for my $kernel (@kernels) {
		my %hash_kernel;
		my $prefix = ${nvprof}.'_'.${kernel};
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
						
						my $filename = $prefix.'_'.${i}.'_'.${j}.'_'.${k}.'_'.$m.'_'.${n_reps}.'_'.${IsStream}.'_'.${nt}.'.txt';
						
						unless (-e $filename) {
							print "$filename does not exist.\n"; 
							next;
						}
						open my $fh, '<', $filename or die "Cannot open $filename: $!";
						print "Processing $filename...\n";
						
						$hash_iter{$nt} = func_process_file($fh);
						#map{print "hash_iter{0} keys: ", $_,"\n"} keys %{$hash_iter{$nt}}; ;
						map { $hash_metrics{$_}++} keys %{$hash_iter{$nt}};						
						
					}elsif($IsStream==1){
						$ntVal = min(ceil(min($i,$j,$k)/$basicBlock),$ntMax);
						
                        for(my $nt=1;$nt<$ntVal+1;++$nt){
							
							my $filename = $prefix.'_'.${i}.'_'.${j}.'_'.${k}.'_'.$m.'_'.${n_reps}.'_'.${IsStream}.'_'.${nt}.'.txt';
							
							unless (-e $filename) {
							print "$filename does not exist.\n"; 
							next;
							}
							open my $fh, '<', $filename or die "Cannot open $filename: $!";
							print "Processing $filename...\n";	
							
							$hash_iter{$nt} = func_process_file($fh);
							map { $hash_metrics{$_}++} keys %{$hash_iter{$nt}};
							
						}
						
					}					

					if(%hash_iter){	
						$hash_IsStream{$m}=\%hash_iter;
					}
				}
				
				#map{print "keys: ", $_,"\n"} keys %hash_metrics ;
				$hash_sz{$IsStream}=\%hash_IsStream;
			}
			if(%hash_sz){
				$hash_kernel{$xyz} = \%hash_sz;
			}		
		}
		$hash_server{$kernel}=\%hash_kernel;
	}
	$hash_table{$sv}=\%hash_server;
}


#map{print "features: ", $_,"\n"} keys %hash_features ;

###################print output##############################
my $postpose = "nvprof_postpose";
my @output_server;

my @out_features;
my @features = keys %hash_features;
push @out_features,"machine","kernel","size(x)","size(y)","size(z)","iteration","IsStream","nTile","metric";
push @out_features,@features;
push @output_server,[@out_features];

#map{print "output features: ", $_,"\n"} @out_features ;


for my $sv(@servers){
	my ($max_x,$max_y,$max_z) = func_find_sz_max_xyz($sv);
	my @sz_zxy = func_calc_sz_zxy($max_x,$max_y,$max_z);
	for my $kernel (@kernels) {
		
		#my @out_kel;
		
	
		#my $postfix = $postpose.'_'.$kernel;
		for my $sz(@sz_zxy){
		
			my ($k,$i,$j) = @$sz;
			my $xyz=$i.'_'.$j.'_'.$k;
		
			next if not exists $hash_table{$sv}->{$kernel}->{$xyz};
		
			for(my $IsStream=0;$IsStream<2;++$IsStream){		
				for(my $m=1;$m<$n_iter;++$m){
				
					if($IsStream==0){
						my $nt = 0;
						for my $metric  (keys %hash_metrics){
							my @out_val;
							push @out_val,$sv,$kernel,$i,$j,$k,$m,$IsStream,$nt,$metric;
							#print "metric:", $metric,"\n";
							for my $fet (@features){
								#print "features:", $fet,"\n";
								my $val;
								
															
								if (exists $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{$nt}->{$metric}->{$fet} ){
									$val = $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{$nt}->{$metric}->{$fet} ;
								}else{
									$val = 0;
								}
								push @out_val, $val;
							}
							push @output_server,[@out_val];
							
						}	
					
					}elsif ($IsStream==1){
						$ntVal = min(ceil(min($i,$j,$k)/$basicBlock),$ntMax);
                        for(my $nt=1;$nt<$ntVal+1;++$nt){

							for my $metric  (keys %hash_metrics){
								my @out_val;
								push @out_val,$sv,$kernel,$i,$j,$k,$m,$IsStream,$nt,$metric;
								for my $fet (@features){
									my $val;
							
									if (exists $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{$nt}->{$metric}->{$fet} ){
										$val = $hash_table{$sv}->{$kernel}->{$xyz}->{$IsStream}->{$m}->{$nt}->{$metric}->{$fet};
									}else{
										$val = 0;
									}
									push @out_val, $val;
								}
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


#my @threads = (
#    [31,1],
##	[15,2]#,
##   [11,4],
##   [05,8]
#);
#
#
#my $n_iter = 11;
#my $n_reps = 5; 
#my $its = 1;
#
#
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
##my $filename = shift @ARGV;
##$filename = "output" unless ($filename);
#
##   Possible kernels:
##   OldLaplaceSolver NewLaplaceSolver SeqLaplaceSolver CompleteLaplaceSolver 
##   NaivePtrLaplaceSolver OMPLaplaceSolver ReallocLaplaceSolver 
##   LaplaceSolverCDLoop NaiveLaplaceSolver
#
##my @kernels = qw/Seq OMP OMPFG OMPFGV2 FineGrain FineGrainM2V2 FineGrainM2 InPlace Naive NaiveTPsPtr InPlaceTPs FineGrainTPs/; # Seq/; 
#my @kernels = qw/StencilCudaGpu37 /; # Seq/; 
##my @kernels = qw/StencilCudaHybrid3 /; # Seq/; 
#my @outputs = ();
#my $sz_start = 1000;
#my $sz_end = 51000;
#my $sz_step = 2000;
#
#my $nvprof = "nvprof_trace";
##for (my $i= $sz_start;$i<$sz_end;$i=$i+$sz_step){
#for my $sz(@sz_zxy){
#    my ($k,$i,$j) = @$sz;
#    for my $thd (@threads) {
#        my ($n_cu,$n_su) = @$thd;
#        $ENV{'DARTS_NUM_CU'}    = $n_cu;
#        $ENV{'DARTS_NUM_SU'}    = $n_su;
#        $ENV{'OMP_PROC_BIND'}   = 'true';
#        $ENV{'OMP_NUM_THREADS'} = ($n_cu+1)*$n_su;
#		my $n_total = ($n_cu+1)*$n_su;
#        #if ($n_total <= 16){
#        #    #$ENV{'GOMP_CPU_AFFINITY'}="0-7 16-23";
#		#	$ENV{'DARTS_AFFINITY'}="0-7 16-23";
#        #}else{
#        #    #$ENV{'GOMP_CPU_AFFINITY'}="0-31";
#		#	$ENV{'DARTS_AFFINITY'}="0-31";
#        #}
#		$ENV{'DARTS_AFFINITY'}="0-${n_cu}";
#        for my $kernel (@kernels) {
#            my $ker = './' . $kernel;
#            if (-e $ker) {
#				for(my $m=1;$m<$n_iter;++$m){
#					open my $fh, '>>', "${kernel}_${i}_${j}_${k}_${n_cu}_${n_su}_${n_reps}_${m}.txt" or die $!;
#					for (my $idx = 0; $idx < $its; ++$idx) {
#						print $fh  `numactl --interleave=0-1 $ker $i $j $k $m $n_reps `;
#						system `nvprof --print-gpu-trace --print-api-trace --csv -u ns $ker $i $j $k $m $n_reps 2> ${nvprof}_${kernel}_${i}_${j}_${k}_${n_cu}_${n_su}_${n_reps}_${m}.txt`; 
#	
#					}
#				}
#            } else {
#                print "Couldn't find kernel named $ker\n";
#            }
#        }
#    }
#}



