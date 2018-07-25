#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;

use Cwd;

#use Cwd qw(abs_path);
my $currpath = Cwd::abs_path();
print "$currpath\n";

my @currpath_parts = split('/', $currpath);
my $grandparent = join "/",map { $_ } @currpath_parts[0..@currpath_parts-3];
print "$grandparent\n";
my $dataFolderName = "data";
my $datapath = "$grandparent/$dataFolderName";
chdir($datapath) or die "Can't chdir to $datapath $!";	
print "$datapath \n";

my @test =qw/gpu_streams gpu_streamNot/;

my @servers = qw/f3 supermicro debian ccsl/;

my @kernels = qw/StencilCudaGpu37 /;
my $k1 = "StencilCudaGpu37";

my @sz_zxy = (
    [50  ,100 ,100 ],
    [50  ,200 ,200 ],
    [50  ,400 ,400 ],
    [50  ,800 ,800 ],
    [50  ,1000,1000],
    [100 ,100 ,100 ],
    [100 ,200 ,200 ],
    [100 ,400 ,400 ],
    [100 ,800 ,800 ],
    [100 ,1000,1000],
    [200 ,100 ,100 ],
    [200 ,200 ,200 ],
    [200 ,400 ,400 ],
    [200 ,800 ,800 ],
    [200 ,1000,1000],	
    [400 ,100 ,100 ],
    [400 ,200 ,200 ],
    [400 ,400 ,400 ],
    [400 ,800 ,800 ],
    [400 ,1000,1000],
    [800 ,100 ,100 ],
    [800 ,200 ,200 ],
    [800 ,400 ,400 ],
    [800 ,800 ,800 ],
    [800 ,1000,1000],	
    [1000,100 ,100 ],
    [1000,200 ,200 ],
    [1000,400 ,400 ],
    [1000,800 ,800 ],
    [1000,1000,1000]

);
my $iter_start = 1;
my $iter_end   = 11;
my $iter_step  = 1;
my $n_rept = 5;
my $startPos=8;
my $endPos = $n_rept+$startPos-1;

my %exe_time_gpu;

for my $server (@servers){

	my %exe_server;
	my $thread;
	my $xEnd;
	my $yEnd;
	
	if ($server eq "f3"){
		$thread = 31;
	}elsif($server eq "supermicro"){
		$thread = 39;
	}elsif($server eq "debian"){
		$thread = 11;
	}elsif($server eq "ccsl"){
		$thread = 7;
	}
	my $nCores = $thread+1;
	my $n_cu = $thread;
	my $n_su = 1;
	for my $tp (@test){
		my %exe_tp;
		my $testgpupath = "$datapath/$server/$tp";
		chdir($testgpupath) or die "Can't chdir to $testgpupath $!";
		print "$testgpupath \n";
		for my $sz (@sz_zxy){
		    my ($z,$x,$y) = @$sz;
			my $ratio = 0;
			my $value = $z.'*'.$x.'*'.$y;
			my $value1 = $x.'_'.$y.'_'.$z;
			
			my $k = $k1;
			my %exe_sz;
			for (my $it=$iter_start; $it<$iter_end;$it=$it+$iter_step){
			
				my $filename = $k . '_' . $value1 . '_' . $n_cu . '_' . $n_su .'_'.$n_rept.'_'.$it.'.txt';
				unless (-e $filename) {
					print "$filename does not exist.\n";
					$exe_sz{$it}=0;
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
				if (@data_t == 0){
					$avg = 0;
				}else{
					$avg /=@data_t;
				}
				$exe_sz{$it}=sprintf ("%.2f",$avg);							
			}
			$exe_tp{$sz}=\%exe_sz;
		}	
		$exe_server{$tp}=\%exe_tp;
	}
	$exe_time_gpu{$server}=\%exe_server;
}


chdir($currpath) or die "Can't chdir to $currpath $!";	
#print results to files (weak file)
my $etype = "weak_execution";

for my $server (@servers){
	for my $tp (@test){
		my @v_exe_time;
		push @v_exe_time, ["#size","volumn" ,($iter_start..$iter_end-1)];
		for my $sz (@sz_zxy){
			my ($z,$x,$y) = @$sz;
			
			my $value = $z.'*'.$x.'*'.$y;
			my $volumn = $z*$x*$y;
			my @line;
			push @line, $value, $volumn;
			for (my $it=$iter_start; $it<$iter_end;$it=$it+$iter_step){
				my $v = $exe_time_gpu{$server}->{$tp}->{$sz}->{$it};
				if ($v==0){				
					next;
				};
				push @line,$v;
			}
			if(@line !=2){
				push @v_exe_time, [@line];
			}
		}
		my $file_out = $etype.'_'.$server.'_'.$tp.'.csv';
		open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
		print $fh join(',',@$_) . "\n" for (@v_exe_time);
	}
}

my $etype1 = "S-NS";
for my $server (@servers){
	my @v_exe_time;
	push @v_exe_time, ["#size","volumn" ,($iter_start..$iter_end-1)];
	for my $sz (@sz_zxy){
		my ($z,$x,$y) = @$sz;
		
		my $value = $z.'*'.$x.'*'.$y;
		my $volumn = $z*$x*$y;
		my @line;
		push @line, $value, $volumn;
		for (my $it=$iter_start; $it<$iter_end;$it=$it+$iter_step){
			my @tmp=();
			for my $tp (@test){
				my $v = $exe_time_gpu{$server}->{$tp}->{$sz}->{$it};
				push @tmp, $v;
			}
			if(!(0~~@tmp)){
				my $diff = $tmp[0]-$tmp[1];
				sprintf ("%.2f",$diff);
				push @line, $diff;
			}
			
		}
		if(@line !=2){
			push @v_exe_time, [@line];
		}
	}
	
	my $file_out = $etype.'_'.$etype1.'_'.$server.'.csv';
	open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
	print $fh join(',',@$_) . "\n" for (@v_exe_time);

}



#my $sw = 'weak';
#my $GC= 'DAWL_2GB';
#my $thd='threads';
#my $suf='compact';
#
##for my $server (@servers){
##	my $fdname = $server.'_'.$sw.'_'.$GC;
##	my $path = "$currpath/$fdname";
##	chdir($path) or die "Can't chdir to $path $!";
##	print "$path\n";
##	chdir($currpath) or die "Can't chdir to $currpath $!";
##}
#
##my @vals = (1 .. 10);
###my @vals = (1,3,5);
#
##my @kernels = qw/StencilCudaCpu StencilCudaGpu StencilCudaHybrid3/;
#my @kernels = qw/StencilSEQ37 StencilCUDA37 StencilCudaCpu37 StencilCudaGpu37 StencilCudaHybrid37/;
#
#
#
#
##my @stats = ([@kernels]);
#
#my @threads=(
#	[31,1],
#	[39,1],
#	[11,1],
#	[11,1],
#	[8,1],
##	[15,2]
#	);
#
#my @sz_zxy = (
#    [50,200,200],
#    [100,200,200],
#    [200,200,200],
#    [800,400,400],
#    [200,800,800],
#    [400,800,800],
#    [800,800,800],
#    [800,1000,1000],
#    [1000,1000,1000]
#);
#	
#	
#my $tc=40;
#my $tc1=$tc-1;
#my $n_rept=10;
#my $startPos=8;
#my $endPos = $n_rept+$startPos-1;
#my $sz_start = 1000;
#my $sz_end = 50000;	
#my $sz_step = 2000; 
#
#
#for my $server (@servers){
#
#	my $thread;
#	my $xEnd;
#	my $yEnd;
#	#for my $thread (@threads) {
#	if ($server eq "f4"){
#		$thread = 31;
#		$sz_end = 50000;
#		$xEnd = 1000;
#		$yEnd = 1000;
#	}elsif($server eq "supermicro"){
#		$thread = 39;
#		$sz_end = 50000;
#		$xEnd = 1000;
#		$yEnd = 1000;
#	}elsif($server eq "debian"){
#		$thread = 11;
#		$sz_end = 37000;
#		$xEnd=1000;
#		$yEnd=1000;
#	}elsif($server eq "hive"){
#		$thread = 11;
#		$sz_end = 37000;
#	}elsif($server eq "ccsl"){
#		$thread = 7;
#		$sz_end = 26000;
#		$xEnd=800;
#		$yEnd=800;
#	}
#	my $nCores = $thread+1;
#	#my $fdname = $server.'_'.$sw.'_'.$GC.'_'.$nCores.$thd.'_'.$suf;
#	my $fdname = $server.'_'.$sw.'_'.$GC;
#	my $path = "$currpath/$fdname";
#	chdir($path) or die "Can't chdir to $path $!";	
#	print "$path\n";
#	
#	my @t_exe_time;
#	my @t_speedup;
#	push @t_exe_time, ["#size",@kernels];
#	push @t_speedup, ["#size",@kernels];
#	
#	#my ($n_cu,$n_su) = @$thread;
#	#for (my $v=$sz_start; $v<$sz_end;$v=$v+$sz_step){
#	#	my $value = $v;
#	for my $sz (@sz_zxy){
#	    my ($z,$x,$y) = @$sz;
#		my $ratio = 0;
#		if(($x>$xEnd)or ($y>$yEnd)){
#			next;
#		}
#		my $value = $z.'*'.$x.'*'.$y;
#		my $value1 = $x.'_'.$y.'_'.$z;
#		my @line  = ($value);
#		my @line_s  = ($value);
#		for my $k (@kernels) {
#			my @files = glob("$k\_$value1*");
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
#
#			
#		}
#		push @t_exe_time, [@line];
#		push @t_speedup, [@line_s];
#		
#	}
#	chdir($currpath) or die "Can't chdir to $currpath $!";	
#	my $file_out = $server.'_'."weak_execution.dat";
#	open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
#	print $fh join(',',@$_) . "\n" for (@t_exe_time);
#	
#	my $file_out_s = $server.'_'.'weak_speedup'.'.dat';
#	open my $fh_s, '>', $file_out_s or die "couldn't open $file_out_s: $!";
#	print $fh_s join(',',@$_) . "\n" for (@t_speedup);
#		
#	#}
#	
#}






#for my $server (@servers){
#	my $fdname = $server.'_'.$sw.'_'.$GC;
#	my $path = "$currpath/$fdname";
#	chdir($path) or die "Can't chdir to $path $!";	
#	print "$path\n";
#	my $thread;
#	#for my $thread (@threads) {
#		if ($server eq "f4"){
#			$thread = [31,1];
#			$sz_end = 50000;
#		}elsif($server eq "supermicro"){
#			$thread = [39,1];
#			$sz_end = 50000;
#		}elsif($server eq "debian"){
#			$thread = [11,1];
#			$sz_end = 37000;
#		}elsif($server eq "hive"){
#			$thread = [11,1];
#			$sz_end = 37000;
#		}elsif($server eq "ccsl"){
#			$thread = [7,1];
#			$sz_end = 26000;
#		}
#		my @t_exe_time;
#		my @t_speedup;
#		push @t_exe_time, ["#size",@kernels];
#		push @t_speedup, ["#size",@kernels];
#		
#		my ($n_cu,$n_su) = @$thread;
#		for (my $v=$sz_start; $v<$sz_end;$v=$v+$sz_step){
#			my $value = $v;
#			my $ratio = 0;
#			my @line  = ($value);
#			my @line_s  = ($value);
#			for my $k (@kernels) {
#				my $filename = $k . '_' . $value . '_' . $n_cu . '_' . $n_su . '.txt';
#				unless (-e $filename) {
#					print "$filename does not exist.\n"; 
#					next;
#				}
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
#				
#			}
#			push @t_exe_time, [@line];
#			push @t_speedup, [@line_s];
#			
#		}
#		chdir($currpath) or die "Can't chdir to $currpath $!";	
#		my $file_out = $server.'_'."weak_execution.dat";
#		open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
#		print $fh join(',',@$_) . "\n" for (@t_exe_time);
#		
#		my $file_out_s = $server.'_'.'weak_speedup'.'.dat';
#		open my $fh_s, '>', $file_out_s or die "couldn't open $file_out_s: $!";
#		print $fh_s join(',',@$_) . "\n" for (@t_speedup);
#		
#	#}
#	
#}
	

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


