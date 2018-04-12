#!/usr/bin/perl

use warnings;
use strict;

use Cwd;

my $currpath = Cwd::abs_path();
print "$currpath\n";

#hw_information
my @hwinfo = qw(Server_Name NumOfSocket CPU_Type CPU_Clock(GHz) CPU_Cores CPU_Threads  CPU_Mem(GB) CPU_L3_Cache(MB) GPU_Type GPU_Clock(GHz)  NumOfSM SP/SM  GPU_Mem(GB) GPU_L2_Cache(MB) CC NumOfCE PCIe);
#print @hwinfo;

my %all_hw_info;

my @servers = qw/fatnode supermicro debian hive ccsl/;

for my $server (@servers){
	my %server_info;
	$server_info{"Server_Name"} = $server;
	if ($server eq "fatnode"){
		$server_info{"NumOfSocket"} = 2 ;
		$server_info{"CPU_Type"} = "Intel";
		$server_info{"CPU_Clock(GHz)"} = 2.6;		
		$server_info{"CPU_Cores"} = 16;
		$server_info{"CPU_Threads"} = 32;
		$server_info{"CPU_Mem(GB)"} = "64(2x32)";
		$server_info{"CPU_L3_Cache(MB)"} = 20;
		$server_info{"GPU_Type"} = "K20";
		$server_info{"GPU_Clock(GHz)"} = 0.71;
		$server_info{"NumOfSM"} =13 ;
		$server_info{"SP/SM"} = 192;
		$server_info{"GPU_Mem(GB)"} = 4.8;
		$server_info{"GPU_L2_Cache(MB)"} = 1.25;		
		$server_info{"CC"} = 3.5;
		$server_info{"NumOfCE"} = 2;	
		$server_info{"PCIe"} = 6.2;		
	
	}elsif ($server eq "supermicro"){
		$server_info{"NumOfSocket"} = 2;
		$server_info{"CPU_Type"} = "Intel";
		$server_info{"CPU_Clock(GHz)"} = 3;		
		$server_info{"CPU_Cores"} = 20;
		$server_info{"CPU_Threads"} = 40;
		$server_info{"CPU_Mem(GB)"} = "256(2x128)";
		$server_info{"CPU_L3_Cache(MB)"} = 25;
		$server_info{"GPU_Type"} = "K20";
		$server_info{"GPU_Clock(GHz)"} = 0.71;
		$server_info{"NumOfSM"} =13 ;
		$server_info{"SP/SM"} = 192;
		$server_info{"GPU_Mem(GB)"} = 4.8;
		$server_info{"GPU_L2_Cache(MB)"} = 1.25 ;		
		$server_info{"CC"} = 3.5;
		$server_info{"NumOfCE"} = 2;	
		$server_info{"PCIe"} = 6.1;	
	
	}elsif ($server eq "debian"){
		$server_info{"NumOfSocket"} = 1;
		$server_info{"CPU_Type"} = "Intel";
		$server_info{"CPU_Clock(GHz)"} = 3.4;		
		$server_info{"CPU_Cores"} = 6;
		$server_info{"CPU_Threads"} = 12;
		$server_info{"CPU_Mem(GB)"} = "31";
		$server_info{"CPU_L3_Cache(MB)"} = 12;
		$server_info{"GPU_Type"} = "Titan" ;
		$server_info{"GPU_Clock(GHz)"} = 0.88;
		$server_info{"NumOfSM"} = 14;
		$server_info{"SP/SM"} = 192;
		$server_info{"GPU_Mem(GB)"} = 6;
		$server_info{"GPU_L2_Cache(MB)"} = 1.5;		
		$server_info{"CC"} = 3.5;
		$server_info{"NumOfCE"} = 2;	
		$server_info{"PCIe"} = 11.5;	
	
	}elsif ($server eq "hive"){
		$server_info{"NumOfSocket"} = 1;
		$server_info{"CPU_Type"} = "Intel";
		$server_info{"CPU_Clock(GHz)"} = 3.2;		
		$server_info{"CPU_Cores"} = 6;
		$server_info{"CPU_Threads"} = 12;
		$server_info{"CPU_Mem(GB)"} = "30";
		$server_info{"CPU_L3_Cache(MB)"} = 12;
		$server_info{"GPU_Type"} = "GTX680";
		$server_info{"GPU_Clock(GHz)"} = 1.06;
		$server_info{"NumOfSM"} = 8;
		$server_info{"SP/SM"} = 192;
		$server_info{"GPU_Mem(GB)"} = 2;
		$server_info{"GPU_L2_Cache(MB)"} = 0.5;		
		$server_info{"CC"} = 3;
		$server_info{"NumOfCE"} = 1;	
		$server_info{"PCIe"} = 5.6;	
	
	}elsif ($server eq "ccsl"){	
		$server_info{"NumOfSocket"} = 1;
		$server_info{"CPU_Type"} = "Intel";
		$server_info{"CPU_Clock(GHz)"} = 3.4;		
		$server_info{"CPU_Cores"} = 4;
		$server_info{"CPU_Threads"} = 8;
		$server_info{"CPU_Mem(GB)"} = "16";
		$server_info{"CPU_L3_Cache(MB)"} = 8;
		$server_info{"GPU_Type"} = "K20";
		$server_info{"GPU_Clock(GHz)"} = 0.75;
		$server_info{"NumOfSM"} = 15;
		$server_info{"SP/SM"} = 192;
		$server_info{"GPU_Mem(GB)"} = 12;
		$server_info{"GPU_L2_Cache(MB)"} =1.5 ;		
		$server_info{"CC"} = 3.5 ;
		$server_info{"NumOfCE"} = 2;	
		$server_info{"PCIe"} = 10.3;		
	}
	$all_hw_info{$server}=\%server_info;
}

my @kernels = qw/StencilCudaHybrid4/;
my @gbase = (1000,2000,4000,8000);
my @ratio = (0.5,1,1.5,2);
my @msize = (17000,19000,21000,23000,25000);


my $strongOrweak = 'weak';
my $platform = 'darts';
my $type= 'hybrid';
my $sufix = 'p';


my %all_exe_time;
my $n_su = 1;


for my $server (@servers){
	my $sv; 
	if ($server eq "fatnode"){
		$sv = "f4";
	}else {
		$sv = $server;
	}
	
	my $fdname = $sv.'_'.$strongOrweak.'_'.$platform.'_'.$type.'_'.$sufix;
	my $path = "$currpath/$fdname";
	chdir($path) or die "Can't chdir to $path $!";	
	print "$path\n";
	
	my $thread;

	if ($sv eq "f4"){
		$thread = [31,1];
	}elsif($sv eq "supermicro"){
		$thread = [39,1];
	}elsif($sv eq "debian"){
		$thread = [11,1];
	}elsif($sv eq "hive"){
		$thread = [11,1];
	}elsif($sv eq "ccsl"){
		$thread = [7,1];
	}
	my $n_cu = $thread->[0];

	my $total_cores=($n_cu+1)*$n_su;
	my $tc=$total_cores;
	my $tc1=$tc-1;
	my $n_rept=20;
	my $startPos=8;
	my $endPos = $n_rept+$startPos-1;
	my %exe_time;
	for my $k (@kernels) {
		my %k_time;
		for my $msz (@msize){
			my %msz_time;
			for my $gb (@gbase){
				my %gb_time;
				for my $ra (@ratio){
								
					my $filename = $k . '_' . $msz . '_' . $n_cu . '_' . $n_su .'_'.$gb.'_'.$ra. '.txt';
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
					$gb_time{$ra} = sprintf ("%.2f",$avg);	
				}
				$msz_time{$gb} = \%gb_time;	
			}
			$k_time{$msz} = \%msz_time;
		}
		$exe_time{$k} = \%k_time;
	}
	$all_exe_time{$server}=\%exe_time;
}

chdir($currpath) or die "Can't chdir to $currpath $!";
print "$currpath\n";
#output files

my @exe_info=();
push @exe_info, [@hwinfo,'Total_Workload','Initial_GPU_Workload','CPU/GPU_Workload_Ratio','Exe_Time'];


for my $server (@servers){

	my $thread;

	if ($server eq "fatnode"){
		$thread = [31,1];
	}elsif($server eq "supermicro"){
		$thread = [39,1];
	}elsif($server eq "debian"){
		$thread = [11,1];
	}elsif($server eq "hive"){
		$thread = [11,1];
	}elsif($server eq "ccsl"){
		$thread = [7,1];
	}
	my $n_cu = $thread->[0];

#	my @v_hw_info =(  $all_hw_info{$server}->{"Server_Name"}
#					  ,$all_hw_info{$server}->{"NumOfSocket"}
#					  ,$all_hw_info{$server}->{"NumOfSocket"}			
#					  ,$all_hw_info{$server}->{"CPU_Type"}			
#					  ,$all_hw_info{$server}->{"CPU_Clock(GHz)"} 				
#					  ,$all_hw_info{$server}->{"CPU_Cores"} 			
#					  ,$all_hw_info{$server}->{"CPU_Threads"} 		
#					  ,$all_hw_info{$server}->{"CPU_Mem(GB)"} 			
#					  ,$all_hw_info{$server}->{"CPU_L3_Cache(MB)"} 			
#					  ,$all_hw_info{$server}->{"GPU_Type"} 			
#					  ,$all_hw_info{$server}->{"GPU_Clock(GHz)"}			
#					  ,$all_hw_info{$server}->{"NumOfSM"} 			
#					  ,$all_hw_info{$server}->{"SP/SM"} 			
#					  ,$all_hw_info{$server}->{"GPU_Mem(GB)"} 			
#					  ,$all_hw_info{$server}->{"GPU_L2_Cache(MB)"} 			
#					  ,$all_hw_info{$server}->{"CC"}			
#					  ,$all_hw_info{$server}->{"NumOfCE"}
#					  ,$all_hw_info{$server}->{"PCIe"}		
#					);
	
	my @v_hw_info;
	for my $hw (@hwinfo){
		push @v_hw_info,$all_hw_info{$server}->{$hw };
	
	}

	for my $k (@kernels) {
	
		for my $msz (@msize){
			for my $gb (@gbase){					
				for my $ra (@ratio){		
					push @exe_info, [@v_hw_info,$msz,$gb,$ra,$all_exe_time{$server}->{$k}->{$msz}->{$gb}->{$ra}];	
				}
			}
		}
	}	
}

my $file_out = 'weak_execution'.'.csv';
open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
print $fh join(',',@$_) . "\n" for (@exe_info);




###
###
###
###
###my @kernels = qw/StencilCudaHybrid4/;
###my @gbase = (1000,2000,4000,8000);
###my @ratio = (0.5,1,1.5,2);
###my @msize = (17000,19000,21000,23000,25000);
###
###my @threads=(
###	[31,1],
###	[39,1],
###	[11,1],
###	[11,1],
###	[8,1],
###	);
###	
###
###	
###	
###my %exe_time;
###my $n_su = 1;
###my $n_cu;
###
###for my $thread (@threads) {	
###	$n_cu = $thread->[0];
###	
###	my $total_cores=($n_cu+1)*$n_su;
###	my $tc=$total_cores;
###	my $tc1=$tc-1;
###	my $n_rept=20;
###	my $startPos=8;
###	my $endPos = $n_rept+$startPos-1;
###	
###	
###	for my $k (@kernels) {
###		my %k_time;
###		for my $msz (@msize){
###			my %msz_time;
###			for my $gb (@gbase){
###				my %gb_time;
###				for my $ra (@ratio){
###								
###					my $filename = $k . '_' . $msz . '_' . $n_cu . '_' . $n_su .'_'.$gb.'_'.$ra. '.txt';
###					unless (-e $filename) {
###						print "$filename does not exist.\n"; 
###						next;
###					}
###					open my $fh, '<', $filename or die "Cannot open $filename: $!";
###					print "Processing $filename...\n";				
###					my @data_t;
###					while (<$fh>) {
###						chomp;
###						my @data = (split /\s*,\s*/,$_)[$startPos..$endPos];
###						my @data_s = sort{$a <=> $b}@data;
###						shift @data_s;
###						pop @data_s;
###						@data_t=(@data_t,@data_s);
###					}
###					my $avg=0;
###					$avg +=$_ for @data_t;
###					$avg /=@data_t;
###					$gb_time{$ra} = sprintf ("%.2f",$avg);	
###				}
###				$msz_time{$gb} = \%gb_time;	
###			}
###			$k_time{$msz} = \%msz_time;
###		}
###		$exe_time{$k} = \%k_time;
###	}
###}
###
####output files
###
###for my $thread (@threads) {	
###	$n_cu = $thread->[0];
###	
###	my $total_cores=($n_cu+1)*$n_su;
###	if($total_cores >$tc){
###		next;
###	}
###	for my $k (@kernels) {
###		
###		for my $msz (@msize){
###			my @v_exe_time=();
###			push @v_exe_time,['#ratio', @gbase];			
###			for my $ra (@ratio){
###				my @line;
###				push @line, $ra,$total_cores;			
###				for my $gb (@gbase){
###					
###					push @line, $exe_time{$k}->{$msz}->{$gb}->{$ra};	
###				}
###				push @v_exe_time, [@line];
###				
###			}
###			my $file_out = $k.'_'.$n_su.'_'.$msz.'_'.'weak_execution'.'.dat';
###			open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
###			print $fh join(',',@$_) . "\n" for (@v_exe_time);
###		}
###		
###		
###	}
###
###}



