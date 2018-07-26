#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;


use Cwd;

#use Cwd qw(abs_path);
my $currpath = Cwd::abs_path();
print "$currpath\n";


my @servers = qw/f3 supermicro debian ccsl/;
my @gputest = qw/gpu_streams gpu_streamNot/;

my $prefix = "nvprofTraces";



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


my %h_total;
my %h_items;
my %h_test;





for my $tp (@gputest){

	for my $sv(@servers){

		my $filename = $prefix . '_' . $sv . '_' . $tp .'.csv';
		unless (-e $filename) {
			print "$filename does not exist.\n";
			next;				
		}
		open my $fh, '<', $filename or die "Cannot open $filename: $!";
		print "Processing $filename...\n";
		
		my @items =(split /\s*,\s*/,<$fh>);
		chomp (@items);
		#print @items;
		my @data;
		while (<$fh>) {
			chomp;
			my @data1= (split /\s*,\s*/,$_);
			push @data, [@data1];	
			my @tt=@data1[4..8];
			$h_total{$tp}{$data1[0]}{$data1[1]}{$data1[2]}{$data1[3]}=\@tt;
		}
		my %seen=();
		for(my $i=0;$i<@items-5;++$i){
			my @tt = grep { ! $seen{$_} ++ } map {$_->[$i]} @data;
			$h_items{$items[$i]}=\@tt;
			#print $items[$i],$h_items{$items[$i]}->[0],"\n";
		}		
	}
}


#print $h_total{"gpu_streams"}->{"f3"}->{"100_100_50"}->{"1"}->{"total_running_time"}->[-1];

###print output


my @checks = ("total_CUDA memcpy HtoD", "total_CUDA memcpy DtoH" ,"total_kernel_execution_time" ,"total_CUDA_functions" ,"total_running_time" );
my @outputNames = ("transHtoD","transDtoH","kernel","cudadriver","gpu_running");

my $out_suf0  = "overlap";

my @out0;

my $postpose = "nvprof_postpose";
if(-d $postpose ){
	unlink glob "./$postpose/*";
	rmdir $postpose;
}
unless(mkdir $postpose) {
    die "Unable to create $postpose\n";
}
my $outpath = "$currpath/$postpose";;

my @tmp =("stream Y/NO", "sever","size","volumn");
my @its = map{"its=".$_} (1..10);
push @out0,[@tmp,@its];

chdir($outpath) or die "Can't chdir to $outpath $!";

for (my $ck =0; $ck<@outputNames;++$ck){
	
	my @out_ck;
	push @out_ck,[@tmp,@its];

	for my $tp (@gputest){
		
		for my $sv(@servers){
			for my $sz (@sz_zxy){
				my ($z,$x,$y) = @$sz;
	
				my $value = $z.'*'.$x.'*'.$y;
				my $value1 = $x.'_'.$y.'_'.$z;
				my $volumn = $z*$x*$y;
				my @line;
			
				push @line,$tp, $sv,$value,$volumn;
			
				
				for (my $it=1;$it<11;++$it){
					
					if (((exists $h_total{$tp}{$sv}{$value1}{$it}))&&($h_total{$tp}{$sv}{$value1}{$it}{$checks[4]}[-1]!=0)){
						my $vl;
						if ($ck < (@outputNames-1)){
							$vl = ($h_total{$tp}{$sv}{$value1}{$it}{$checks[0]}[-1]);
						}else{
							$vl = ($h_total{$tp}{$sv}{$value1}{$it}{$checks[0]}[-1])
							+	($h_total{$tp}{$sv}{$value1}{$it}{$checks[1]}[-1])
							+	($h_total{$tp}{$sv}{$value1}{$it}{$checks[2]}[-1])
							+	($h_total{$tp}{$sv}{$value1}{$it}{$checks[3]}[-1])
							-	($h_total{$tp}{$sv}{$value1}{$it}{$checks[4]}[-1]);
						}
											
		
						
						push @line,$vl;
						
					}
	
				}
				if(@line>5){
					push @out_ck,[@line];	
				}
				
			}
		
		}
	
	}
	
	my $file_out = $prefix.'_'.$outputNames[$ck].'.csv';
	open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
	print $fh join(',',@$_) . "\n" for (@out_ck);

}

chdir($currpath) or die "Can't chdir to $outpath $!";

#for my $tp (@gputest){
#	
#	for my $sv(@servers){
#		for my $sz (@sz_zxy){
#		    my ($z,$x,$y) = @$sz;
#
#			my $value = $z.'*'.$x.'*'.$y;
#			my $value1 = $x.'_'.$y.'_'.$z;
#			my $volumn = $z*$x*$y;
#			my @line0;
#		
#			push @line0,$tp, $sv,$value,$volumn;
#		
#		
#			for (my $it=1;$it<11;++$it){
#				if (((exists $h_total{$tp}{$sv}{$value1}{$it}))&&($h_total{$tp}{$sv}{$value1}{$it}{$checks[4]}[-1]!=0)){
#
#					my $overlap = ($h_total{$tp}{$sv}{$value1}{$it}{$checks[0]}[-1])
#					+	($h_total{$tp}{$sv}{$value1}{$it}{$checks[1]}[-1])
#					+	($h_total{$tp}{$sv}{$value1}{$it}{$checks[2]}[-1])
#					+	($h_total{$tp}{$sv}{$value1}{$it}{$checks[3]}[-1])
#					-	($h_total{$tp}{$sv}{$value1}{$it}{$checks[4]}[-1]);
#					
#					push @line0,$overlap;
#					
#				}
#
#			}
#			if(@line0>5){
#				push @out0,[@line0];	
#			}
#			
#		}
#	
#	}
#
#}
#	my $file_out = $prefix.'_'.$out_suf0.'.csv';
#	open my $fh, '>', $file_out or die "couldn't open $file_out: $!";
#	print $fh join(',',@$_) . "\n" for (@out0);

