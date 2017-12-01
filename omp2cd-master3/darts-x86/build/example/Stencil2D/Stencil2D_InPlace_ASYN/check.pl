#!/usr/bin/perl

use warnings;
use strict;


my $fn = $ARGV[0];
unless (-e $fn){
	print "$fn does not exist.\n";
	next;
}

#open my $fh, '<',$ARGV[0] or die $!;

#my $val=3;
my $Iter = $ARGV[1];
my $Repeat = ($ARGV[2])?($ARGV[2]):1;
my $val = $Iter*$Repeat;
my $thrd = ($ARGV[3])?($ARGV[3]):30;
my $thrdMus1 = $thrd -1;
my @csSet=("U","D","C");
my @idxSet=(0..$thrdMus1);
my %allHash;
my @tsSet=(1..$Iter);
for my $csVar (@csSet){
	for my $csIdx(@idxSet){
		for my $tsIdx (@tsSet){
			my $k = $csVar.$csIdx.",ts:".$tsIdx;
			$allHash{$k}=0;
		}
	}
}

open my $fh, '<',$fn or die "can't open $fn: $!";
print "Processing $fn...\n";
my @total;
my @totalR;
#my %ttHash;
#map{chomp($_);push @total,(split /\s*,\s*/,$_); ++$ttHash{$_};++$allHash{$_};} <$fh>;
#map{chomp($_);push @total,(split /\s*,\s*/,$_); ++$allHash{$_};} <$fh>;
map{chomp($_); push @total,(split /\s*,\s*/,$_);push @totalR,[(split /\s*,\s*/,$_)] ; $allHash{$_}++; } <$fh>;
close ($fh);

my $idx=0;
#my @cs = grep {++$idx %2} @total; 
my @cs = map {$_->[0]} @totalR; 
my %csHash;
map {$csHash{$_}++}@cs;

#my $size = keys %allHash;
#print $size, "\n";
#for my $k (keys(%allHash)){
#	print "$k => $allHash{$k};";
#}
for my $k (keys(%csHash)){
	if ($csHash{$k}!=$val){
		print "wrong: $k => $csHash{$k}\n";
	}
}
foreach my $k (keys %allHash){
	if($allHash{$k}>$Repeat){
		print "repeat: $k => $allHash{$k} \n";
	}
}

foreach my $k(keys %allHash){
	if($allHash{$k}!=$Repeat){
		print "miss: $k \n";
	}
}
