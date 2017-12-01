#!/usr/bin/perl

use warnings;
use strict;

my $fn='1'.'.txt';
unless (-e $fn){
	print "$fn does not exist.\n";
	next;
}

my $val=5;

my @csSet=("C","S");
my @idxSet=(0..30);
my %allHash;
my @tsSet=(1..$val);
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
#my %ttHash;
#map{chomp($_);push @total,(split /\s*,\s*/,$_); ++$ttHash{$_};++$allHash{$_};} <$fh>;
map{chomp($_);push @total,(split /\s*,\s*/,$_); ++$allHash{$_};} <$fh>;
close ($fh);
my $idx=0;
my @cs = grep {++$idx %2} @total; 
my %csHash;
map {$csHash{$_}++}@cs;

for my $k (keys(%csHash)){
	if ($csHash{$k}!=$val){
		print "lost: $k => $csHash{$k}\n";
	}
}
foreach my $k (keys %allHash){
	if($allHash{$k}>1){
		print "repeat: $k => $allHash{$k} \n";
	}
}

foreach my $k(keys %allHash){
	if($allHash{$k}==0){
		print "miss: $k \n";
	}
}
