use warnings;
use strict;

my $genus = 'Sapro';
my $dir = '/Users/singhal/thesisWork/exomeCapture/preRepeat_' . $genus . '/';
my $start = $dir . $genus . '_preRepeatTargets.fa';
my $self = $dir . 'self.blast.out';
my $out =  $dir . 'repeat.blast.out';
my $rm = $dir . 'repeatMasker.out';

my %dup;

open(IN, "<$self");
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/, $line);
	$d[0] =~ s/_\S+$//i;
	$d[1] =~ s/_\S+$//i;
	
	if ($d[0] ne $d[1]) {
		$dup{$d[0]}++;
		$dup{$d[1]}++;
		}
	}
close(IN);	
	
open(IN, "<$out");
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/, $line);
	$d[0] =~ s/_\S+$//i;
	$dup{$d[0]}++;
	}
close(IN);	

open(IN, "<$rm");
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\s+/, $line);
	if ($line =~ m/ENS/) {
		if ($d[1] > 500) {
			$d[5] =~ s/_\S+$//i;
			$dup{$d[5]}++;
			}
		}	
	}
close(IN);

open(IN, "<$start");
while(<IN>) {
	chomp(my $line = $_);
	if ($line =~ m/>(\S+)/) {
		my $id = $1; 
		my $short = $id;
		$short =~ s/_\S+$//i;
		unless ($dup{$short}) {
			print $line, "\n";
			my $seq = <IN>;
			print $seq;
			}
		}
	}
close(IN);	