use warnings;
use strict;

my $dir = '/Users/singhal/Desktop/activeWork/exomeCapture/';
my $file = $dir . 'preRepeat/divergence.out';

my %exons;
open(IN, "<$file");
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/, $line);
	$exons{$d[0]}++;
	}
close(IN);	

foreach my $exon (keys %exons) {
	my $aln = $dir . 'results/' . $exon . '.fa.aln';
	open(IN, "<$aln");
	
	my %s; my $id;
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			$id = $1;
			}
		else {
			$s{$id} .= $line;
			}
		}
	close(IN);	
	
	foreach my $s (keys %s) {
		my $seq = $s{$s};
		$seq =~ s/-//g;
		$s{$s} = $seq;
		}
			
	my %g;
	foreach my $s (keys %s) {
		my $g = $1 if $s =~ m/(\S+)_/;
		if ($g{$g}) {
			$g{$g} = $s{$s} if length($s{$s}) > length($g{$g});
			}
		else {
			$g{$g} = $s{$s};
			}
		}
	close(IN);
	
	if (scalar(keys %g) > 2) {
		foreach my $g (keys %g) {
			print ">", $exon, "_", $g, "\n";
			print $g{$g}, "\n";
			}
		}
	}	