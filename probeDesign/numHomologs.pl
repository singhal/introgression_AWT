use warnings;
use strict;

my @species = qw(Sapro Lampro1 Lampro2 Carlia);
my $dir = '/Users/singhal/thesisWork/exomeCapture/';
my %homologs;
foreach my $species (@species) {
	my $out = $dir . 'preRepeat_' . $species . '/' . $species . '_preRepeatTargets.fa';	
	open(IN, "<$out");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/(ENS\S+)_[1|0]/) {
			my $id = $1;
			$homologs{$id}++;
			}
		}
	close(IN);	
	}	
foreach my $h (keys %homologs) {
	unless ($homologs{$h} == 8) {
		delete($homologs{$h});
		}
	}
print $_, "\n" foreach keys %homologs;	