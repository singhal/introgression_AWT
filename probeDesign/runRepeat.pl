use warnings;
use strict;

my @species = qw(Sapro Lampro1 Lampro2 Carlia);
my $dir = '/Users/singhal/thesisWork/exomeCapture/';
my %homologs;
foreach my $species (@species) {
	my $out = $dir . 'preRepeat_' . $species . '/' . $species . '_preRepeatTargets.fa';	
	
	my $call1 = system("formatdb -i $out -p F");
	my $self = $dir . 'preRepeat_' . $species . '/self.blast.out'; 
	my $call2 = system("blastall -p blastn -d $out -i $out -m 8 -a 4 -o $self -e 1e-20");
	my $rep = $dir . 'preRepeat_' . $species . '/repeat.blast.out'; 
	my $repeat = $dir . 'repeatDatabase/repeatVertebrate.fa';
	my $call3 = system("blastall -p blastn -d $repeat -i $out -m 8 -a 4 -o $rep -e 1e-20");
	}