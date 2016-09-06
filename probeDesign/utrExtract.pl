use warnings;
use strict;

#a script to get the UTRs

my $dir = '/Users/singhal/Desktop/activeWork/exomeCapture/';
my @rflp = qw(carliaRFLP gilliesRFLP nbmbRFLP saproRFLP);
my @seq = qw(Carlia_N Lampro Lampro_N Sapro_C);

for (my $i = 0; $i < scalar(@rflp); $i++) {
	my $rflp = $dir . 'rflp/' . $rflp[$i];
	my $seq = $dir . 'seqFiles/' . $seq[$i] . '_trinity.fa.final.annotated';
	my $id = $rflp[$i];
	$id =~ s/RFLP//;
	
	my %rflp;
	open(IN, "<$rflp");
	my $junk = <IN>;
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		$rflp{$d[0]}{'pos'} = $d[1];
		$rflp{$d[0]}{'match'} = $d[2];
		}
	close(IN);
	
	open(IN, "<$seq");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $c = $1;
			if ($rflp{$c}) {
				my $ge = $1 if $line =~ m/ge(\d+)/;
				my $match = $1 if $line =~ m/(ENS\S+)/;
				chomp(my $s = <IN>);
				print ">", $rflp{$c}{'match'}, "_utr_", $id, "\n";
				if ($rflp{$c}{'pos'} > $ge) {
					my $potStart = $rflp{$c}{'pos'} - 150;
				 	if ($potStart > $ge) {
				 		my $subseq = substr $s, $potStart, 300;
				 		print $subseq, "\n";
				 		}
				 	else {
				 		my $subseq = substr $s, $ge, 300;
				 		print $subseq, "\n";
				 		}
				 	}
				 else {
				 	my $start = $rflp{$c}{'pos'} - 150;
				 	my $subseq = substr $s, $start, 300;
				 	print $subseq, "\n";
				 	}
				 }
			}
		}
	}	
		