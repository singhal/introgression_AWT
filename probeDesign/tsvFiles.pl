use warnings;
use strict;

my $seq = shift;
open(IN, "<$seq");
my $out1 = $seq . ".renamed";
my $out2 = $seq . ".tsv";
open(OUT1, ">$out1");
open(OUT2, ">$out2");
my $tracker = 1;
while(<IN>) {
	chomp(my $line = $_);
	if ($line =~ m/>(\S+)/) {
		print OUT1 ">chr", $tracker, "_$1\n";
		chomp(my $seq = <IN>);
		print OUT1 $seq, "\n";
		my $length = length($seq);
		print OUT2 "$tracker\tchr" . $tracker . "_$1\t1\t$length\n";
		$tracker++;
		}
	}
close(OUT1); close(OUT2); close(IN);	