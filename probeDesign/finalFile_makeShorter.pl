use warnings;
use strict;

my $start = shift;
my $min = 300;
my $max = 5000;
my $lMax = 2.99e6;

my %seq;
open(IN, "<$start");
while(<IN>) {
	chomp(my $line = $_);
	if ($line =~ m/>(\S+)/) {
		my $id = $1;
		my $gene = $1 if $id =~ m/(\S+)_/;
		chomp(my $line = <IN>);
		$seq{$gene}{$id} = $line;
		}
	}	
	
my $length = 0;	
foreach my $c (keys %seq) 	{
	if ($length < $lMax) {
		my @id = keys %{$seq{$c}};
		my $id = $id[0];
		if (length($seq{$c}{$id}) > $min && length($seq{$c}{$id}) < $max) {
			foreach my $i (keys %{$seq{$c}}) {
				print ">", $i, "\n", $seq{$c}{$i}, "\n";
				$length += length($seq{$c}{$i});
				}
			}
		}
	}	