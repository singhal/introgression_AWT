use warnings;
use strict;

#1. 2 columns are required. 1st is ID of the probe, 2nd is the sequence of the probe
#2. the probes should be randomized before uploading.

my $file = shift;
my %probes;

open(IN, "<$file");
while(<IN>) {
	chomp(my $line = $_);
	my @d = split(/\t/,$line);
	$probes{$d[0]} = $d[3];
	}

my @probes = sort {$a cmp $b} keys %probes;
fisher_yates_shuffle(\@probes);

foreach my $id (@probes) {
	print $id, "\t", $probes{$id}, "\n";
	}

sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
	    }
	}	