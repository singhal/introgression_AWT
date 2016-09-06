use warnings;
use strict;

my $file = shift;
my @seqfiles = </Users/singhal/Desktop/activeWork/exomeCapture/seqFiles/*annotated>;

my %seq;
open(IN, "<$file");
while(<IN>) {
	chomp(my $line = $_);
	if ($line =~ m/>(ENS\S+)/) {
		$seq{$1}++;
		}
	}
close(IN);
	
my %match;	
foreach my $seq (@seqfiles) {
	open(IN, "<$seq");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/(ENS\S+)/) {
			if ($seq{$1}) {
				my $match = $1;
				my $gs = $1 if $line =~ m/gs(\d+)/;
				my $ge = $1 if $line =~ m/ge(\d+)/;
				my $length = $ge - $gs + 1;
				$gs = $gs - 1;
				
				my $id = $1 if $seq =~ m/([A-Z|\_]+)_trinity/i;
				my $old_id = $1 if $line =~ m/(>\S+)/;
				my $new_id = $old_id . '_' . $id;
				$line =~ s/$old_id/$new_id/;
				
				
				chomp(my $s = <IN>);
				$s = substr $s, $gs, $length;
			
				$match{$match}{$line} = $s;
				}
			}
		}
	}	
	
foreach my $match (keys %match) {
	foreach my $id (keys %{$match{$match}}) {
		print $id, "\n", $match{$match}{$id}, "\n";
		}
	}	