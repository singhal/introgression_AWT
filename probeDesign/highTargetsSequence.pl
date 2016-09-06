use warnings;
use strict;

my @files = </Users/singhal/Desktop/activeWork/exomeCapture/finalResults/high*out>;
my @seqfiles = </Users/singhal/Desktop/activeWork/exomeCapture/seqFiles/*annotated>;


my %seq;
foreach my $file (@files) {
	open(IN, "<$file");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/(ENS\S+)/) {
			$seq{$1}++;
			}
		}
	close(IN);
	}
	
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
				
				my $id = $1 if $seq =~ m/([A-Z|\_]+)_[a-z]_trinity/i;
				my $old_id = $1 if $line =~ m/(>\S+)/;
				my $new_id = $old_id . '_' . $id;
				$line =~ s/$old_id/$new_id/;
				
				
				chomp(my $s = <IN>);
				$s = substr $s, $gs, $length;
				
				if ($match{$match}{$id}) {
					$match{$match}{$id} = $s if length($s) > length($match{$match}{$id});
					}
				else {
					$match{$match}{$id} = $s;
					}
				}
			}
		}
	}	
	
foreach my $match (keys %match) {
	my $num = scalar(keys %{$match{$match}});
	foreach my $id (keys %{$match{$match}}) {
		print ">", $match, "_", $id, "\n", $match{$match}{$id}, "\n";
		}
	}		