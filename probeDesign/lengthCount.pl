use warnings;
use strict;

my $file = shift;

my $length;
open(IN, "<$file");
while(<IN>) {
	unless ($_ =~ m/>/) {
		chomp(my $line = $_);
		$length += length($line);
		}
	}
print $length, "\n";	