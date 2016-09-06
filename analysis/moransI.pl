use warnings;
use strict;

my @contact = qw(sjo nBMB gillies carlia);
#my @contact = qw(carlia);
my @dist = qw(0 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000);

foreach my $contact (@contact) {
	my $file = "/Users/sonal/thesisWork/introgression/LD/LD_" . $contact . ".out";

	my %d;
	open(IN,"<$file");
	while(<IN>) {
		chomp(my $line = $_);
		
		my @d = split(/\t/,$line);
		push(@{$d{$d[6]}},\@d); 
		# push(@{$d{$d[2]}},\@d); 
	
		}
	close(IN);	

	for (my $n = 0; $n < scalar(@dist) - 1; $n++) {
	
		my $nLoci = 0;
		my $sumNumerator = 0;
		my @numerator = 0;
		my $nComparisons = 0;
		my $sumDenominator = 0;
		my @denominator = 0;
	
		my $dist1 = $dist[$n];
		my $dist2 = $dist[$n + 1];
	
		foreach my $c (keys %d) {
			my @tmp = @{$d{$c}};
			if (scalar(@tmp) > 1) {
				$nLoci += scalar(@tmp);  
			
				for (my $i = 0; $i < scalar(@tmp); $i++) {
					my $value = $tmp[$i][4] ** 2;
					$sumDenominator += $value;
					push(@denominator, $value);
					}
			
				for (my $j = 0; $j < scalar(@tmp); $j++) {
					for (my $k = $j + 1; $k < scalar(@tmp); $k++) {
						# switch back to 3 for more local version
						if ( abs($tmp[$j][7] - $tmp[$k][7]) >= $dist1) {
							if ( abs($tmp[$j][7] - $tmp[$k][7]) <= $dist2) {
								$nComparisons++;
								my $value = $tmp[$j][4] * $tmp[$k][4];
								$sumNumerator += $value;
								push(@numerator, $value);
								}
							}	
						}			
					}	
				}
			}	
		if ($nComparisons > 0) {
			my $moranI = ($nLoci * $sumNumerator) / ($nComparisons * $sumDenominator);
			print $contact, "\t", $dist1, "\t", $moranI, "\t$nComparisons\n";
			}
		}
	}	
	