use warnings;
use strict;

my $dir = '/Users/singhal/Desktop/activeWork/exomeCaptureArray/';
my $depth = $dir . 'lineageSNPs/';
my @species = qw(Carlia_N Sapro_C);
my $qual = 50;

foreach my $lineage (@species) {
	my $seq = $depth . $lineage . '_trinity.fa.final.annotated';
	my $genus = $1 if $lineage =~ m/(\S+)_/;
	
	my $seqInfo = parseSeq($seq);
	my %seq = %{$seqInfo};
		
	my $out = $dir . $genus . '.fst.out';
	my $vcf = $depth . 'joint' . $genus . '.vcf.out';
	open(OUT, ">$out");
	print OUT "gene\ttype\tFst\n";
	
	my %var;
	open(IN, "<$vcf");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/contig/) {
			my @d = split(/\t/,$line);
			if ($d[5] >= $qual) {
				for (my $i = 9; $i < 14; $i++) {
					$var{$d[0]}{$d[1]}{'1'}++ if $d[$i] =~ m/^1\//;
					$var{$d[0]}{$d[1]}{'1'}++ if $d[$i] =~ m/^[0|1]\/1/;
					}
				for (my $i = 14; $i < 19; $i++) {
					$var{$d[0]}{$d[1]}{'2'}++ if $d[$i] =~ m/^1\//;
					$var{$d[0]}{$d[1]}{'2'}++ if $d[$i] =~ m/^[0|1]\/1/;
					}	
				}
			}
		}	
		
	foreach my $c (keys %var) {
		my @fst;
		foreach my $pos (keys %{$var{$c}}) {
			if ($pos >= $seq{$c}{'ge'}) {
				my ($a1, $a2);
				if ($var{$c}{$pos}{'1'}) {
					$a1 = $var{$c}{$pos}{'1'};
					}
				else {
					$a1 = 0;
					}
				if ($var{$c}{$pos}{'2'}) {
					$a2 = $var{$c}{$pos}{'2'};
					}
				else {
					$a2 = 0;
					}
				my $fis1 = 2 * ( $a1 / 10 ) * ( (10 - $a1) / 10 );	
				my $fis2 = 2 * ( $a2 / 10 ) * ( (10 - $a2) / 10 );	
				my $fis = ($fis1 + $fis2) / 2;
				my $a = $a1 + $a2;
				my $fit  = 2 * ( $a / 20 ) * ( (20 - $a) / 20 );
				my $fst;
				if ($fit == 0) {
					$fst = 0;
					}
				else {	
					$fst = ($fit - $fis) / $fit;
					}	
				push(@fst,$fst);
				}
			}
		my $final_fst = sprintf("%.3f", average(\@fst));
			
		print OUT $seq{$c}{'match'}, "\t", "UTR", "\t", $final_fst, "\n";
		
		@fst = '';
		foreach my $pos (keys %{$var{$c}}) {
			if ($pos <= $seq{$c}{'ge'} && $pos >= $seq{$c}{'gs'}) {
				my ($a1, $a2);
				if ($var{$c}{$pos}{'1'}) {
					$a1 = $var{$c}{$pos}{'1'};
					}
				else {
					$a1 = 0;
					}
				if ($var{$c}{$pos}{'2'}) {
					$a2 = $var{$c}{$pos}{'2'};
					}
				else {
					$a2 = 0;
					}
				my $fis1 = 2 * ( $a1 / 10 ) * ( (10 - $a1) / 10 );	
				my $fis2 = 2 * ( $a2 / 10 ) * ( (10 - $a2) / 10 );	
				my $fis = ($fis1 + $fis2) / 2;
				my $a = $a1 + $a2;
				my $fit  = 2 * ( $a / 20 ) * ( (20 - $a) / 20 );
				my $fst;
				if ($fit == 0) {
					$fst = 0;
					}
				else {	
					$fst = ($fit - $fis) / $fit;
					}	
				push(@fst,$fst);
				}
			}
		$final_fst = sprintf("%.3f", average(\@fst));
			
		print OUT $seq{$c}{'match'}, "\t", "CDS", "\t", $final_fst, "\n";
		}	
	close(OUT);	
	}

sub average {
	my ($a) = @_;
	my @a = @{$a};
	my $avg;
	if (scalar(@a) > 0) {
		my $sum = 0; my $num = scalar(@a);
		foreach my $var (@a) {
			if ($var) {
				$sum += $var;
				}
			}	
		$avg = $sum/$num;
		}
	else {
		$avg = 0;
		}
	return($avg);
	}
sub parseSeq {
	my ($s) = @_;
	my %s;
	open(IN, "<$s");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+).*ENS.*/) {
			my $c = $1; my @d = split(/\t/,$line);
			chomp(my $seq = <IN>);
			$s{$c}{'seq'} = $seq;
			$s{$c}{'length'} = length($seq);
			$s{$c}{'match'} = $d[2];
			
			if ($d[1] =~ m/5u(\d+)/) {
				$s{$c}{'5u'} = $1;
				}
			if ($d[1] =~ m/gs(\d+)/) {
				$s{$c}{'gs'} = $1;
				}
			if ($d[1] =~ m/ge(\d+)/) {
				$s{$c}{'ge'} = $1;
				}
			if ($d[1] =~ m/3u(\d+)/) {
				$s{$c}{'3u'} = $1;	
				}	
			$s{$c}{'ulength'} = $s{$c}{'length'} - $s{$c}{'ge'};	
			$s{$c}{'cdslength'} = 	$s{$c}{'ge'} - $s{$c}{'gs'};
			}
		}
	close(IN);	
	return(\%s);	
	}	