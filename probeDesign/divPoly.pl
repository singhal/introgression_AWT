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
		
	my $out = $dir . $genus . '.divpoly.out';
	my $vcf = $depth . 'joint' . $genus . '.vcf.out';
	open(OUT, ">$out");
	print OUT "gene\ttype\tpoly1\tpoly2\trawdiv\tnetdiv\tlength\n";
	
	my %var;
	open(IN, "<$vcf");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/contig/) {
			my @d = split(/\t/,$line);
			if ($d[5] >= $qual) {
				for (my $i = 9; $i < 19; $i++) {
					my $geno1 = $1 if $d[$i] =~ m/([0|1])\//;
					my $geno2 = $1 if $d[$i] =~ m/\/([0|1])/;
					$var{$d[0]}{$d[1]}{$i} = [$geno1,$geno2];
					}
				}
			}
		}	
		
	foreach my $c (keys %var) {
		my @haplo;
		for (my $i = 9; $i < 19; $i++) {
			for (my $j = 0; $j < 2; $j++) {
				my $haplo;
				foreach my $pos (keys %{$var{$c}}) {
					if ($pos >= $seq{$c}{'ge'}) {
						my $rand = rand();
						if ($rand > 0.5) {
							$haplo .= $var{$c}{$pos}{$i}[1];
							}
						else {
							$haplo .= $var{$c}{$pos}{$i}[0];
							}
						}	
					}
				push(@haplo,$haplo);
				}
			}
			
		if ($haplo[0]) {	
			my @h1 = @haplo[0..9];	
			my @h2 = @haplo[10..19];
			#now time to calculate polymorphism and divergence
			my $within1 =  sprintf("%.4f",varwithin(\@h1)/$seq{$c}{'ulength'});
			my $within2 = sprintf("%.4f",varwithin(\@h2)/$seq{$c}{'ulength'});
			my $raw = sprintf("%.4f",varbtn(\@h1,\@h2)/$seq{$c}{'ulength'});
			my $net = $raw - ($within1 + $within2) / 2;
		
			print OUT $seq{$c}{'match'}, "\tUTR\t", $within1, "\t", $within2, "\t", $raw, "\t", $net, "\t$seq{$c}{'ulength'}\n";
			}
		else {
			print OUT $seq{$c}{'match'}, "\tUTR\t", "0", "\t", "0", "\t", "0", "\t", "0", "\t$seq{$c}{'ulength'}\n";
			}
		}
		
	foreach my $c (keys %var) {
		my @haplo;
		for (my $i = 9; $i < 19; $i++) {
			for (my $j = 0; $j < 2; $j++) {
				my $haplo;
				foreach my $pos (keys %{$var{$c}}) {
					if ($pos <= $seq{$c}{'ge'} && $pos >= $seq{$c}{'gs'}) {
						my $rand = rand();
						if ($rand > 0.5) {
							$haplo .= $var{$c}{$pos}{$i}[1];
							}
						else {
							$haplo .= $var{$c}{$pos}{$i}[0];
							}
						}	
					}
				push(@haplo,$haplo);
				}
			}
			
		if ($haplo[0]) {	
			my @h1 = @haplo[0..9];	
			my @h2 = @haplo[10..19];
			#now time to calculate polymorphism and divergence
			my $within1 =  sprintf("%.4f",varwithin(\@h1)/$seq{$c}{'cdslength'});
			my $within2 = sprintf("%.4f",varwithin(\@h2)/$seq{$c}{'cdslength'});
			my $raw = sprintf("%.4f",varbtn(\@h1,\@h2)/$seq{$c}{'cdslength'});
			my $net = $raw - ($within1 + $within2) / 2;
		
			print OUT $seq{$c}{'match'}, "\tCDS\t", $within1, "\t", $within2, "\t", $raw, "\t", $net, "\t$seq{$c}{'cdslength'}\n";
			}
		else {
			print OUT $seq{$c}{'match'}, "\tCDS\t", "0", "\t", "0", "\t", "0", "\t", "0", "\t$seq{$c}{'cdslength'}\n";
			}
		}
	close(OUT);	
	}

sub varwithin {
	my ($haplo) = @_;
	my @h = @{$haplo};
	my $sum; my $compare;
	for (my $i = 0; $i  < scalar(@h); $i++) {
		for (my $j = 0; $j  < scalar(@h); $j++) {
			my $diff = 0;
			my $mask = $h[$i] ^ $h[$j];
			$diff++ while ($mask =~ /[^\0]/g);
			$sum += $diff;
			$compare++;
			}
		}
	$sum = $sum/$compare;	
	return($sum);
	}


sub varbtn {
	my ($haplo1,$haplo2) = @_;
	my @h1 = @{$haplo1};
	my @h2 = @{$haplo2};
	my $sum; my $compare;
	for (my $i = 0; $i  < scalar(@h1); $i++) {
		for (my $j = 0; $j  < scalar(@h2); $j++) {
			my $diff = 0;
			my $mask = $h1[$i] ^ $h2[$j];
			$diff++ while ($mask =~ /[^\0]/g);
			$sum += $diff;
			$compare++;
			}
		}
	$sum = $sum/$compare;	
	return($sum);
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