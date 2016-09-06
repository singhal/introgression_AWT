use warnings;
use strict;

my @species = qw(Sapro Lampro1 Lampro2 Carlia);
my $dir = '/Users/singhal/thesisWork/exomeCapture/';
my $targetL = 3.6e6; #gives 1.8 unique MB
my %compare = (	'Sapro' => ['Sapro_C', 'Sapro_S'],
				'Lampro1' => ['Lampro_N', 'Lampro_C'], 
				'Lampro2' => ['Lampro_C', 'Lampro_S'], 
				'Carlia' => ['Carlia_N', 'Carlia_S']);

my %homologs;
foreach my $species (@species) {
	my $div = $dir . 'preRepeat_' . $species . '/' . $species . '_divergence.out';
	open(IN, "<$div");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/(ENS\S+)/) {
			my $id = $1;
			unless ($line =~ m/exon/) {
				$homologs{$id}++;
				}
			}
		}
	close(IN);	
	}	
foreach my $h (keys %homologs) {
	unless ($homologs{$h} == 4) {
		delete($homologs{$h});
		}
	}
my @homologs = sort {$a cmp $b} keys %homologs;
fisher_yates_shuffle(\@homologs);

foreach my $species (@species) {
	my $length = 0;
	my %exons;
	my $outliers = $dir . 'preRepeat_' . $species . '/' . $species . '.outliers.out';
	my $candidate = $dir . 'preRepeat_' . $species . '/candidateGenes.fa';
	my $mito = $dir . 'preRepeat_' . $species . '/mitochondrialExons.fa';
	my $utr = $dir . 'preRepeat_' . $species . '/finalUTR.fa';
	my $div = $dir . 'preRepeat_' . $species . '/' . $species . '_divergence.out';
	
	my %cds;
	foreach my $seq (@{$compare{$species}}) {
		my $seqFile = $dir . 'seqFiles/' . $seq . '_trinity.fa.final.annotated';
		open(IN, "<$seqFile");
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/(ENS\S+)/) {
				my $match = $1;
				my $gs = $1 if $line =~ m/gs(\d+)/;
				my $ge = $1 if $line =~ m/ge(\d+)/;
				my $length = $ge - $gs + 1;
				$gs = $gs - 1;
				
				chomp(my $s = <IN>);
				$s = substr $s, $gs, $length;
				
				$cds{$match}{$seq} = $s;		
				}
			}
		close(IN);	
		}
		
	open(IN, "<$utr");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $id = $1;
			chomp(my $seq = <IN>);
			push(@{$exons{$id}}, $seq);
			$length += length($seq);
			}
		}
	close(IN);	
	
	open(IN, "<$mito");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $id = $1;
			chomp(my $seq = <IN>);
			push(@{$exons{$id}}, $seq);
			$length += length($seq);
			}
		}
	close(IN);
	
	my %targets;
	open(IN, "<$outliers");
	while(<IN>) {
		chomp(my $id = $_);
		$targets{$id}++;
		}
	close(IN);	
	
	open(IN, "<$candidate");
	while(<IN>) {
		chomp(my $id = $_);
		if ($id =~ m/(ENSAC\S+)/) {
			$targets{$1}++;
			}
		}
	close(IN);
	
	foreach my $target (keys %targets) {
		my $aln = $dir . 'results_' . $species . '/' . $target . '.fa.aln';
	
		my %seq; my $id;
		
		if (-f $aln) {		
			open(IN, "<$aln");
			while(<IN>) {
				chomp(my $line = $_);
				if ($line =~ m/>(\S+)/) {
					$id = $1;
					}
				else {
					$line =~ s/-//g;
					$seq{$id} .= $line;
					}
				}
			close(IN);	
			}
		else {
			if ($cds{$target}) {
				foreach my $i (keys %{$cds{$target}}) {
					$seq{$i} = $cds{$target}{$i};
					}
				}
			}
			
		if (scalar(keys %seq) == 1) {
			foreach my $s (keys %seq) {
				push(@{$exons{$target}}, $seq{$s});
				$length += length($seq{$s});
				}
			foreach my $s (keys %seq) {
				push(@{$exons{$target}}, $seq{$s});
				$length += length($seq{$s});
				}			
			}
		elsif (scalar(keys %seq) == 2) {
			foreach my $s (keys %seq) {
				push(@{$exons{$target}}, $seq{$s});
				$length += length($seq{$s});
				}
			}		
		}
		
	open(IN, "<$div");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/(ENS\S+exon\S+)/) {
			my $target = $1;
			my $aln = $dir . 'results_' . $species . '/' . $target . '.fa.aln';
				
			my %seq; my $id;
		
			if (-f $aln) {
				open(ALN, "<$aln");
				while(<ALN>) {
					chomp(my $line = $_);
					if ($line =~ m/>(\S+)/) {
						$id = $1;
						}
					else {
						$line =~ s/-//g;
						$seq{$id} .= $line;
						}
					}
				close(ALN);	
				}
			else {
				if ($cds{$target}) {
					foreach my $i (keys %{$cds{$target}}) {
						$seq{$i} = $cds{$target}{$i};
						}
					}
				}
							
			my %n;
			$n{$_}++ foreach @{$compare{$species}};	
			foreach my $c (keys %seq) {
				unless($n{$c}) {
					delete($seq{$c});
					}
				}	
			
			my $shortTarget = $1;
			$shortTarget =~ s/_\S+$//;
			unless($exons{$shortTarget}) {
				if (scalar(keys %seq) == 1) {
					foreach my $s (keys %seq) {
						push(@{$exons{$target}}, $seq{$s});
						$length += length($seq{$s});
						}
					foreach my $s (keys %seq) {
						push(@{$exons{$target}}, $seq{$s});
						$length += length($seq{$s});
						}			
					}
				elsif (scalar(keys %seq) == 2) {
					foreach my $s (keys %seq) {
						push(@{$exons{$target}}, $seq{$s});
						$length += length($seq{$s});
						}
					}	
				}
			}	
		}
	close(IN);
		
	foreach my $target (@homologs) {
		my $shortTarget = $target;
		$shortTarget =~ s/_\S+$//;
		unless($exons{$shortTarget}) {
			if ($length < $targetL) {
		
				
				my $aln = $dir . 'results_' . $species . '/' . $target . '.fa.aln';						
				my %seq; my $id;
		
				if (-f $aln) {
					open(ALN, "<$aln");
					while(<ALN>) {
						chomp(my $line = $_);
						if ($line =~ m/>(\S+)/) {
							$id = $1;
							}
						else {
							$line =~ s/-//g;
							$seq{$id} .= $line;
							}
						}
					close(ALN);	
					}
				else {
					if ($cds{$target}) {
						foreach my $i (keys %{$cds{$target}}) {
							$seq{$i} = $cds{$target}{$i};
							}
						}
					}
							
				my %n;
				$n{$_}++ foreach @{$compare{$species}};	
				foreach my $c (keys %seq) {
					unless($n{$c}) {
						delete($seq{$c});
						}
					}	
			
				if (scalar(keys %seq) == 1) {
					foreach my $s (keys %seq) {
						push(@{$exons{$target}}, $seq{$s});
						$length += length($seq{$s});
						}
					foreach my $s (keys %seq) {
						push(@{$exons{$target}}, $seq{$s});
						$length += length($seq{$s});
						}			
					}
				elsif (scalar(keys %seq) == 2) {
					foreach my $s (keys %seq) {
						push(@{$exons{$target}}, $seq{$s});
						$length += length($seq{$s});
						}
					}	
				}	
			}
		}

	my $out = $dir . 'preRepeat_' . $species . '/' . $species . '_preRepeatTargets.fa';		
	open(OUT, ">$out");
	foreach my $target (sort {$a cmp $b} keys %exons) {
		for (my $i = 0; $i < scalar(@{$exons{$target}}); $i++) {
			print OUT ">", $target, "_", $i, "\n", $exons{$target}[$i], "\n";
			}
		}	
	close(OUT);	
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