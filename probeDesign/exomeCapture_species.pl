use warnings;
use strict;
use Bio::DB::Fasta;

my $dir = '/Users/singhal/Desktop/activeWork/exomeCapture/';
my $results_dir = $dir . 'results_Sapro/';
#consider all exons of this length or greater
my $exonLength = 200;
#consider all exons which are matched at least this much or greater
my $exonMatch = 0.25; 
my $anolisExons = $dir . 'anolisExons.fa';
my %compare2 = ('Sapro_C' => 'Sapro_S');
my $genus = 'Sapro';
my $minGC = 0.3;
my $maxGC = 0.8;
my $np = 4;
my $eval = 1e-20;
my $parse = 1;
my $minDiv = 0.001;
my $maxDiv = 0.05;

##########################
# start the subroutines  #
##########################

my $divpoly = $dir . 'cladeResults/' . $genus . ".divpoly.out";
my $anolisOut = $anolisExons . ".unique";
mkdir($results_dir) unless(-d $results_dir);
#cdsBlast();
#alignHomologs();
my $div = divpoly();
my $info = calculateDiv($div);
finalDiv($info,$div);

##########################
# behold the subroutines #
##########################

sub alignHomologs {
	my @seq = <$results_dir*ENS*fa>;
	
	foreach my $file (@seq) {
		my $out = $file . ".aln";
		my $call1 = system("muscle -in $file -out $out");
		my $call2 = system("rm $file") if (-f $out);
		}		
	}
	
sub finalDiv {
	my ($final,$divpoly) = @_;
	my %final = %{$final};
	my %divpoly = %{$divpoly};
	foreach my $lineage (keys %compare2) {
		my $out1 = $dir . 'seqFiles/' . $lineage . '_cds.fa';
		
		my %seq;
		open(IN, "<$out1");
		while(<IN>) {
			if ($_ =~ m/>(\S+)/) {
				my $id = $1;
				chomp(my $seq = <IN>);
				$seq{$id} = $seq;
				}
			}
		close(IN);	
		
		my @exons = <$results_dir*exon*aln>;
		
		foreach my $c (keys %seq) {
			my $new_c = $1 if $c =~ m/^([A-Z|0-9]+)_/ig;
			unless($final{$new_c}) {
				my $seq = $seq{$c};
				my $gc = ($seq =~ s/[gc]//ig);
				$gc = sprintf("%.3f",$gc/length($seq{$c}));
				if ($divpoly{$new_c}) {
					$final{$new_c} = {'length' => length($seq{$c}), 'gc' => $gc, 'div' => 'NA', 'rawdiv' => $divpoly{$new_c}{'raw'}, 'netdiv' => $divpoly{$new_c}{'net'}}; 
					}
				else {
					$final{$new_c} = {'length' => length($seq{$c}), 'gc' => $gc, 'div' => 'NA', 'rawdiv' => 'NA', 'netdiv' => 'NA'}; 
					}
				}
			}
			
		foreach my $exon (@exons) {
			my $name = $1 if $exon =~ m/(ENS\S+exon\d+)/;
			unless($final{$name}) {
				open(IN, "<$exon");
				
				my %species; my $id;				
				while(<IN>) {
					chomp(my $line = $_);
					if ($line =~ m/>(\S+)/) {
						$id = $1;
						}
					else {
						$species{$id} .= $line;
						}
					}
			
				foreach my $id (keys %species) {
					unless($id =~ m/$genus/) {
						delete($species{$id});
						}
					}
					
				foreach my $c (keys %species) {
					my $seq = $species{$c};
					my $gc = ($seq =~ s/[gc]//ig);
					$gc = sprintf("%.3f",$gc/length($species{$c}));
					
					my $short = $name;
					$short = $1 if $name =~ m/^(\S+)_/;
					
					if ($divpoly{$short}) {
						$final{$name} = {'length' => length($species{$c}), 'gc' => $gc, 'div' => 'NA', 'rawdiv' => $divpoly{$short}{'raw'}, 'netdiv' => $divpoly{$short}{'net'}}; 
						}
					else {
						$final{$name} = {'length' => length($species{$c}), 'gc' => $gc, 'div' => 'NA', 'rawdiv' => 'NA', 'netdiv' => 'NA'};  
						}
					}
				}
			}
		}
		
	my @keys = ('length', 'gc', 'div', 'rawdiv', 'netdiv');
	print "name\t"; print join("\t",@keys); print "\n";
		
	foreach my $name (sort {$a cmp $b} keys %final) {
		if ($final{$name}{'length'} > $exonLength) {
			if ($final{$name}{'gc'} > $minGC && $final{$name}{'gc'} < $maxGC) {
				if ($final{$name}{'div'} ne 'NA') {
					if ($final{$name}{'div'} > $minDiv && $final{$name}{'div'} < $maxDiv) {
						}
					else {
						delete($final{$name});
						}
					}
				else {
					delete($final{$name});
					}
				}
			else {
				delete($final{$name});
				}
			}
		else {
			delete($final{$name});
			}
		}	
					
	my %i;					
	foreach my $name (keys %final) {					
		if ($name =~ m/exon/) {
			my $short = $1 if $name =~ m/(\S+)_/;
			$i{$short}++;
			}
		}	
	
	foreach my $name (sort {$a cmp $b} keys %final) {
		unless($i{$name}) {
			print $name, "\t";
			foreach my $key (@keys) {
				print $final{$name}{$key}, "\t";
				}
			print "\n";	
			}
		}	
	
	}	
				
sub cdsBlast {
	foreach my $lineage (keys %compare2) {
		my $seq = $dir . 'seqFiles/' . $lineage . '_trinity.fa.final.annotated';
		my $out1 = $dir . 'seqFiles/' . $lineage . '_cds.fa';
		open(OUT, ">$out1");
		open(IN, "<$seq");
		my %seq1;
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/ENS/) {
				chomp(my $seq = <IN>);
				my $match = $1 if $line =~ m/(ENS\S+)/;
				my $gs = $1 if $line =~ m/gs(\d+)/;
				my $ge = $1 if $line =~ m/ge(\d+)/;
				my $l = $ge - $gs + 1;
				my $cds = substr $seq, $gs - 1, $l;
				$seq1{$match . "_" . $lineage} = $cds;
				print OUT ">", $match, "_", $lineage, "\n", $cds, "\n";
				}
			}
		close(OUT);	
		close(IN);
		
		my $out2 = $dir . 'seqFiles/' . $compare2{$lineage} . '_trinity.fa.final.annotated.blast';
		my %seq2;
		open(IN, "<$out2");
		while(<IN>) {
			if ($_ =~ m/>(\S+)/) {
				my $id = $1;
				chomp(my $seq = <IN>);
				$seq2{$id} = $seq;
				}
			}
		close(IN);		
		
		my $call1 = system("formatdb -i $out1 -p F") unless(-f $out1 . ".nin");
		my $call2 = system("formatdb -i $out2 -p F") unless(-f $out2 . ".nin");
		my $bout1 = $out1 . ".blast.out";
		my $bout2 = $out1 . ".recipblast.out";
		my $call3 = system("blastall -p blastn -d $out1 -i $out2 -m 8 -a $np -b 1 -o $bout1 -e $eval") unless(-f $bout1);
		my $call4 = system("blastall -p blastn -d $out2 -i $out1 -m 8 -b 1 -o $bout2 -a $np -e $eval") unless (-f $bout2);
		
		my %matches;
		open(IN, "<$bout2");
		while(<IN>) {
			chomp(my $line = $_);
			my @d = split(/\t/,$line);
			$matches{$d[0]}{$d[1]} = 'no';
			}
		close(IN);
		
		open(IN, "<$bout1");
		while(<IN>) {
			chomp(my $line = $_);
			my @d = split(/\t/,$line);
			if ($matches{$d[1]}{$d[0]}) {
				$matches{$d[1]}{$d[0]} = 'yes';
				}
			}
		close(IN);
		
		foreach my $c (keys %matches) {
			foreach my $match (keys %{$matches{$c}}) {
				if ($matches{$c}{$match} eq 'yes') {
					my $out1 = 'exonHomolog.fa';
					open(OUT1, ">$out1");
					print OUT1 ">", $c, "\n" , $seq1{$c}, "\n";
					my $elength = length($seq1{$c});
					close(OUT1);
		
					my $out2 = 'exonHomolog2.fa';
					open(OUT2, ">$out2");
					print OUT2 ">", $match, "\n" , $seq2{$match}, "\n";
					close(OUT2);
		
					my @call = `exonerate -m coding2coding $out1 $out2 --showvulgar yes --showalignment no`;
		
					if ($call[3]) {
						if ($call[3] =~ m/vulgar/) {
							my @d = split(/\s/,$call[3]);
							my $length = abs($d[6] - $d[7]);					
							my $sub; 
							my $start;
			
							if ($d[4] =~ m/\+/ && $d[8] =~ m/\+/) {
								$start = $d[6];
								$sub = substr $seq2{$match}, $start, $length;
								}
							elsif ($d[4] =~ m/\+/ && $d[8] =~ m/\-/) {
								$start = $d[7];
								$sub = substr $seq2{$match}, $start, $length;
								$sub = reverse($sub);
								$sub =~ tr/ATGCatgc/TACGtacg/;
								}
							elsif ($d[4] =~ m/\-/ && $d[8] =~ m/\-/) {
								$start = $d[7];
								$sub = substr $seq2{$match}, $start, $length;
								}
							elsif ($d[4] =~ m/\-/ && $d[8] =~ m/\+/) {
								$start = $d[6];
								$sub = substr $seq2{$match}, $start, $length;
								$sub = reverse($sub);
								$sub =~ tr/ATGCatgc/TACGtacg/;
								}
				
							if ($length/$elength >= $exonMatch) {
								my $name = $1 if $c =~ m/^([A-Z|0-9]+)/;
								my $out = $results_dir . $name . ".fa";
								open(FINAL, ">>$out");
								print FINAL ">", $lineage, "\n", $seq1{$c}, "\n";
								print FINAL ">", $compare2{$lineage}, "\n", $sub, "\n";
								close(FINAL);
								}
							}				
						}	
					}
				}
			}
		my $cleanup = system("rm exonHomolog.fa exonHomolog2.fa");
		}
	}

sub divpoly {
	my %div;
	open(IN,"<$divpoly");
	while(<IN>) {
		my @d = split(/\t/,$_);
		if ($d[1] eq 'UTR') {
			$div{$d[0]}{'net'} = $d[5];
			$div{$d[0]}{'raw'} = $d[4];
			}
		}
	close(IN);
	return(\%div);
	}
	
sub calculateDiv {
	my ($div) = @_;
	my %final;
	my %div = %{$div}; 
	my @aln = <$results_dir*aln>;
	
	foreach my $aln (@aln) {
		open(IN, "<$aln");
		my %species; my $id;
				
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+)/) {
				$id = $1;
				}
			else {
				$species{$id} .= $line;
				}
			}
			
		foreach my $id (keys %species) {
			unless($id =~ m/$genus/) {
				delete($species{$id});
				}
			}	
				
		if (scalar(keys %species) > 0) {
			my $div = 'NA';			
			foreach my $c1 (keys %compare2) {
				if ($species{$c1} && $species{$compare2{$c1}}) {
					$div = div($species{$c1},$species{$compare2{$c1}});		
					unless ($div eq 'NA') {
						$div = sprintf("%.4f",$div);					
						}
					}
				}	
							
			my %info;
			foreach my $c (keys %species) {
				my $seq = $species{$c};
				my $l = ($seq =~ s/[atgc]//ig);
				$seq = $species{$c};
				my $gc = ($seq =~ s/[gc]//ig);
				$gc = sprintf("%.3f",$gc/$l);
				$info{$c} = {'length' => $l, 'gc' => $gc};
				}
			
			my @c = sort {$info{$a}{'length'} <=> $info{$b}{'length'}} keys %info;
			
			my $gc = $info{$c[$#c]}{'gc'};
			my $maxLength = $info{$c[$#c]}{'length'};	
			my $name = $1 if $aln =~ m/(ENSA\S+)\.fa/;
		
			if ($maxLength > $exonLength) {
				if ($gc > $minGC && $gc < $maxGC) {
					my $short = $name;
					if ($name =~ m/(\S+)_/) {
						$short = $1;
						}
					$final{$name} = {'length' => $maxLength, 'gc' => $gc, 'div' => $div};
					if ($div{$short}) {
						$final{$name}{'rawdiv'} = $div{$short}{'raw'};
						$final{$name}{'netdiv'} = $div{$short}{'net'};
						}
					else {
						$final{$name}{'rawdiv'} = 'NA';
						$final{$name}{'netdiv'} = 'NA';
						}
					}
				}	
			}
		}	
	close(OUT);
	return(\%final);
	}
		
sub div {	
	my ($seq1,$seq2) = @_;
	my @b1 = split(//,$seq1);
	my @b2 = split(//,$seq2);
	
	my %ti = ('A' => 'G', 'G' => 'A', 'C' => 'T', 'T' => 'C');
	
	my $gc; my $ti = 0; my $tv = 0; my $l = 0;

	for (my $j = 0; $j < scalar(@b1); $j++) {
		if (uc($b1[$j]) =~ m/[A|T|G|C]/i) {
			if (uc($b1[$j]) =~ m/[G|C]/i) {
				$gc++;
				}	
			if (uc($b2[$j]) =~ m/[A|T|G|C]/) {	
				$l++;
				unless (uc($b1[$j]) eq uc($b2[$j])) {
					#this is a mutation!
					if ( $ti{uc($b1[$j])} eq uc($b2[$j]) ) {
						$ti++;
						}
					else {	
						$tv++;
						}
					}	
				}	
			}
		}
				
	my $div;	
	if ($l > 0) {	
		my $p = $ti/$l; my $q = $tv/$l; my $w = $gc/$l;
				
		if ($w == 1) {
			$div = 'NA';
			}
		else {	
			my $a = 1 - ( $p/ (2 * $w * (1 - $w) ) ) - $q;
			my $b = 1 - 2 * $q ;
			if ($a <= 0 || $b <= 0 ) {
				$div = 'NA';
				}
			else {	
				$div = (-2*$w) * ( 1 - $w) *  log( $a ) - 0.5 * ( 1 - 2 * $w * ( 1 - $w ) ) * log($b);
				}
			}	
		}
	else {
		$div = 'NA';
		}
	return($div);	
	}	

			