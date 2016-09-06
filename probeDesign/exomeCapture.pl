use warnings;
use strict;
use Bio::DB::Fasta;

my $dir = '/Users/singhal/Desktop/activeWork/exomeCapture/';
my $results_dir = $dir . 'results/';
my $seq1 = $dir . 'seqFiles/Carlia_N_trinity.fa.final.annotated';
my $seq2 = $dir . 'seqFiles/Carlia_S_trinity.fa.final.annotated';
my $seq3 = $dir . 'seqFiles/Lampro_N_trinity.fa.final.annotated';
my $seq4 = $dir . 'seqFiles/Lampro_C_trinity.fa.final.annotated';
my $seq5 = $dir . 'seqFiles/Lampro_S_trinity.fa.final.annotated';
my $seq6 = $dir . 'seqFiles/Sapro_C_trinity.fa.final.annotated';
my $seq7 = $dir . 'seqFiles/Sapro_S_trinity.fa.final.annotated';
my $dbP = $dir . "genomes/anolis_Prot.fa";
my $dbG = $dir . "genomes/anolis_DNA.fa";
#consider all exons of this length or greater
my $exonLength = 200;
#consider all exons which are matched at least this much or greater
my $exonMatch = 0.25; 
my $anolisExons = $dir . 'anolisExons.fa';
my $exonFile = $dir . 'uniqueExons.fa';
my $np = 4;
my $eval = 1e-10; #evalue required for matches
my $divout = $dir . 'divergence.out';
my %compare1 = ('Sapro' => 'Lampro',
				'Lampro' => 'Carlia',
				'Carlia' => 'Sapro');
my %compare2 = ('Sapro_C' => 'Sapro_S',
				'Carlia_N' => 'Carlia_S',
				'Lampro_N' => 'Lampro_C',
				'Lampro_C' => 'Lampro_S');
my $minGC = 0.3;
my $maxGC = 0.8;

##########################
# start the subroutines  #
##########################

my $anolisOut = $anolisExons . ".unique";
mkdir($results_dir) unless(-d $results_dir);
#findUniqueExons($anolisExons,$dbP,$dbG);
#alignment($seq1, 'Carlia_N');
#alignment($seq2, 'Carlia_S');
#alignment($seq3, 'Lampro_N');
#alignment($seq4, 'Lampro_C');
#alignment($seq5, 'Lampro_S');
#alignment($seq6, 'Sapro_C');
#alignment($seq7, 'Sapro_S');
#alignHomologs();
calculateDiv();

##########################
# behold the subroutines #
##########################

sub calculateDiv {
	my @aln = <$results_dir*aln>;
	open(OUT, ">$divout");
	
	foreach my $aln (@aln) {
		open(IN, "<$aln");
		my %species; my %genus; my $id; my $genus; 
				
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+)/) {
				$id = $1;
				$genus = $1 if $id =~ m/(\S+)_/;
				}
			else {
				$species{$id} .= $line;
				$genus{$genus} .= $line;
				}
			}
			
		my @div;	
		my $zero = 0;
		foreach my $c1 (keys %compare1) {
			my $div = div($genus{$c1},$genus{$compare1{$c1}});
			$div = sprintf("%.4f",$div) unless ($div eq 'NA');
			unless ($div eq 'NA') {
				$zero++ if $div == 0;
				}
			push(@div,$div);
			}
			
			
		foreach my $c1 (keys %compare2) {
			my $div = 'NA';
			if ($species{$c1} && $species{$compare2{$c1}}) {
				$div = div($species{$c1},$species{$compare2{$c1}});		
				unless ($div eq 'NA') {
					$div = sprintf("%.4f",$div);					
					}
				}
			push(@div,$div);	
			}	
			
		my @length; my @gc; my $length = 0;
		foreach my $c (keys %species) {
			my $seq = $species{$c};
			my $l = ($seq =~ s/[atgc]//ig);
			$seq = $species{$c};
			my $gc = ($seq =~ s/[gc]//ig);
			$length++ if $l < $exonLength;
			push(@length,$l);
			push(@gc,$gc);
			}
		my $avgLength = sprintf("%d",average(\@length));	
		my $avgDiv = sprintf("%.4f",average(\@div));
		my $avgGC = sprintf("%.4f",average(\@gc)/$avgLength);
		my $var = sprintf("%.4f",variance($avgDiv,\@div));
		my $name = $1 if $aln =~ m/(ENSA.*exon\d+)/;
		
		unless ($length) {
			unless ($zero) {
				unless ($var > 0.002) {
					if ($avgLength > $exonLength) {
						if ($avgGC > $minGC && $avgGC < $maxGC) {
							print OUT $name, "\t", $avgLength, "\t", $avgDiv, "\t", $var, "\t", $avgGC, "\t";
							print OUT join("\t",@div);
							print OUT "\n";
							}
						}
					}
				}
			}	
		}
	close(OUT);	
	}
	
sub variance {
	my ($avg, $a) = @_;
	my @a = @{$a};
	my $var;
	if (scalar(@a) > 0) {
		my $sum = 0; my $num;
		foreach my $var (@a) {
			unless ($var eq 'NA') {
				$sum += ($avg - $var) ** 2;
				$num++;
				}
			}	
		$var = $sum/$num;
		}
	else {
		$var = 0;
		}
	return($var);
	}		
	
sub average {
	my ($a) = @_;
	my @a = @{$a};
	my $avg;
	if (scalar(@a) > 0) {
		my $sum = 0; my $num;
		foreach my $var (@a) {
			unless ($var eq 'NA') {
				$sum += $var;
				$num++;
				}
			}	
		$avg = $sum/$num;
		}
	else {
		$avg = 0;
		}
	return($avg);
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

sub alignHomologs {
	my @seq = <$results_dir*ENS*fa>;
	
	foreach my $file (@seq) {
		open(IN, "<$file");
		my %g; my %s;
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+)/) {
				my $genus = $1 if $line =~ m/>(\S+)_/;
				my $species = $1 if $line =~ m/>(\S+)/;
			
				$g{$genus}++; $s{$species}++;
				}
			}
		close(IN);
		if (scalar(keys %g) > 2) {
			my $out = $file . ".aln";
			my $call1 = system("muscle -in $file -out $out");
			my $call2 = system("rm $file") if (-f $out);
			}
		}
	}
	
sub alignment {
	my ($seq,$id) = @_;
		
	my %seq;
	open(SEQ, "<$seq");
	my $seqout = $seq . ".blast";
	open(SEQOUT, ">$seqout");
	while(<SEQ>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $id = $1;
			chomp(my $seq = <SEQ>);
			$seq{$id} = $seq;
			print SEQOUT ">", $id, "\n$seq\n";
			}
		}
	close(SEQ); close(SEQOUT);
	
	my %exons;
	my $anolisBlast = $anolisOut . ".blast";
	open(ANOLEIN, "<$anolisOut");
	open(ANOLEOUT, ">$anolisBlast");
	while(<ANOLEIN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $id = $1;
			chomp(my $seq = <ANOLEIN>);
			$exons{$id} = $seq;
			print ANOLEOUT ">", $id, "\n$seq\n";
			}
		}
	close(ANOLEOUT); close(ANOLEIN);
	
	my $call1 = system("formatdb -i $seqout -p F") unless(-f $seqout . ".nin");
	my $call2 = system("formatdb -i $anolisBlast -p F") unless(-f $anolisBlast . ".nin");
	my $out1 = $seqout . ".blast.out";
	my $out2 = $seqout . ".recipblast.out";
	my $call3 = system("blastall -p blastn -d $seqout -i $anolisBlast -m 8 -a $np -b 1 -o $out1 -e $eval") unless(-f $out1);
	my $call4 = system("blastall -p blastn -d $anolisBlast -i $seqout -m 8 -b 1 -o $out2 -a $np -e $eval") unless (-f $out2);
		
	my %matches;
	open(IN, "<$out1");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		$matches{$d[0]}{$d[1]} = 'no';
		}
	close(IN);
		
	open(IN, "<$out2");
	while(<IN>) {
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		if ($matches{$d[1]}{$d[0]}) {
			$matches{$d[1]}{$d[0]} = 'yes';
			}
		}
	close(IN);
	
	foreach my $exon (keys %exons) {
		if ($matches{$exon}) {
			foreach my $match (keys %{$matches{$exon}}) {
				if ($matches{$exon}{$match} eq 'yes') {
					my $out1 = 'exonHomolog.fa';
					open(OUT1, ">$out1");
					print OUT1 ">", $exon, "\n" , $exons{$exon}, "\n";
					my $elength = length($exons{$exon});
					close(OUT1);
		
					my $out2 = 'exonHomolog2.fa';
					open(OUT2, ">$out2");
					print OUT2 ">", $match, "\n" , $seq{$match}, "\n";
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
								$sub = substr $seq{$match}, $start, $length;
								}
							elsif ($d[4] =~ m/\+/ && $d[8] =~ m/\-/) {
								$start = $d[7];
								$sub = substr $seq{$match}, $start, $length;
								$sub = reverse($sub);
								$sub =~ tr/ATGCatgc/TACGtacg/;
								}
							elsif ($d[4] =~ m/\-/ && $d[8] =~ m/\-/) {
								$start = $d[7];
								$sub = substr $seq{$match}, $start, $length;
								}
							elsif ($d[4] =~ m/\-/ && $d[8] =~ m/\+/) {
								$start = $d[6];
								$sub = substr $seq{$match}, $start, $length;
								$sub = reverse($sub);
								$sub =~ tr/ATGCatgc/TACGtacg/;
								}
				
							if ($length/$elength >= $exonMatch) {
								my $out = $results_dir . $exon . ".fa";
								open(FINAL, ">>$out");
								my $end = $start + $length - 1;
								my $loc = $match . '_' . $start . "_" . $end;
								print FINAL ">$id\t$loc\n", $sub, "\n";
								close(FINAL);
								}
							}				
						}	
					}
				}
			}
		}		
	my $cleanup = system("rm exonHomolog.fa exonHomolog2.fa");
	}


sub findUniqueExons {
	my ($anolisExons,$databaseP,$databaseG) = @_;
		
	open(OUT, ">$anolisExons");	
		
	my $anno = makeHash($dbP); 	
	#need to use exonerate to define utr etc; call to external program
	my $dbP = Bio::DB::Fasta->new($databaseP);
	my $dbG = Bio::DB::Fasta->new($databaseG);

	my %anno = %$anno;
	
	my %seq;
	open(IN, "<$databaseP");
	while(<IN>) {
		if ($_ =~ m/(ENSACAP\S+).*(ENSACAG\S+)/) {
			my $gene = $1;
			my $prot = $1;
			chomp(my $seq = <IN>);
			if ($seq{$gene}) {
				if ($seq{$gene}{'seq'} < length($seq)) {
					$seq{$gene}{'prot'} = $prot;
					$seq{$gene}{'seq'} = $seq;
					}
				}
			else {	
				$seq{$gene}{'prot'} = $prot;
				$seq{$gene}{'seq'} = $seq;
				}
			}
		}
	close(IN);
				
	foreach my $geneID (keys %seq) {			
		my %final;	
		my $anolisID = $seq{$geneID}{'prot'};
		if ($dbG->get_Seq_by_id($anno{$anolisID}{'contig'})) {
			#need to run exonerate			
			my $DNA = $dbG->get_Seq_by_id($anno{$anolisID}{'contig'})->subseq($anno{$anolisID}{'start'} => $anno{$anolisID}{'end'});
			my $seq = $dbP->get_Seq_by_id($anolisID)->seq;	

			my $target2 = "target2.fa";
			my $query2 = "query2.fa";
		    	
			open(T2, ">$target2");
			open(Q2, ">$query2");
			print T2 ">$anolisID\n$DNA\n";
			print Q2 ">protein\n$seq\n";
			my @call = `exonerate --model protein2genome $query2 $target2 --showvulgar no --showalignment no --showtargetgff yes`;
			my %cds;
			foreach(@call){
				if ($_ =~ m/orientation/) {
				   	$cds{'orient'} = $1 if $_ =~ m/([\+|\-])\s+\./;
					}
				elsif ($_ =~ m/exon\s+/) {
					my $bounds = $1 if $_ =~ m/exon\t([0-9]+\t[0-9]+)/;;
					my @bounds = split("\t", $bounds);
	    			$cds{'exon'}{$bounds[0]}=$bounds[1];
					}
				}
				
			#need to allow for a failure to find cds 
			if (keys %cds){
				my $tracker = 1;
				if ($cds{'orient'} eq '+') {
					foreach(sort {$a <=> $b} keys %{$cds{'exon'}}){
				  	  my $start = $_;
	   					my $end = $cds{exon}{$start};
	   					my $length = $end - $start + 1;
	   					$start = $start - 1;
	   					my $sub = substr $DNA, $start, $length;
			    		my $loc = $anno{$anolisID}{'contig'};
			    		my $loc_s = $anno{$anolisID}{'start'} + $start;
			    		my $loc_e = $loc_s + $length - 1;
			    		$loc = $loc . '_' . $loc_s . '_' . $loc_e;

	   					my $id = ">" . $anolisID . "_exon" . $tracker;
						$final{$id}{'seq'} = $sub;
						$final{$id}{'loc'} = $loc;
    					$tracker++;
						}
						}
				else {
					foreach(sort {$b <=> $a} keys %{$cds{'exon'}}){
					    my $start = $_;
	    				my $end = $cds{exon}{$start};
	    				my $length = $end - $start + 1;
	    				$start = $start - 1;
	    				my $sub = substr $DNA, $start, $length;
	    				$sub = reverse($sub);
	    				$sub =~ tr/ATGCatgc/TACGtacg/;

				    	my $loc = $anno{$anolisID}{'contig'};
				    	my $loc_s =$anno{$anolisID}{'start'} + $start;
						my $loc_e = $loc_s + $length - 1;                       
						$loc = $loc. '_' .$loc_s . '_' . $loc_e;
	
						my $id = ">" . $anolisID . "_exon" . $tracker;
						$final{$id}{'seq'} = $sub;
						$final{$id}{'loc'} = $loc;
						}
					}
				foreach my $id (keys %final) {
					print OUT $id, "\t$final{$id}{'loc'}\n", $final{$id}{'seq'}, "\n";
					}
				}	
			else {
				print "couldn't define CDS for $anolisID\n";
				}
			}
		}
	close(OUT);	
	my $call1 = system("cd-hit-est -i $anolisExons -o $anolisOut -c 1.00 -l $exonLength");
	my $call2 = system("rm $anolisOut" . ".*");
	}	
		
sub makeHash {
	my ($database) = @_;
		
	my %anno;
	open(ANNO, "<$database");
	while(<ANNO>){
		chomp(my $line = $_);
		if ($line =~ m/^>(\S+)/) {
			my $id = $1;
			my @a = split(/\s/,$line);
			my @d = split(/:/,$a[2]);
			my $contig = $d[2];
			my $start = $d[3];
			my $end = $d[4];
			$anno{$id} = {'contig'=>$contig, 'start'=>$start, 'end'=>$end};
			}
		}	
	close(ANNO);

	return(\%anno)
	}			