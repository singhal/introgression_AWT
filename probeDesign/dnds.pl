use strict;
use warnings;

my @contacts = (["Carlia_N","Carlia_S"], ["Sapro_C","Sapro_S"], ["Lampro_N", "Lampro_C", "Lampro_S"]);

my @seqfiles = 
</Users/singhal/Desktop/activeWork/speciesTree/seqfiles/*annotated>;
my @blastfiles = 
</Users/singhal/Desktop/activeWork/speciesTree/seqfiles/*blast>;
my %seq;
my $resultsDir = '/Users/singhal/Desktop/activeWork/dnds/';
my $seqDir = $resultsDir . 'align/';	
my $seqfile = $resultsDir . 'referenceCodingSeq.fa';
my $dndsout = $resultsDir . 'dnds.out';
my $np = 4;
my $chi = 5.99;
my $perCut = 0.4;
my $lengthCut = 300;

mkdir($seqDir) unless(-d $seqDir);
mkdir($resultsDir) unless(-d $resultsDir);
makeSeq();
recipBlast();
alignHomologs();
calculateDnDsALL();
foreach my $c (@contacts) {
	calculateDnDsPaired($c);
	}
foreach my $c (@contacts) {
	calculateDnDsBranch($c);
	}	
	
sub calculateDnDsBranch {
	my ($c) = @_;
	
	my @c = @{$c};
	my %c;
	$c{$_}++ foreach @c;
	my $contact = join("_",@c);
	my $out = $resultsDir . 'dndsBranch_' . $contact . '.out';
	
	my @aln = <$seqDir*aln>;

	open(DNDSOUT, ">$out");
	
	foreach my $aln (@aln) {
		open(IN, "<$aln");
		my %seq; my $id; my @seq; 
				
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+)/) {
				$id = $1;
				}
			else {
				$seq{$id} .= $line;
				}
			}
			
		my $length1 = length($seq{$id});
		foreach my $id (keys %seq) {
			my $seq = $seq{$id};
			$seq =~ s/\-//g;
			unless (length($seq)/$length1 >= $perCut || length($seq) >= $lengthCut) {
				delete($seq{$id});
				}
			}
		
		#need to figure out the frame
		my $length;
		my @stops = ();
		for (my $i = 0; $i < 3; $i++) {
			my $stop = 0;
			foreach my $id (keys %seq) {
				my $subseq = substr $seq{$id}, $i;
				$length = length($subseq);
				my $aa = translate($subseq);
				$stop++ if $aa =~ m/[A-Z]\*[A-Z]/;
				}
			push(@stops,$i) if $stop == 0;
			}
			
		my $name = $1 if $aln =~ m/(ENSA.*)\.fa/;		
		if (scalar(@stops) == 0) {
			print DNDSOUT $name, "\t", $length, "\t", "NA\n";
			}
		elsif (scalar(@stops) > 1) {
			my @omega; 
			my @cdslength;
			foreach my $stop (@stops) {
				my ($omega,$cdslength) = runPaml_branch(\%c,\%seq,$stop);
				push(@omega,$omega) if $omega;
				push(@cdslength,$cdslength) if $cdslength; 
				}
			@omega = sort {$a <=> $b} @omega;
			print DNDSOUT $name, "\t", $cdslength[0], "\t", $omega[0], "\n";	
			}
		else {	
			my ($omega,$cdslength) = runPaml_branch(\%c,\%seq,$stops[0]);
			print DNDSOUT $name, "\t", $cdslength, "\t", $omega, "\n";	
			}
		}
	close(DNDSOUT);
	}	

sub runPaml_branch {
	my ($c,$seq,$stop) = @_;
	my $out = $resultsDir . "eugong_branch.phy";
	
	my %c = %{$c};
	my %seq = %{$seq};
	my @id = keys %seq;	
	my $cdslength;
	my $tree = $resultsDir . "treefile_branch";
	open(TREE, ">$tree");
	print TREE "(";
	for (my $i = 0; $i < scalar(@id); $i++) {
		my $id = $id[$i];
		$cdslength = length(substr $seq{$id}, $stop);
		$cdslength = int($cdslength/3) * 3;
		$seq{$id} = substr $seq{$id}, $stop, $cdslength;
		if ($c{$id}) {
			print TREE $id, " #1";
			}
		else {
			print TREE $id;
			}
		if ($i == (scalar(@id) - 1) ) {
			print TREE ");";
			}
		else {
			print TREE ",";
			}
		}
	close(TREE);	
						
	open(OUT, ">$out");
	my $num = scalar(@id);
	print OUT " $num $cdslength\n";	
	foreach my $id (keys %seq) {
		print OUT "$id  $seq{$id}\n";
		}
	close(OUT);
								
	my @lnl;		
	my $call = system("codeml codeml_branch.ctl");
	open(IN, "<eugong_branch.out");
	my $omega = 'NA';
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/^site class/) {
			chomp(my $line1 = <IN>);
			chomp(my $line2 = <IN>);
			chomp(my $line3 = <IN>);
			
			my @d1 = split(/\s+/,$line1);
			my @d2 = split(/\s+/,$line2);
			my @d3 = split(/\s+/,$line3);
			
			if ($d1[3] > 0 && $d1[4] > 0) {
				$omega = $d3[5];
				}
			else {
				$omega = $d2[4];
				}
			$omega = 'NA' if $omega =~ m/^999\.0+$/;			
			}
		}
	close(IN);	
	return($omega,$cdslength);	
	}	
	
sub calculateDnDsPaired {
	my ($c) = @_;
	
	my @c = @{$c};
	my %c;
	$c{$_}++ foreach @c;
	my $contact = join("_",@c);
	my $out = $resultsDir . 'dnds_' . $contact . '.out';
	
	my @aln = <$seqDir*aln>;

	open(DNDSOUT, ">$out");
	
	foreach my $aln (@aln) {
		open(IN, "<$aln");
		my %seq; my $id; my @seq; 
				
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+)/) {
				$id = $1;
				}
			else {
				$seq{$id} .= $line;
				}
			}
			
		my $length1 = length($seq{$id});
		foreach my $id (keys %seq) {
			if ($c{$id}) {	
				my $seq = $seq{$id};
				$seq =~ s/\-//g;
				unless (length($seq)/$length1 >= $perCut || length($seq) >= $lengthCut) {
					delete($seq{$id});
					}
				}
			else {
				delete($seq{$id});
				}
			}
		
		#need to figure out the frame
		my $length;
		my @stops = ();
		for (my $i = 0; $i < 3; $i++) {
			my $stop = 0;
			foreach my $id (keys %seq) {
				my $subseq = substr $seq{$id}, $i;
				$length = length($subseq);
				my $aa = translate($subseq);
				$stop++ if $aa =~ m/[A-Z]\*[A-Z]/;
				}
			push(@stops,$i) if $stop == 0;
			}
			
		my $name = $1 if $aln =~ m/(ENSA.*)\.fa/;		
		if (scalar(@stops) == 0) {
			print DNDSOUT $name, "\t", $length, "\t", "NA\n";
			}
		elsif (scalar(@stops) > 1) {
			my @omega; 
			my @cdslength;
			foreach my $stop (@stops) {
				my ($omega,$cdslength) = runPaml_pair(\%seq,$stop);
				push(@omega,$omega) if $omega;
				push(@cdslength,$cdslength) if $cdslength; 
				}
			@omega = sort {$a <=> $b} @omega;
			print DNDSOUT $name, "\t", $cdslength[0], "\t", $omega[0], "\n";	
			}
		else {	
			my ($omega,$cdslength) = runPaml_pair(\%seq,$stops[0]);
			print DNDSOUT $name, "\t", $cdslength, "\t", $omega, "\n";	
			}
		}
	close(DNDSOUT);
	}	

sub runPaml_pair {
	my ($seq,$stop) = @_;
	my $out = $resultsDir . "eugong_pair.phy";
							
	my %seq = %{$seq};
	my @id = keys %seq;	
	my $cdslength;
	my $tree = $resultsDir . "treefile_pair";
	open(TREE, ">$tree");
	print TREE "(";
	for (my $i = 0; $i < scalar(@id); $i++) {
		my $id = $id[$i];
		$cdslength = length(substr $seq{$id}, $stop);
		$cdslength = int($cdslength/3) * 3;
		$seq{$id} = substr $seq{$id}, $stop, $cdslength;
		print TREE $id;
		if ($i == (scalar(@id) - 1) ) {
			print TREE ");";
			}
		else {
			print TREE ",";
			}
		}
	close(TREE);	
						
	open(OUT, ">$out");
	my $num = scalar(@id);
	print OUT " $num $cdslength\n";	
	foreach my $id (keys %seq) {
		print OUT "$id  $seq{$id}\n";
		}
	close(OUT);
			
	my @lnl;		
	my $call = system("codeml codeml_pair.ctl");
	open(IN, "<eugong_pair.out");
	my $omega = 'NA';
	while(<IN>) {
		if ($_ =~ m/omega/) {
			$omega = $1 if $_ =~ m/([0-9|\.]+)/;
			$omega = 'NA' if $omega =~ m/^999\.0+$/;			
			}
		}
	close(IN);	
	return($omega,$cdslength);	
	}

sub calculateDnDsALL {
	my @aln = <$seqDir*aln>;
	open(DNDSOUT, ">$dndsout");
	
	foreach my $aln (@aln) {
		open(IN, "<$aln");
		my %seq; my $id; my @seq; 
				
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+)/) {
				$id = $1;
				}
			else {
				$seq{$id} .= $line;
				}
			}
			
		my $length1 = length($seq{$id});
		foreach my $id (keys %seq) {
			my $seq = $seq{$id};
			$seq =~ s/\-//g;
			unless (length($seq)/$length1 >= $perCut || length($seq) >= $lengthCut) {
				delete($seq{$id});
				}
			}
		
		#need to figure out the frame
		my $length;
		my @stops = ();
		for (my $i = 0; $i < 3; $i++) {
			my $stop = 0;
			foreach my $id (keys %seq) {
				my $subseq = substr $seq{$id}, $i;
				$length = length($subseq);
				my $aa = translate($subseq);
				$stop++ if $aa =~ m/[A-Z]\*[A-Z]/;
				}
			push(@stops,$i) if $stop == 0;
			}
			
		my $name = $1 if $aln =~ m/(ENSA.*)\.fa/;		
		if (scalar(@stops) == 0) {
			print DNDSOUT $name, "\t", $length, "\t", "NA", "\tNA\n";
			}
		elsif (scalar(@stops) > 1) {
			my @omega1; my @omega2; 
			my @cdslength;
			foreach my $stop (@stops) {
				my ($omega1,$omega2,$cdslength) = runPaml(\%seq,$stop);
				push(@omega1,$omega1) if $omega1;
				push(@omega2,$omega2) if $omega2;
				push(@cdslength,$cdslength) if $cdslength; 
				}
			@omega1 = sort {$a <=> $b} @omega1;
			@omega2 = sort {$a <=> $b} @omega2;
			print DNDSOUT $name, "\t", $cdslength[0], "\t", $omega1[0], "\t", $omega2[0], "\n";	
			}
		else {	
			my ($omega1,$omega2,$cdslength) = runPaml(\%seq,$stops[0]);
			print DNDSOUT $name, "\t", $cdslength, "\t", $omega1, "\t", $omega2, "\n";	
			}
		}
	close(DNDSOUT);
	}
	
sub runPaml {
	my ($seq,$stop) = @_;
	my $out = $resultsDir . "eugong.phy";
							
	my %seq = %{$seq};
	my @id = keys %seq;	
	my $cdslength;
	my $tree = $resultsDir . "treefile";
	open(TREE, ">$tree");
	print TREE "(";
	for (my $i = 0; $i < scalar(@id); $i++) {
		my $id = $id[$i];
		$cdslength = length(substr $seq{$id}, $stop);
		$cdslength = int($cdslength/3) * 3;
		$seq{$id} = substr $seq{$id}, $stop, $cdslength;
		print TREE $id;
		if ($i == (scalar(@id) - 1) ) {
			print TREE ");";
			}
		else {
			print TREE ",";
			}
		}
	close(TREE);	
						
	open(OUT, ">$out");
	my $num = scalar(@id);
	print OUT " $num $cdslength\n";	
	foreach my $id (keys %seq) {
		print OUT "$id  $seq{$id}\n";
		}
	close(OUT);
			
	my @lnl;		
	my $call = system("codeml codeml.ctl");
	open(IN, "<eugong.out");
	my $omega1; my $omega2;
	while(<IN>) {
		if ($_ =~ m/omega/) {
			$omega1 = $1 if $_ =~ m/([0-9|\.]+)/;
			$omega1 = 'NA' if $omega1 =~ m/^999\.0+$/;			
			}
		elsif ($_ =~ m/lnL/) {
			push(@lnl, $1) if $_ =~ m/\):\s+([0-9|\.|\-]+)/;
			}
		elsif ($_ =~ m/dN\/dS for site classes \(K\=3\)/) {
			if ( (2 * abs($lnl[2] - $lnl[1])) >= $chi) {	
				my $junk = <IN>; my $junk2 = <IN>; my $data = <IN>;
				my @d = split(/\s+/,$data);
				$omega2 = $d[3];
				}
			else {
				$omega2 = 'NA';
				}
			}
		}
	close(IN);	
	return($omega1,$omega2,$cdslength);	
	}

sub translate {
	my $string = shift;
	$string = uc($string);
	my @codons = $string =~ m/(\S\S\S)/g;
	my %codons = (	'ATG'=>'M','ACG'=>'T','CTG'=>'L','CCG'=>'P','GTG'=>'V','GCG'=>'A','TTG'=>'L','TCG'=>'S',
					'ATA'=>'I','ACA'=>'T','CTA'=>'L','CCA'=>'P','GTA'=>'V','GCA'=>'A','TTA'=>'L','TCA'=>'S',
					'ATC'=>'I','ACC'=>'T','CTC'=>'L','CCC'=>'P','GTC'=>'V','GCC'=>'A','TTC'=>'F','TCC'=>'S',
					'ATT'=>'I','ACT'=>'T','CTT'=>'L','CCT'=>'P','GTT'=>'V','GCT'=>'A','TTT'=>'F','TCT'=>'S',
					'AGG'=>'R','AAG'=>'K','CGG'=>'R','CAG'=>'Q','GGG'=>'G','GAG'=>'E','TGG'=>'W','TAG'=>'*',
					'AGA'=>'R','AAA'=>'K','CGA'=>'R','CAA'=>'Q','GGA'=>'G','GAA'=>'E','TGA'=>'*','TAA'=>'*',
					'AGC'=>'S','AAC'=>'N','CGC'=>'R','CAC'=>'H','GGC'=>'G','GAC'=>'D','TGC'=>'C','TAC'=>'Y',
					'AGT'=>'S','AAT'=>'N','CGT'=>'R','CAT'=>'H','GGT'=>'G','GAT'=>'D','TGT'=>'C','TAT'=>'Y');
	my $translate;
	foreach(@codons) {
		if ($codons{$_}) {
			$translate = $translate . $codons{$_};
			}
		else {
#			print "ERROR: ILLEGAL PASS TO CODON TRANSLATION: $_ is not a codon!\n";
			$translate = $translate . 'X';
			}
		}
	return($translate);
	}	


sub alignHomologs {
	my @seq = <$seqDir*ENS*fa>;
	foreach my $seq (@seq) {
		open(IN, "<$seq");
		my $out = $seq . ".aln";
		my $call1 = system("muscle -in $seq -out $out") unless (-f $out);
		my $call2 = system("rm $seq") if (-f $out);
		}
	}
	
sub recipBlast {
	my $call1 = system("formatdb -i $seqfile -p F") unless (-f $seqfile . ".nin");
	
	my %anolisSeq;
	open(IN, "<$seqfile");
	while(<IN>) {
		chomp(my $line = $_);
		if ($line =~ m/>(\S+)/) {
			my $id = $1;
			chomp(my $seq = <IN>);
			$anolisSeq{$id} = $seq;
			}
		}
	close(IN);	
		
	foreach my $seq (@blastfiles) {
		my $call2 = system("formatdb -i $seq -p F") unless (-f $seq . ".nin");
		my $species = $1 if $seq =~ m/([a-z]+_[a-z])_trinity/i;
		my $out1 = $resultsDir . $species . "toReference.blast.out";
		my $out2 = $resultsDir . "referenceTo" . $species . ".blast.out";
		my $call3 = system("blastall -p blastn -d $seqfile -i $seq -m 8 -b 1 -o $out1 -a $np -e 1e-100") unless (-f $out1);
		my $call4 = system("blastall -p blastn -d $seq -i $seqfile -m 8 -b 1 -o $out2 -a $np -e 1e-100") unless (-f $out2);
		
		my %contig;
		open(IN, "<$seq");
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+)/) {
				my $id = $1;
				chomp(my $seq = <IN>);
				$contig{$id} = $seq;
				}
			}
		close(IN);	
		
		my %matches;
		open(IN, "<$out2");
		while(<IN>) {
			chomp(my $line = $_);
			my @d = split(/\t/,$line);
			$matches{$d[0]}{$d[1]} = 'no';
			}
		close(IN);
		
		open(IN, "<$out1");
		while(<IN>) {
			chomp(my $line = $_);
			my @d = split(/\t/,$line);
			if ($matches{$d[1]}{$d[0]}) {
				$matches{$d[1]}{$d[0]} = 'yes';
				}
			}
		close(IN);
				
		foreach my $anoleID (keys %matches) {
			my $outfile = $seqDir . $anoleID . '.fa';
			open(HOM1, ">homolog1.fa");
			print HOM1 ">$anoleID\n$anolisSeq{$anoleID}\n";
			close(HOM1);
			foreach my $gene (keys %{$matches{$anoleID}}) {
				if ($matches{$anoleID}{$gene} eq 'yes') {
					#exonerate time!
					open(HOM2, ">homolog2.fa");
					print HOM2 ">$gene\n$contig{$gene}\n";
					close(HOM2);
					
					my @call = `exonerate -m coding2coding homolog1.fa homolog2.fa --showvulgar yes --showalignment no --showtargetgff no`;
		
					if ($call[3]) {
						if ($call[3] =~ m/vulgar/) {
							my @d = split(/\s/,$call[3]);
							my $length = abs($d[6] - $d[7]);					
							my $sub; 
			
							if ($d[4] =~ m/\+/ && $d[8] =~ m/\+/) {
								my $start = $d[6];
								$sub = substr $contig{$gene}, $start, $length;
								}
							elsif ($d[4] =~ m/\+/ && $d[8] =~ m/\-/) {
								my $start = $d[7];
								$sub = substr $contig{$gene}, $start, $length;
								$sub = reverse($sub);
								$sub =~ tr/ATGCatgc/TACGtacg/;
								}
							elsif ($d[4] =~ m/\-/ && $d[8] =~ m/\-/) {
								my $start = $d[7];
								$sub = substr $contig{$gene}, $start, $length;
								}
							elsif ($d[4] =~ m/\-/ && $d[8] =~ m/\+/) {
								my $start = $d[6];
								$sub = substr $contig{$gene}, $start, $length;
								$sub = reverse($sub);
								$sub =~ tr/ATGCatgc/TACGtacg/;
								}
								
							open(OUT, ">>$outfile");	
							print OUT ">", $species, "\n", $sub, "\n";
							close(OUT);
							}										
						}									
					}
				}	
			}
		}
	my $call5 = system("rm homolog1.fa homolog2.fa");	
	}

sub makeSeq {
	foreach my $seq (@seqfiles) {
		open(IN, "<$seq");
		while(<IN>) {
			chomp(my $line = $_);
			if ($line =~ m/>(\S+).*ENS.*/) {
				my $id = $1;
				my @d = split(/\s+/,$line);
				chomp(my $s = <IN>);
			
				my $complete = 'incl';			
				my $start = $1 if $d[1] =~ m/gs(\d+)/;
				my $end = $1 if $d[1] =~ m/ge(\d+)/;
				my $length = $end - $start + 1;
			
				if ($d[1] =~ m/5u/ && $d[1] =~ m/3u/) {
					$complete = 'compl';
					}
				
				$s = substr $s, $start - 1, $length;
				
				if ($seq{$d[2]}) {
					unless ($seq{$d[2]}{'eval'} > $d[-1]) {
						if ($seq{$d[2]}{'full'} eq 'incl' && $complete eq 'compl') {
							$seq{$d[2]} = {'full' => $complete, 'eval' => $d[-1], 'length' => $length, 'seq' => $s};
							}
						unless ($seq{$d[2]}{'full'} eq 'compl') {	
							if ($seq{$d[2]}{'length'} < $length) {
								$seq{$d[2]} = {'full' => $complete, 'eval' => $d[-1], 'length' => $length, 'seq' => $s};
								}
							if ($seq{$d[2]}{'eval'} > $d[-1]) {
								$seq{$d[2]} = {'full' => $complete, 'eval' => $d[-1], 'length' => $length, 'seq' => $s};
								}
							}	
						}	
					}
				else {
					$seq{$d[2]} = {'full' => $complete, 'eval' => $d[-1], 'length' => $length, 'seq' => $s};
					}
				}
			}
		close(IN);	
		}	
	
	open(OUT, ">$seqfile");
	foreach my $c (keys %seq) {
		print OUT ">", $c, "\n";
		print OUT $seq{$c}{'seq'}, "\n";
		}
	close(OUT);
	}	