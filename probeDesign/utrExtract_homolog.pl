use warnings;
use strict;

my @id = qw(Carlia Lampro1 Lampro2 Sapro);
my @seq = qw(Carlia_S Lampro_C Lampro_S Sapro_S);

for (my $i = 0; $i < scalar(@id); $i++) {
	my $file = '/Users/singhal/Desktop/activeWork/exomeCapture/prerepeat_' . $id[$i] . '/utr.fa';
	
	my $seq = '/Users/singhal/Desktop/activeWork/exomeCapture/seqFiles/' . $seq[$i] . '_trinity.fa.final.annotated.blast';
	
	my @d = `blastall -p blastn -d $seq -i $file -m 8 -a 4 -b 1 -e 1e-20`;
	
	my %matches;
	foreach my $d (@d) {
		my @data = split(/\s+/,$d);
		$matches{$data[0]} = $data[1];
		}
		
	my %seq;	
	open(IN, "<$seq");
	while(<IN>) {
		if ($_ =~ m/>(\S+)/) {
			my $id = $1;
			chomp(my $s = <IN>);
			$seq{$id} = $s;
			}
		}
	close(IN);
	
	my %utr;	
	open(IN, "<$file");
	while(<IN>) {
		if ($_ =~ m/>(\S+)/) {
			my $id = $1;
			chomp(my $s = <IN>);
			$utr{$id} = $s;
			}
		}
	close(IN);
	
	my $out = '/Users/singhal/Desktop/activeWork/exomeCapture/prerepeat_' . $id[$i] . '/finalUTR.fa';
	open(OUT, ">$out");
	
	foreach my $c (keys %matches) {
		my $out1 = 'exonHomolog1.fa';
		open(OUT1, ">$out1");
		print OUT1 ">", $c, "\n" , $utr{$c}, "\n";
		close(OUT1);
	
		my $out2 = 'exonHomolog2.fa';
		open(OUT2, ">$out2");
		print OUT2 ">", $matches{$c}, "\n" , $seq{$matches{$c}}, "\n";
		close(OUT2);
		
		my @call = `exonerate -m affine:local $out1 $out2 --showvulgar yes --showalignment no`;
		
		if ($call[2]) {
			if ($call[2] =~ m/vulgar/) {
				my @d = split(/\s/,$call[2]);
				my $length = abs($d[6] - $d[7]);					
				my $sub; 
				my $start;
			
				if ($d[4] =~ m/\+/ && $d[8] =~ m/\+/) {
					$start = $d[6];
					$sub = substr $seq{$matches{$c}}, $start, $length;
					}
				elsif ($d[4] =~ m/\+/ && $d[8] =~ m/\-/) {
					$start = $d[7];
					$sub = substr $seq{$matches{$c}}, $start, $length;
					$sub = reverse($sub);
					$sub =~ tr/ATGCatgc/TACGtacg/;
					}
				elsif ($d[4] =~ m/\-/ && $d[8] =~ m/\-/) {
					$start = $d[7];
					$sub = substr $seq{$matches{$c}}, $start, $length;
					}
				elsif ($d[4] =~ m/\-/ && $d[8] =~ m/\+/) {
					$start = $d[6];
					$sub = substr $seq{$matches{$c}}, $start, $length;
					$sub = reverse($sub);
					$sub =~ tr/ATGCatgc/TACGtacg/;
					}
					
				print OUT ">", $c, "_1\n", $utr{$c}, "\n";
				print OUT ">", $c, "_2\n", $sub, "\n";
				delete($utr{$c});
				}
			}	
		}	
		
	print OUT ">", $_, "_1\n", $utr{$_}, "\n" foreach keys %utr;	
	close(OUT);	
	my $cleanup = system("rm exonHomolog1.fa exonHomolog2.fa");
	}	