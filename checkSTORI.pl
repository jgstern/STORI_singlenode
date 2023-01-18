#command: perl /path/to/checkSTORI.pl /path/to/file-A /path/to/file-B [-q | -c | -org <file>] [-s]
#-q means quiet mode
#-c means concise, ie omit the taxaIds, -1s, and ?s. NB: If this option is used, gis on the same line will not necessarily be from the same taxon.
#-org <file> means organize the output according to a file, which specifies on each line <taxID> <species.name> <phylum>
#-s means print the score for each family member in addition to the family members

use lib '/nv/hp10/jstern7/perl5reinstall/lib';
use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';
use Data::Dumper;
use List::MoreUtils qw / uniq /;

my $inputA=shift(@ARGV);
my $inputB=shift(@ARGV);
my $quietOrConcise=shift(@ARGV);

my $quietFlag=0; my $conciseFlag=0; my $orgFile=-1; my $orgFlag=0;
if ($quietOrConcise =~ m/-q/) {
	$quietFlag = 1;
}
elsif ($quietOrConcise =~ m/\-c/) {
	$conciseFlag = 1;
}
elsif ($quietOrConcise =~ m/\-org/) {
	$orgFile = shift(@ARGV);
	$orgFlag=1;
}

my $score = shift(@ARGV);
my $scoreFlag=0;
$scoreFlag=1 if ($score =~ m/\-s/);

my @taxaArr=(); my $finishedCount=0;
my @temp = @{GetFams($inputA)}; #read the two datasets into memory
my %famsA = %{$temp[0]};
my %scoresA = %{$temp[1]}; 
my $runtime_a = $temp[2];
#print Dumper \%famsA;
#print Dumper \%scoresA;

my @temp = @{GetFams($inputB)}; #read the two datasets into memory
my %famsB = %{$temp[0]};
my %scoresB = %{$temp[1]}; 
my $runtime_b = $temp[2];
#print Dumper \%scoresB;
@taxaArr = uniq @taxaArr;
my $taxaCount=($#taxaArr + 1);

#Find families between the two datasets that have common sequences
my @pairScores=();
foreach my $famFromA (keys %famsA) {
	my %consensusFam = ();
	my %chosenFam = ();
	my $chosenFamName=-1;
	my $consensusSize = -1;
	foreach my $famFromB (keys %famsB) {
		my $agreementScore = GetAgreementScore($famsA{$famFromA},$famsB{$famFromB});
		my @temp = ($famFromA, $famFromB, $agreementScore);
		push @pairScores, [@temp];
	}
}
my @pairScores_sorted = sort { $b->[2] <=> $a->[2] } @pairScores;  #sort descending
#print Dumper \@pairScores_sorted;

my @famNames=();
my %consensusFams=();
my %consensusScores=();
foreach my $pair_ref (@pairScores_sorted) {
	if ($pair_ref->[2] > 0) {
		my $famA = $pair_ref->[0];
		my $famB = $pair_ref->[1];
		if (exists $famsA{$famA}) {
			if (exists $famsB{$famB}) {
				my @famsAndScores = @{GetConsensusFam($famsA{$famA},$famsB{$famB},$scoresA{$famA},$scoresB{$famB})};
				my %consensus = %{$famsAndScores[0]};
				my %consensus_s = %{$famsAndScores[1]};
				my $consensusName = $famA . $famB;
				push @famNames, $consensusName;
				delete $famsA{$famA};
				delete $famsB{$famB};
				
				#print $consensusName . "\n";
				#print Dumper \%consensus;
				
				$consensusFams{$consensusName} = {%consensus};
				$consensusScores{$consensusName} = {%consensus_s}
			}
		}
	}
}

my %famAgreement=();
foreach my $family (@famNames) {
	my $disagreements=0;
	foreach my $taxon (@taxaArr) {
		if ($consensusFams{$family}{$taxon} eq "?") {
			$disagreements++;
		}
	}
	my $agreement = ($taxaCount - $disagreements);
	my $score = ($agreement/$taxaCount);
	$famAgreement{$family}=$score;
}

#print "quietFlag= $quietFlag\n";
if ($quietFlag==1) {
	#print "finishedCount= $finishedCount\n";
	print "runtime_a: $runtime_a\n";
	print "runtime_b: $runtime_b\n";
	if ($finishedCount==2) { print "bothfinished\n"; }
}

if ($orgFlag==1) {
	print "-1 -1 ";
}

print "score ";
foreach my $family (@famNames) {
	print $famAgreement{$family} . " ";
}
if ($#famNames == -1) {
	print "-1\n"; }

if (($quietFlag==0) && ($conciseFlag==0) && ($orgFlag==0)) {
	print "\n";
	print "family ";
	foreach my $family (@famNames) {
		print $family . " ";
	}
	print "\n";
	
	foreach my $taxon (@taxaArr) {
		print $taxon . " ";
		foreach my $family (@famNames) {	
			if (exists $consensusFams{$family}{$taxon}) {
				print $consensusFams{$family}{$taxon} . " ";
			}
			else {
				print "-1 ";
			}
		}
		print "\n";
	}
	print "\n\n";
	foreach my $taxon (@taxaArr) {
		print $taxon . " " if ($scoreFlag==1);
		foreach my $family (@famNames) {	
			if (exists $consensusScores{$family}{$taxon}) {
				print $consensusScores{$family}{$taxon} . " " if ($scoreFlag==1);
			}
			else {
				print "-1 " if ($scoreFlag==1);
			}
		}
		print "\n" if ($scoreFlag==1);
	}
}
elsif (($quietFlag==0) && ($conciseFlag==1) && ($orgFlag==0)) {
	print "\n";
	print "family ";
	foreach my $family (@famNames) {
		print $family . " ";
	}
	print "\n";
	
	my @conciseGis=();
	
	my $row=0; my $col=0;
	foreach my $family (@famNames) {
		my %temp = %{$consensusFams{$family}};
		foreach my $taxon (keys %temp) {
			if (!($temp{$taxon} =~ m/\?/)) {
				$conciseGis[$row][$col] = $temp{$taxon};
				$row++;
			}
		}
		$col++;
		$row=0;
	}
	
	foreach my $row_ref (@conciseGis) {
		my @row = @{$row_ref};
		my $c=0;
		print "taxon_na ";
		foreach my $fam (@famNames) {
			if (defined $row[$c]) {
				print $row[$c] . " "; }
			else {
				print "na "; }
			$c++;
		}
		print "\n";
	}
	
}
elsif ($orgFlag==1) {
	print "\n";
	print "taxID name phylum ";
	foreach my $family (@famNames) {
		print $family . " ";
	}
	print "\n";
	
	#print "orgFile\= $orgFile\n";
	
	if (!($orgFile eq "-1")) {
		open org, $orgFile;
	}
	
	my @taxids_ordered=();
	my @names_ordered=();
	my @phyla_ordered=();
	
	while (<org>) {
		my $line = $_;
		chomp($line);
		if ($line =~ m/(\d+)\t(.+)\t(.+)/) {
			push @taxids_ordered, $1;
			push @names_ordered, $2;
			push @phyla_ordered, $3;
			#print "file\: $1 $2 $3\n";
		}
	}
	
	my $r=0;
	foreach my $taxon (@taxids_ordered) {
		print $taxon . " " . $names_ordered[$r] . " " . $phyla_ordered[$r] . " ";
		foreach my $family (@famNames) {	
			if (exists $consensusFams{$family}{$taxon}) {
				print $consensusFams{$family}{$taxon} . " ";
			}
			else {
				print "-1 ";
			}
		}
		print "\n";
		$r++;
	}
	print "\n\n";
	my $r=0;
	foreach my $taxon (@taxids_ordered) {
		print $taxon . " " . $names_ordered[$r] . " " . $phyla_ordered[$r] . " " if ($scoreFlag==1);
		foreach my $family (@famNames) {	
			if (exists $consensusScores{$family}{$taxon}) {
				print $consensusScores{$family}{$taxon} . " "  if ($scoreFlag==1);
			}
			else {
				print "-1 " if ($scoreFlag==1);
			}
		}
		print "\n" if ($scoreFlag==1);
		$r++;
	}
	close org;
}

print "\n";

sub GetAgreementScore {
	my %fam1 = %{shift(@_)};
	my %fam2 = %{shift(@_)};
	my $score=0;
	foreach my $taxon (keys %fam1) {
		if (exists $fam2{$taxon}) {
			print "will the regex work correctly on the following?\n";
			print $fam1{$taxon};
			print "\n";
			sleep 30;
			$fam1{$taxon} =~ m/\{*(\d+)\}*/;
			my $gi1 = $1;
			$fam2{$taxon} =~ m/\{*(\d+)\}*/;
			my $gi2 = $1;
			if ($gi1 eq $gi2) {
				$score++;
			}
		}
	}
	return $score;
}

sub GetConsensusFam {
	my %fam1 = %{shift(@_)};
	my %fam2 = %{shift(@_)};
	my %scores1 = %{shift(@_)};
	my %scores2 = %{shift(@_)};
	my %consensus=();
	my %scores=();
	foreach my $taxon (keys %fam1) {
		if (exists $fam2{$taxon}) {
			$fam1{$taxon} =~ m/\{*(\d+)\}*/;
			my $gi1 = $1;
			$fam2{$taxon} =~ m/\{*(\d+)\}*/;
			my $gi2 = $1;
			if ($gi1 eq $gi2) {
				$consensus{$taxon} = $gi1;
				$scores{$taxon} = ($scores1{$taxon} + $scores2{$taxon});
			}
			else {
				$consensus{$taxon} = "?";
				$scores{$taxon} = "?";
			}
		}
		else {
			$consensus{$taxon} = "?";
			$scores{$taxon} = "?";
		}
	}
	foreach my $taxon (keys %fam2) {
		if (exists $fam1{$taxon}) {
		}
		else {
			$consensus{$taxon} = "?";
			$scores{$taxon} = "?";
		}
	}
	my @temp = (\%consensus, \%scores);
	return \@temp;
}

sub GetFams {
	#print "entered GetFams\n";
	my $inputFile = shift(@_);
	print "inputFile is $inputFile\n";
	open (input, $inputFile);
	my @fileArr=<input>;
	my $startLine=-1;
	my $startLine_score=-1;
	my $runtime=-1;
	for (my $r=$#fileArr; $r>=0; $r--) {
		if ($fileArr[$r] =~ m/final\stable\smain/) {
			$startLine = ($r + 2);
			$r=-1;
		}
		elsif ($fileArr[$r] =~ m/score\stable/) {
			$startLine_score = ($r + 2);
			#needn't set r to -1 because the score table always follows "final table main"
		}
	}
	my $finishedFlag=0;
	for (my $r=$#fileArr; $r>=0; $r--) {
		if ($fileArr[$r] =~ m/runtime\s(.+)\s\[total-runtime\](.*)/) {
			#print "match\n";
			$runtime=$1; $theRest = $2;
			if ($theRest =~ m/finished/) { $finishedFlag=1; }
			$r=-1;
		}
	}
	if ($finishedFlag==1) {
		$finishedCount++;
		if ($quietFlag==0) {
			print "$inputFile has finished.\n";
		}
	}
	
	my %fam = ();
	my @famsArr=();
	my %scores = ();
	
	if ($startLine>0) {
		#print "startLine is $startLine\n";
		for (my $r=$startLine; $r<=$#fileArr; $r++) {
			my $line=$fileArr[$r];
			chomp($line);
			#print $line . "\n";
			if ($line=~ m/^taxon/) {
				$line =~ m/taxon\s(.+)/;
				my $famsTemp = $1;
				chomp($famsTemp);
				@famsArr = split / /, $famsTemp;
			}
			elsif ($line =~ m/^\d+/) {
				$line =~ m/(\d+)\s\s(.+)/;
				my $gisTemp = $2;
				chomp($gisTemp);
				my @giArr = split / /, $gisTemp;
				my $taxon = $1;
				push @taxaArr, $taxon;
				foreach my $family (@famsArr) {
					if ($family =~ m/.+/) {
						my $gi = shift(@giArr);
						if ($gi ne -1) {
							$fam{$family}{$taxon} = $gi;
							#print "fam{$family}{$taxon}=$gi\n";
						}
					}
				}
			}
			else {
				$r=($#fileArr + 1);
			}
		}
	}
	if ($startLine_score>0) {
		for (my $r=$startLine_score; $r<=$#fileArr; $r++) {
			my $line=$fileArr[$r];
			chomp($line);
			if ($line=~ m/^taxon/) {
				$line =~ m/taxon\s(.+)/;
				my $famsTemp = $1;
				chomp($famsTemp);
				@famsArr = split / /, $famsTemp;
			}
			elsif ($line =~ m/^\d+/) {
				$line =~ m/(\d+)\s\s(.+)/;
				my $scoresTemp = $2;
				chomp($scoresTemp);
				my @scoresArr = split / /, $scoresTemp;
				my $taxon = $1;
				foreach my $family (@famsArr) {
					if ($family =~ m/.+/) {
						my $score = shift(@scoresArr);
						if ($score != -1) {
							$scores{$family}{$taxon} = $score;
							#print "fam{$family}{$taxon}=$gi\n";
						}
					}
				}
			}
			else {
				$r=($#fileArr + 1);
			}
		}
	}
	
	my @temp = (\%fam, \%scores, $runtime);
	return \@temp;
}
