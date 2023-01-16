#STORI: Selectable Taxa Ortholog Retrieval Iteratively
#command: /path/to/STORI_test.pl runNumber /path/to/source/files/dir windowSize finalMaxFams

#Each iterator (each STORI.pl instance) BLASTs the seeds against taxon databases.
#Within each iterator,
#the results of each search become queries for subsequent searches. Each iterator
#dynamically assigns sequences to different families as the results accumulate.
#By keeping track of the frequency with which queries from each family choose a
#particular sequence as a best hit, STORI assigns each hit to a family.

#!/usr/bin/perl -wT
use lib '/nv/hp10/jstern7/perl5reinstall/lib';
use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';

use strict;
use Fcntl ':flock';
use Data::Dumper;
use List::MoreUtils qw / uniq /;
use Statistics::Descriptive;
use Time::Elapse;


my $hitDir="/home/ec2-user/STORI/universal20150110/hits";            #directory containing taxon-taxon blast results
my $blastdbDir="/home/ec2-user/STORI/universal20150110/blast";
my $blastdbcmdPath = "/home/ec2-user/STORI/blastdbcmd";
my $blastpPath = "/home/ec2-user/STORI/blastp"; 
my $getParentTaxaPath = "/home/ec2-user/STORI/getParentTaxa.pl";

my $runNumber = shift(@ARGV);		
my $sourceFilesDir = shift(@ARGV);
my $windowSize = shift(@ARGV);				#number of taxa in the sliding window = windowSize. 
my $finalMaxFams = shift(@ARGV);			#the # of fams can reach as much finalMaxFams families.
my $sc = shift(@ARGV);

my $blastTemp = $sourceFilesDir . "/tempBlastDB"; #in a few lines we are going to copy the relevant data
								#files from $blastdbDir to $blastTemp, and then set
								#$blastdbDir = $blastTemp. This copying will mean that
								#STORI will access the local volume on the compute node
								#rather than over the network. (or, if running multiple instances of
								#STORI.pl on the same node, this step will prevent attempts to access
								#the same file simultaneously by 2 different processes)
					
my $cmd = "mkdir $blastTemp";
system($cmd);


my $taxaFile = $sourceFilesDir . "/taxa-master[" . $runNumber . "].txt";
my $batchQueryFilePath="$sourceFilesDir/query-gis[" . $runNumber . "-iter0].txt";  #file listing gis that you want orthologs for
my $outputFile="$sourceFilesDir/STORI_out_sc" . $sc . "_" . $runNumber . ".txt";

my $parentTaxaDataPath = "$sourceFilesDir/parentTaxaData.txt";

my %taxon_gi_assigned=();
my $tripleCheck = 3;
my $superDuperRepeatFlag=1;
my %orphans=();
my %familyCensus=();

my $queryFastaPath = $sourceFilesDir . "/temp-query.fasta";
my $blastResultsPath = $sourceFilesDir . "/temp-blast.results";

Time::Elapse->lapse(my $nowTot = "total-runtime");
open (OUT, ">$outputFile");

# assign each gi in the query-gi file to its parent taxon; put it %gi_lookup_hash
# (calls an outside script.)
my %gi_lookup_hash;
my $cmd = "touch $parentTaxaDataPath";
system($cmd);
my $cmd = "rm $parentTaxaDataPath";
system($cmd);

my $cmd = "perl $getParentTaxaPath $batchQueryFilePath $hitDir $taxaFile $sourceFilesDir";
print OUT "cmd is $cmd\n";
system($cmd);
while (!(-e $parentTaxaDataPath)) {
	print OUT "Looking up parent taxa for query gis...";
	sleep 10;
}
sleep 10;
open (parentTaxaData, $parentTaxaDataPath);
while (<parentTaxaData>) {
	my $line = $_;
	my @temp = split / /, $line;
	#gi_lookup_hash{<gi>} = <taxa>;
	chomp($temp[0]); chomp($temp[1]); 
	$gi_lookup_hash{$temp[0]} = $temp[1];
}
print OUT Dumper \%gi_lookup_hash;


#The seed sequences come from a random draw from the results pool of a
#user-initiated keyword search of the natural language annotations of every
#taxon’s protein sequence database. (beginSTORI.pl)


#2: load query gis into %taxon_gi_assigned
open (BATCH_QUERY, "$batchQueryFilePath");
my @queryArr_unparsed = <BATCH_QUERY>; my %tempHash=();
close BATCH_QUERY;

foreach my $line (@queryArr_unparsed) {		
	chomp($line);
	if (!($line =~ m/h\s.+/)) {
		if ($line =~ m/(.+?)\s(.+?)/) {
			#tempHash{<taxID>}{<gi>} = <score>
			$tempHash{$gi_lookup_hash{$1}}{$1} = $2;
		}
	}
	elsif ($line =~ m/h\s(.+)/) {
		$taxon_gi_assigned{$1} = {%tempHash};
		%tempHash=();
	}
}


ShowNumFams(\%taxon_gi_assigned, "after initialization");

#load master taxa list
open (taxa, $taxaFile);
my @taxaArr=<taxa>;
my $count=0;
foreach my $taxon (@taxaArr)
{
	chomp($taxon);  #removes any newline characters at the end of the string
	$taxaArr[$count] = $taxon;
	$count++;
	
	#ok here goes the copying
	my $tempTaxonPIN = $blastTemp . "/$taxon\.pin";
	if (!(-e $tempTaxonPIN)) {
		my $cmd = "cp $blastdbDir/$taxon.* $blastTemp";
		print "$cmd\n";
		system($cmd);
	}
}
print "\n";
$blastdbDir = $blastTemp;
print "using $blastdbDir as database location\n";



my $searchOrderPlaceMarker =0;
my %small_txgi=();
my $avgVisitCount=0;
my $seedDecay=1000;
while ($superDuperRepeatFlag == 1) {		#this loop causes multiple iterations where the taxon search order is reshuffled after each iteration
											#superDuperRepeatFlag is determined by ShuffleTaxaArr - when set to 0, the families have stabilized
											#Iteration repeats until a user-defined time limit expires, or the iterator stops finding new sequences.
	$searchOrderPlaceMarker =0;
	while ($searchOrderPlaceMarker <= ($#taxaArr-($windowSize-1))) {
#							At the beginning of an iteration, STORI chooses a small, user defined number of
#						taxa (usually 4). STORI makes this choice by sliding a window of user-defined
#						size (usually 4) down a list of all the taxa. At iteration 1, the window is at
#						the top of the full taxa list, and the window will advance by one list element
#						with each subsequent iteration.
		%small_txgi = %{GetNextWindow()};
#						 After selecting the small list of taxa, STORI
#						cycles through each quasi-family and BLASTs any sequences assigned to that
#						family, within the small taxa window, against each taxon in the window.
#						STORI
#						parses the BLAST results and assigns any best hits to their parent taxa within
#						the family.
		GetSeqs();
#						  After using BLAST to retrieve best-hits for each quasi-family, STORI
#						identifies best hits assigned to multiple families and prunes all but the most
#						popular.
		PruneAndReassignIntermediate();
#						 After pruning, the small taxa window slides down the full taxa list by
#						one element.
		
		$searchOrderPlaceMarker++;
	}

#   PruneAndReassign checks for "orphan" sequences with a score of 1, possibly
#   indicating the sequence is a pseudo-ortholog (see comments in sub).

	PruneAndReassign();
	
#			STORI takes a few additional steps to produce reasonable orthology predictions.
#		 On the one hand, these steps usually prevent capture by local optima,
#		and on the other, they prevent orthology predictions from careening off to an
#		irrelevant part of family space.
		
#		To prevent family scattering, Merge() may execute once the sliding window
#hits the bottom of the taxa list. This step compares every family with every
#other family, and merges all families above a similarity threshold. To avoid
#undermining the boosted seed scores, merges only occur when the seed scores and
#non-seed scores have similar magnitudes.
		
	print OUT "seedDecay is $seedDecay\n";
	if ($seedDecay <= 1.4) {
		print OUT "merging\n";
		Merge(0.8); }
	
#	RemoveExtraFamilies sorts the families by number of member sequences. If the number of
#families is larger than $finalMaxFams, then RemoveExtraFamilies
#deletes the smallest families until the number of families does not exceed the
#maximum allowable.
	
	RemoveExtraFamilies($finalMaxFams);
	
	ShuffleTaxaArr();
	
	OutputTable();
	
	AddOrphansBackIn();
	
	$seedDecay = ResetAllScores();	#we want to keep scores down, so that we are not relying on inertia to make decisions. otherwise the stochasticity specific to a particular run would cause some sequences to have higher scores than others, biasing the results
									#we set to 2 because any lower integer would risk breaking apart existing families (a score of 1 at the end of the traversal means that sequence becomes its own family)
									#NB: it is OK to reset all scores only if there is one hit per taxon. Since we call ResetAllScores() after PruneAndReassign(), this is the case. However, resetting all scores if there were multiple hits per taxon would mess up the families. 
									#Families must consist of the highest scoring hit per taxon only; equalizing the scores would erase this information. So be careful.
	
	print OUT "\nruntime $nowTot (HH:MM:SS)\n";
	
	$avgVisitCount++;
}

print OUT "\nruntime $nowTot (HH:MM:SS) finished\n";

close OUT;


sub AddOrphansBackIn {
	print OUT "adding orphans back in\n";
	foreach my $name (keys %orphans) {
		my %temp = %{$orphans{$name}};
		$taxon_gi_assigned{$name} = {%temp};
	}
	%orphans=();
}

sub Max {
	my($ref_arr) = @_;
	my @array = @{$ref_arr};
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	my $answer = $stat->max();
	return $answer;
}


sub ResetAllScores {		#Allows seed sequences to stick around as long or as short as the user wishes.
							#Higher initial scores for the seed sequences means they will have a garaunteed presence 
							#for more traversals and the algorithm will be more sensitive to identifying these families if they are sparsely represented (because of, eg, gene loss).
							#	ResetAllScores decreases the
							#seed scores by the value of the highest non-seed score, and resets the non-seed
							#scores to 2. This reset enables non-seed sequences to move between families in
							#the next taxa list traversal.
	
	my $reduction = DecideReduction();
	my @scores=();
	
	foreach my $fam (keys %taxon_gi_assigned) {
		my %fam_spec = %{$taxon_gi_assigned{$fam}};
		foreach my $taxon (keys %fam_spec) {
			my %tax_spec = %{$fam_spec{$taxon}};
			foreach my $hit (keys %tax_spec) {
				$taxon_gi_assigned{$fam}{$taxon}{$hit} -= $reduction;
				if ($taxon_gi_assigned{$fam}{$taxon}{$hit} < 2) {
					$taxon_gi_assigned{$fam}{$taxon}{$hit} = 2; }
				push @scores, $taxon_gi_assigned{$fam}{$taxon}{$hit};
			}
		}
	}
	
	#also report back how close we are to when the seed sequences will lose priority
	# max(@scores)/$reduction
	# this number should initially be >> 1 and with each iteration drop closer to 1
	my $maxSc = Max(\@scores);
	my $ansMain = ($maxSc / $reduction);
	print OUT "ansMain \= maxSc \/ reduction\n";
	print OUT "$ansMain \= $maxSc \/ $reduction\n";
	return $ansMain;
	
	sub DecideReduction {		#Find the highest score that is not from a seed sequence. Seed sequences will have scores that 
								#are substantially higher (eg 500) than a score typically gnerated by the traversal (eg 60). 
								#As we complete more traversals, the seed sequences' scores will decrease until they cannot be distinguished
								#from scores for sequences retrieved by STORI.
		my @scoresArr = ();
		#print OUT "\(DecideReduction\) here is taxon_gi_assigned\n";
		#print OUT Dumper \%taxon_gi_assigned;
		foreach my $fam (keys %taxon_gi_assigned) {
			my %fam_spec = %{$taxon_gi_assigned{$fam}};
			foreach my $taxon (keys %fam_spec) {
				my %tax_spec = %{$fam_spec{$taxon}};
				foreach my $hit (keys %tax_spec) {
					push @scoresArr, $taxon_gi_assigned{$fam}{$taxon}{$hit};
				}
			}
		}
		
		my @scoresArr_s = sort { $b <=> $a } @scoresArr; #descending sort
		my $sd = Sd(\@scoresArr_s);
		my $ultimate1 = shift(@scoresArr_s);
		my $ultimate2 = shift(@scoresArr_s);
		my $ultimate3 = shift(@scoresArr_s);
		my @tempArr = ($ultimate1,$ultimate2,$ultimate3);
		my $ultimate = Mean(\@tempArr);
		my $penultimate1 = shift(@scoresArr_s);
		my $penultimate2 = $scoresArr_s[0];
		my $penultimate3 = $scoresArr_s[1];
		@tempArr = ($penultimate1,$penultimate2,$penultimate3);
		my $penultimate = Mean(\@tempArr);
		my $ans=1;
		
		while ((($ultimate - $penultimate) < $sd) && (@scoresArr_s)) {
			$penultimate1 = shift(@scoresArr_s);
			$penultimate2 = $scoresArr_s[0];
			$penultimate3 = $scoresArr_s[1];
			my @tempArr = ($penultimate1,$penultimate2,$penultimate3);
			$penultimate = Mean(\@tempArr);
		}
		if (($ultimate - $penultimate) > $sd) {
			$ans = $penultimate;
			print OUT "ultimate is $ultimate     penultimate is $penultimate\n";
		}
		else {
			print OUT "ultimate close to penultimate\n";
			$ans = $ultimate;
		}
		if ($ans == 0) { $ans = 2; }	#returning 0 would crash STORI, and 2 is fine, the point is that the reduction should be small
		return $ans;
	}
	
}


sub Sd { #exclude -1
	my($ref_arr) = @_;
	my @tempArr = @{$ref_arr};
	my @array;
	foreach my $elt (@tempArr) {
		if ($elt != -1) {
			push @array, $elt;
		}
	}
	
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	my $answer = $stat->standard_deviation();
	return $answer;
}

sub Mean { #exclude -1
	my($ref_arr) = @_;
	my @tempArr = @{$ref_arr};
	my @array;
	foreach my $elt (@tempArr) {
		if ($elt != -1) {
			push @array, $elt;
		}
	}
	
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	my $answer = $stat->mean();
	return $answer;
}


sub OutputTable {
	print OUT "\n\n\n final table main\n\n";
	my @familyNames = keys %taxon_gi_assigned;
	print OUT "taxon ";
	foreach my $family (@familyNames) {
		print OUT $family . " "; }
	print OUT "\n";
	foreach my $taxon (@taxaArr) {
		print OUT $taxon . "  ";
		foreach my $family (@familyNames) {
			if (exists $taxon_gi_assigned{$family}{$taxon}) {
				print OUT GetTopHit($taxon_gi_assigned{$family}{$taxon}) . " "; }
			else {
				print OUT "-1 "; }
		}
		print OUT "\n";
	}
	
	print OUT "\n\n\n score table\n\n";
	my @familyNames = keys %taxon_gi_assigned;
	print OUT "taxon ";
	foreach my $family (@familyNames) {
		print OUT $family . " "; }
	print OUT "\n";
	foreach my $taxon (@taxaArr) {
		print OUT $taxon . "  ";
		foreach my $family (@familyNames) {
			if (exists $taxon_gi_assigned{$family}{$taxon}) {
				print OUT GetScoreForTopHit($taxon_gi_assigned{$family}{$taxon}) . " "; }
			else {
				print OUT "-1 "; }
		}
		print OUT "\n";
	}
}

sub GetNextWindow {
	my %small_txgi=();
	
	my @someTaxa = @taxaArr[$searchOrderPlaceMarker .. ($searchOrderPlaceMarker+($windowSize-1))];
	
	foreach my $taxon (@someTaxa) {
		foreach my $fam (keys %taxon_gi_assigned) {
			if (exists $taxon_gi_assigned{$fam}{$taxon}) {
				$small_txgi{$fam}{$taxon} = GetTopHit($taxon_gi_assigned{$fam}{$taxon});
				#print OUT "(GetNextWindow) \$small_txgi\{$fam\}\{$taxon\} \= GetTopHit\(see below\)\n";
				#print OUT Dumper $taxon_gi_assigned{$fam}{$taxon}; 
				}
			else {
				$small_txgi{$fam}{$taxon} = -1;
			}
		}
	}
	
	#print OUT "\(GetNextWindow\) " . join(" ", @someTaxa) . "\n";
	
	return \%small_txgi;
}


sub GetTopHit {
	my $hits_ref = shift(@_);
	#print OUT "(GetTopHit) hits_ref= " . $hits_ref . "\n";
	my %hits = %{$hits_ref};
	#print OUT Dumper \%hits;
	my @arr=();
	foreach my $hit (keys %hits) {
		my $score = $hits{$hit};
		my @temp = ($hit, $score);
		push @arr, [@temp];
	}
	my @sorted = sort { $b->[1] <=> $a->[1] } @arr; #descending
	my $ans = $sorted[0][0];
	return $ans;
}

sub GetScoreForTopHit {
	my $hits_ref = shift(@_);
	my %hits = %{$hits_ref};
	my @arr=();
	foreach my $hit (keys %hits) {
		my $score = $hits{$hit};
		my @temp = ($hit, $score);
		push @arr, [@temp];
	}
	my @sorted = sort { $b->[1] <=> $a->[1] } @arr; #descending
	my $ans = $sorted[0][1];
	return $ans;
}



sub GetSeqs {
	#print OUT "(GetSeqs) here is small_txgi \n";
	#print OUT Dumper \%small_txgi;
	
	foreach my $fam (keys %small_txgi) {
		MakeQueryFasta($small_txgi{$fam}, $queryFastaPath);
		my %temp = %{$small_txgi{$fam}};
		
		foreach my $taxon (keys %temp) {
			Blast($queryFastaPath, $taxon, $blastResultsPath);
			ParseAndAssign($blastResultsPath, $taxon, $fam);
		}	
	}	
}


#		In addition to replicating PruneAndReassignIntermediate,
#   PruneAndReassign checks for “orphan” sequences with a score of 1.
# When it identifies an orphan, PruneAndReassignIntermediate
#   moves the sequence to a new
#	family. A score of 1 means that only one of that sequence’s presumed orthologs
#	chose it as a best hit; in this scenario paralogy may be more probable than
#	orthology.
sub PruneAndReassign {
	my $countx=0;
	foreach my $taxon (@taxaArr) {
		my %forConsolidation=();
		foreach my $fam (keys %taxon_gi_assigned) {
			if (exists $taxon_gi_assigned{$fam}{$taxon}) {
				my $gi = GetTopHit($taxon_gi_assigned{$fam}{$taxon});
				my $score = GetScoreForTopHit($taxon_gi_assigned{$fam}{$taxon});
				$forConsolidation{$gi}{$fam} = $score;
				delete $taxon_gi_assigned{$fam}{$taxon};
			}
		}
		my $county=0;
		foreach my $gi (keys %forConsolidation) {
			my $famChoice = GetTopHit($forConsolidation{$gi});
			my $score = GetScoreForTopHit($forConsolidation{$gi});
			if ($score == 1) {	#ie, if after the full traversal the top hit has only been chosen once by surrounding probes, then it might not be a family member and we should make it into its own family
				my $orphname= "orph" . $countx . $county;
				$taxon_gi_assigned{$orphname}{$taxon}{$gi} = 0;
			}
			else {
				$taxon_gi_assigned{$famChoice}{$taxon}{$gi} = $score; }
			$county++;
		}
		$countx++;
	}
}


sub PruneAndReassignIntermediate {		#same as PruneAndReassign, except we do not create any new families in this version- 
										#just do reassignments/pruning within the existing families. we create new families 
										#after the traversal of the taxonomic search order has completed.
	foreach my $taxon (@taxaArr) {
		my %forConsolidation=();
		foreach my $fam (keys %taxon_gi_assigned) {
			if (exists $taxon_gi_assigned{$fam}{$taxon}) {
				my $gi = GetTopHit($taxon_gi_assigned{$fam}{$taxon});
				my $score = GetScoreForTopHit($taxon_gi_assigned{$fam}{$taxon});
				$forConsolidation{$gi}{$fam} = $score;
				delete $taxon_gi_assigned{$fam}{$taxon};
			}
		}
		foreach my $gi (keys %forConsolidation) {
			my $famChoice = GetTopHit($forConsolidation{$gi});
			my $score = GetScoreForTopHit($forConsolidation{$gi});
			$taxon_gi_assigned{$famChoice}{$taxon}{$gi} = $score;
		}
	}
}


sub MakeQueryFasta {
	
	my $hashref=shift(@_);
	my %txgi = %{$hashref};
	
	#print OUT "(begin MakeQueryFasta) here is txgi\:\n";
	#print OUT Dumper \%txgi;
	
	my $queryFastaPath = shift(@_);
	
	my $cmd = "touch $queryFastaPath";
	system($cmd);
	my $cmd = "rm $queryFastaPath";
	system($cmd);
	
	#my @smallTaxaArr = keys %txgi;
	my %lookup=();
	
	foreach my $taxon (keys %txgi) {
		#push @smallTaxaArr, $taxon;
		if ($txgi{$taxon} > -1) {
			$lookup{$txgi{$taxon}} = $taxon;
		}
	}
	
	open (fastaFile, ">$queryFastaPath");
	
	foreach my $gi (keys %lookup) {
		my $taxon = $lookup{$gi};
		if ($gi > 0) {
			my $cmd="$blastdbcmdPath -entry $gi -db " . $blastdbDir . "/$taxon -outfmt \"\%f\"";
			#print OUT "cmd is: " . $cmd . "\n";
			my $seq = qx($cmd);
			#print OUT $seq;
			print fastaFile $seq;
		}
	}
	
	close fastaFile;
	#print OUT "(MakeQueryFasta)\n";
}

sub MakeQueryFasta_v {
	
	my $hashref=shift(@_);
	my %txgi = %{$hashref};
	
	#print OUT "(begin MakeQueryFasta_v) here is txgi\:\n";
	#print OUT Dumper \%txgi;
	
	my $queryFastaPath = shift(@_);
	
	my $cmd = "touch $queryFastaPath";
	system($cmd);
	my $cmd = "rm $queryFastaPath";
	system($cmd);
	
	#my @smallTaxaArr = keys %txgi;
	my %lookup=();
	
	foreach my $taxon (keys %txgi) {
		#push @smallTaxaArr, $taxon;
		if ($txgi{$taxon} > -1) {
			$lookup{$txgi{$taxon}} = $taxon;
		}
	}
	
	open (fastaFile, ">$queryFastaPath");
	
	foreach my $gi (keys %lookup) {
		my $taxon = $lookup{$gi};
		if ($gi > 0) {
			my $cmd="$blastdbcmdPath -entry $gi -db " . $blastdbDir . "/$taxon -outfmt \"\%f\"";
			#print OUT "cmd is: " . $cmd . "\n";
			my $seq = qx($cmd);
			#print OUT $seq;
			print fastaFile $seq;
		}
	}
	
	close fastaFile;
	#print OUT "(MakeQueryFasta)\n";
}

sub Blast {
	my $queryFastaPath = shift(@_);
	my $taxon = shift(@_);
	my $blastResultsPath = shift(@_);
	my $cmd = $blastpPath . " -db " . $blastdbDir . "/$taxon -query $queryFastaPath -evalue 0.05 -outfmt \"6 qgi sgi\" -num_descriptions 10 -num_alignments 10 -parse_deflines";
	#print OUT "cmd is: $cmd\n";
	my $output = qx($cmd);
	open (bres, ">$blastResultsPath");
	print bres $output;
	#print OUT "\(Blast sub result\:\)\n";
	#print OUT $output;
	close bres;
	#$cmd = "cat " . $blastResultsPath . " >> $outputFile";
	#system($cmd);
}

sub ParseAndAssign {
	my $blastResultsPath=shift(@_);
	my $taxon=shift(@_);
	my $fam = shift(@_);
	open (results, $blastResultsPath);
	
	my %blastHash=();
	while (<results>) {
		my @tmp = split;
		my $query = $tmp[0];
		my $subject = $tmp[1];
		if (!(exists $blastHash{$query})) {
			$blastHash{$query} = $subject;
		}
	}
	close results;
	
	my %hits=();
	my @hitsArr = values %blastHash;
	my %hitFreq=();
	foreach my $hit (@hitsArr) {
		if (!(exists $hitFreq{$hit})) {
			$hitFreq{$hit} = 1; }
		else {
			$hitFreq{$hit}++; }
	}
	
	my @hitFreqArr=();
	foreach my $hit (keys %hitFreq) {
		my @temp = ($hit, $hitFreq{$hit});
		push @hitFreqArr, [@temp];
	}
	
	my @hitFreqArr_sorted = sort { $b->[1] <=> $a->[1] } @hitFreqArr;  #sort @hitFreqArr descending by the freq column
	my $gi = $hitFreqArr_sorted[0][0];
	my $score = $hitFreqArr_sorted[0][1];
	
	if (($gi =~ m/.+?/) && ($score =~ m/\d+/)) {
		if (!(exists $taxon_gi_assigned{$fam}{$taxon}{$gi})) {
			$taxon_gi_assigned{$fam}{$taxon}{$gi} = $score;
			#print OUT "\$taxon_gi_assigned\{$fam\}\{$taxon\}\{$gi\} \= $score \n";
			}
		else {
			$taxon_gi_assigned{$fam}{$taxon}{$gi} += $score;
			#print OUT "\$taxon_gi_assigned\{$fam\}\{$taxon\}\{$gi\} \+\= $score \n";
			}
	}
}

sub RemoveExtraFamilies {
	my $maxSize = shift(@_);
	my @famSizes=();
	print OUT "families before deletion: " . join(" ", (keys %taxon_gi_assigned)) . "\n";
	while (my($fam,$fam_tax_ref) = each(%taxon_gi_assigned)) {
		my $memberCount=0;
		while (my($tax,$gi) = each(%{$fam_tax_ref})) {
			$memberCount++;
		}
		my @temp=($fam,$memberCount);
		push @famSizes, [@temp];
	}
	my @famSizes_sorted = sort { $b->[1] <=> $a->[1] } @famSizes;  #sort @famSizes descending by the count column
	my @extraFams;
	if ($#famSizes_sorted >= $maxSize) {
		@extraFams=@famSizes_sorted[$maxSize..$#famSizes_sorted];
	}
	foreach my $fam (@extraFams) {
		delete $taxon_gi_assigned{$fam->[0]};
		#print OUT "deleting family " . $fam->[0] . " ";
	}
	#print OUT "\nremoved $#extraFams extra families \n";
	#print OUT "families after deletion: " . join(" ", (keys %taxon_gi_assigned)) . "\n";
}

sub ShowNumFams {
	my $hash_ref = shift(@_);
	my $comment = shift(@_);
	my %temp = %{$hash_ref};
	my @tempArr = keys %temp;
	my $num = $#tempArr;
	print OUT "\n" . $num . " FAMILIES PRESENT at TIMEPOINT: " . $comment . "\n";
}


sub ShuffleTaxaArr {
	my %taxon_fullness=();
	my @taxa_temp=();
	my @taxaArr_new=();
	
	#step 0- survey the fullness of each family and compare to prev iter size to judge convergence
	my @famSizes=();
	my $tempChange=0;
	my $shuf_choice="na";
	while (my($fam,$fam_tax_ref) = each(%taxon_gi_assigned)) {
		my $memberCount=0;
		while (my($tax,$gi) = each(%{$fam_tax_ref})) {
			$memberCount++;
			if (exists $taxon_fullness{$tax}) { #also keep track of how many seqs are in each taxon, so that this can be used to determine the search order
				$taxon_fullness{$tax}++; }
			else {
				$taxon_fullness{$tax}=1; }
		}
		#print " (ShuffleTaxaArr) fam= $fam and memberCount= $memberCount \n";
		my @temp=($fam,$memberCount);
		push @famSizes, [@temp];
		if (exists $familyCensus{$fam}) {
			$familyCensus{$fam}{"oldCount"} = $familyCensus{$fam}{"newCount"};
			#print "familyCensus{$fam}{oldCount} = familyCensus{$fam}{newCount} = " . $familyCensus{$fam}{"newCount"} . "\n";
		}
		else {
			my %tempCounts=();
			$tempCounts{"oldCount"} = 0.1;
			$tempCounts{"newCount"} = 0.1;
			$tempCounts{"change"} = 1;
			$familyCensus{$fam} = {%tempCounts};
		}
		$familyCensus{$fam}{"newCount"} = $memberCount;
		my $t = abs($familyCensus{$fam}{"newCount"} - $familyCensus{$fam}{"oldCount"});
		if ($familyCensus{$fam}{"oldCount"} == 0) {
			$familyCensus{$fam}{"oldCount"}=0.1; }
		$familyCensus{$fam}{"change"} = ($t / $familyCensus{$fam}{"oldCount"});
		if ($familyCensus{$fam}{"change"} > $tempChange) {
			$tempChange = $familyCensus{$fam}{"change"};
			$shuf_choice = $fam;
		}
	}
	
	print OUT "old taxaArr: " . join(" ", @taxaArr) . "\n";
	
	RandomizeTheReshuffling();
	
	print OUT "new taxaArr: " . join(" ", @taxaArr) . "\n";
	
	#determine whether or not the population of sequences in all families has stabilized
	my $convergenceCriterion = (1 / ($#taxaArr - 1)); #eg, suppose there are 350 taxa. if I went from 349 filled to 350 filled, then I have converged. but if I went from 249 to 250, then I have not converged. if it's 249 to 249, then I have converged. so it's a pretty stringent criterion.
	print OUT "\$tempChange \<\=\? \$convergenceCriterion\n";
	print OUT "$tempChange       $convergenceCriterion\n";
	if ($tempChange <= $convergenceCriterion) {
		#then we are almost ready to stop iteration
		$tripleCheck--;
		print OUT "convergence criterion met; tripleCheck = $tripleCheck \n";
		if ($tripleCheck == 0) {
			$superDuperRepeatFlag = 0;
			print OUT "convergence achieved. ballin outa control. $tempChange <= $convergenceCriterion \n";
		}
	}
	else {					#we require no change 3 times in a row. More precisely, we require change low enough to meet the convergence criterion 3 times in a row.
		$tripleCheck = 3;
	}
	
	sub RandomizeTheReshuffling {
		my @tax_unsorted=();
		my %grouped_by_fullness=();
		foreach my $tax (@taxaArr) {
			if (exists $taxon_fullness{$tax}) {
				if (exists $grouped_by_fullness{$taxon_fullness{$tax}}) {
					my @temp= @{$grouped_by_fullness{$taxon_fullness{$tax}}}; 
					push @temp, $tax;
					$grouped_by_fullness{$taxon_fullness{$tax}} = [@temp];
				}
				else {
					my @temp= ($tax);
					$grouped_by_fullness{$taxon_fullness{$tax}} = [@temp];
				}
			}
			else {
				if (exists $grouped_by_fullness{0}) {
					my @temp= @{$grouped_by_fullness{0}}; 
					push @temp, $tax;
					$grouped_by_fullness{0} = [@temp];
				}
				else {
					my @temp= ($tax);
					$grouped_by_fullness{0} = [@temp];
				}
			}
		}
		foreach my $arr_ref (values %grouped_by_fullness) {
			fisher_yates_shuffle($arr_ref);
		}
		my @fullness_arr= (keys %grouped_by_fullness);
		my @fullness_arr_sorted = sort { $b <=> $a } @fullness_arr; #descending
		
		my @tax_sorted=(); my @taxaArr_new=();
		foreach my $fullness (@fullness_arr_sorted) {
			my @tempArr = @{$grouped_by_fullness{$fullness}};
			foreach my $tax (@tempArr) {
				push @tax_sorted, $tax;
			}
		}
		
		my @tax_sorted2 = @tax_sorted;
		
		#this is the fullish-emptyish method
		my $midpoint = int($#tax_sorted2 / 2);
		my $midpoint2 = ($midpoint + 1);
		my @fullish = @tax_sorted2[0..$midpoint];
		my @emptyish = @tax_sorted2[$midpoint2..$#tax_sorted2];
		fisher_yates_shuffle(\@fullish);
		fisher_yates_shuffle(\@emptyish);
		
		my @taxaArr_new2=();
		while($#taxaArr_new2 != $#tax_sorted2) {
			my $full = shift(@fullish);
			my $empty = shift(@emptyish);
			if (defined $full) {
				push @taxaArr_new2, $full; }
			if (defined $empty) {
				push @taxaArr_new2, $empty; }
		}
		@taxaArr = @taxaArr_new2;
		#end fullish-emptyish method
	}
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


sub CheckObvious {
	my $famA = shift(@_);
	my $famB = shift(@_);
	my $threshold = shift(@_);
	my $ans = "no";
	
	if (!($famA eq $famB)) {
		my $numAgreements=0;
		my $numDisagreements=0;
		
		COMPARE: foreach my $taxon (@taxaArr) {
			my $giA = -1;
			my $giB = -1;
			if (exists $taxon_gi_assigned{$famA}{$taxon}) {
				$giA = GetTopHit($taxon_gi_assigned{$famA}{$taxon});
			}
			if (exists $taxon_gi_assigned{$famB}{$taxon}) {
				$giB = GetTopHit($taxon_gi_assigned{$famB}{$taxon});
			}
			
			if ( ($giA > -1) && ($giB > -1) ) {
				if ($giA eq $giB) {
					$numAgreements++;
				}
				else {
					$numDisagreements++;
				}
			}
		}
		
		my $tot = ($numAgreements + $numDisagreements);
		
		if ( $tot > 0 ) {
			my $s = ($numAgreements / $tot);
			if ($s > $threshold) {
				$ans = "merge";
			}
			elsif ($s < (1 - $threshold)) {
				$ans = "no";
			}
			else {
				$ans = "blast";
			}
		}
		else {
			$ans = "blast";
		}
	}
	return $ans;
}	
	
sub Merge {

	my $cutoff = shift(@_);
	my $mergeFlag = 1;
	
	my %shortlist=();
	
	my @famsList = keys %taxon_gi_assigned;
	#print OUT "number of fams $#famsList\n";
	
	my $tot = $#famsList * $#famsList;
	my $c=0;
	print OUT "Merge first pass\: $tot checks...\n";
	foreach my $famA (@famsList) {
		foreach my $famB (@famsList) {
			my $timesaver = CheckObvious($famA, $famB, $cutoff); #there will be numerous situations where
														#extensive blasting is a waste of time. 
														#see function for details...
			#print OUT "\$timesaver\= $timesaver\n";
			#print OUT "Merge first pass $c\/$tot\n";
			if ($timesaver eq "blast") {
				
				my %famA_hash = %{GetSubset($famA)};
				my %famB_hash = %{GetSubset($famB)};
				
				my @agreementArr1 = @{BlastForAgreement(\%famA_hash, \%famB_hash)};
				my @agreementArr2 = @{BlastForAgreement(\%famB_hash, \%famA_hash)};
				
				my $numAgree = ($agreementArr1[0] + $agreementArr2[0]);
				my $numDisagree = ($agreementArr1[1] + $agreementArr2[1]);
				my $fractionAgree = 0;
				my $total = ($numAgree + $numDisagree);
				
				$fractionAgree = ($numAgree/$total) if ($total > 0);
				
				#print OUT "$numAgree\/$total\=$fractionAgree vs $cutoff\n";
				
				if ($fractionAgree > $cutoff) {
					#MergePair($famA, $famB);
					$shortlist{$famA} = 1;
					$shortlist{$famB} = 1;
				}
			}
			elsif ($timesaver eq "merge") {
				#MergePair($famA, $famB);
				$shortlist{$famA} = 1;
				$shortlist{$famB} = 1;
			}
			$c++;
		}
	}
	MergeRecursive(\%shortlist, $cutoff);
}


sub MergeRecursive {

	my $ref = shift(@_);
	my $cutoff = shift(@_);
	my $mergeFlag = 1;
	
	my %shortlist=%{$ref};
	
	while ($mergeFlag == 1) {
		my @famsList = keys %shortlist;
		#print OUT "number of fams $#famsList\n";
		
		my $mergeTemp = 0;
		
		my $tot = $#famsList * $#famsList;
		my $c=0;
		PAIRWISE:  foreach my $famA (@famsList) {
						foreach my $famB (@famsList) {
							my $timesaver = CheckObvious($famA, $famB, $cutoff); #there will be numerous situations where
																		#extensive blasting is a waste of time. 
																		#see function for details...
							#print OUT "\$timesaver\= $timesaver\n";
							print OUT "MergeRecursive $c\/$tot\n";
							if ($timesaver eq "blast") {
								
								my %famA_hash = %{GetSubset($famA)};
								my %famB_hash = %{GetSubset($famB)};
								
								my @agreementArr1 = @{BlastForAgreement(\%famA_hash, \%famB_hash)};
								my @agreementArr2 = @{BlastForAgreement(\%famB_hash, \%famA_hash)};
								
								my $numAgree = ($agreementArr1[0] + $agreementArr2[0]);
								my $numDisagree = ($agreementArr1[1] + $agreementArr2[1]);
								my $fractionAgree = 0;
								my $total = ($numAgree + $numDisagree);
								
								$fractionAgree = ($numAgree/$total) if ($total > 0);
								
								#print OUT "$numAgree\/$total\=$fractionAgree vs $cutoff\n";
								
								if ($fractionAgree > $cutoff) {
									my $merged = MergePair($famA, $famB);
									if (!($merged eq -1)) {
										delete $shortlist{$famA}; delete $shortlist{$famB}; $shortlist{$merged} = 1;
										$mergeTemp = 1;
										last PAIRWISE;
									}						
								}
							}
							elsif ($timesaver eq "merge") {
								my $merged = MergePair($famA, $famB);
								if (!($merged eq -1)) {
									delete $shortlist{$famA}; delete $shortlist{$famB}; $shortlist{$merged} = 1;
									$mergeTemp = 1;
									last PAIRWISE;
								}	
							}
							$c++;
						}
					}
		
		if ($mergeTemp == 1) {
			$mergeFlag = 1;
		}
		else {
			$mergeFlag = 0;
		}
	}
}



sub MergePair {
	my $famA = shift(@_);
	my $famB = shift(@_);
	
	my %newFam = ();
	my $newFamName = -1;
	
	if ((exists $taxon_gi_assigned{$famA}) && (exists $taxon_gi_assigned{$famB})) {
		foreach my $taxon (@taxaArr) {
								#if the two families being merged either agree on the sequence for a taxon
								#or if one of the families has not chosen a sequence, then use the consensus
								#choice for the merged family. However, if the two families disagree, then make
								#two new orphan families, each with one of the disagreeing sequences.
			if (exists $taxon_gi_assigned{$famA}{$taxon}) {
				my $giA = GetTopHit($taxon_gi_assigned{$famA}{$taxon});
				my $scoreA = GetScoreForTopHit($taxon_gi_assigned{$famA}{$taxon});
				
				if (!(exists $taxon_gi_assigned{$famB}{$taxon})) {
					$newFam{$taxon}{$giA} = $scoreA;
				}
				else {			
					my $giB = GetTopHit($taxon_gi_assigned{$famB}{$taxon});
					if ($giB == $giA) {
						$newFam{$taxon}{$giA} = $scoreA;
					}
					else {
						my $orphname= "orph" . $giB;
						$orphans{$orphname}{$taxon}{$giB} = 0;
						my $orphname= "orph" . $giA;
						$orphans{$orphname}{$taxon}{$giA} = 0;
					}
				}
			}
			
			elsif (exists $taxon_gi_assigned{$famB}{$taxon}) {
				my $giB = GetTopHit($taxon_gi_assigned{$famB}{$taxon});
				my $scoreB = GetScoreForTopHit($taxon_gi_assigned{$famB}{$taxon});
				
				if (!(exists $taxon_gi_assigned{$famA}{$taxon})) {
					$newFam{$taxon}{$giB} = $scoreB;
				}
				else {
					my $giA = GetTopHit($taxon_gi_assigned{$famA}{$taxon});
					if ($giA == $giB) {
						$newFam{$taxon}{$giB} = $scoreB;
					}
					else {
						my $orphname= "orph" . $giB;
						$orphans{$orphname}{$taxon}{$giB} = 0;
						my $orphname= "orph" . $giA;
						$orphans{$orphname}{$taxon}{$giA} = 0;
					}
				}
			}
			
		}
	
		#print OUT "here is famA\n";
		#print OUT Dumper $taxon_gi_assigned{$famA};
		#print OUT "here is famB\n";
		#print OUT Dumper $taxon_gi_assigned{$famB};
	
		delete $taxon_gi_assigned{$famA};
		delete $taxon_gi_assigned{$famB};
		$famA .= "m";
		$taxon_gi_assigned{$famA} = {%newFam};
		$newFamName = $famA;
		#print OUT "here is the merged fam\n";
		#print OUT Dumper $taxon_gi_assigned{$famA};
	}
	return $famA;
}

sub BlastForAgreement {
	my $Aref = shift(@_);
	my $Bref = shift(@_);
	
	my %query = %{$Aref};
	my %subject = %{$Bref};
	#print OUT "here is query hash\n";
	#print OUT Dumper \%query;
	#print OUT "here is subject hash\n";
	#print OUT Dumper \%subject;
	
	MakeQueryFasta_v(\%query, $queryFastaPath);
	
	my $agreements = 0; my $disagreements = 0;
	
	foreach my $taxon (keys %subject) {
		Blast($queryFastaPath, $taxon, $blastResultsPath);
		my @hits = @{GetHits($blastResultsPath)};
		#print OUT "here are the hits: " . join(" ", @hits) . "\n";
		my $subj = $subject{$taxon};
		foreach my $hit (@hits) {
			if ($hit == $subj) {
				$agreements++; }
			else {
				$disagreements++; }
		}
	}
	
	my @ans = ($agreements, $disagreements);
	return \@ans;

	sub GetHits {
		my $blastResultsPath=shift(@_);	
		open (GETHITSR, $blastResultsPath);
		
		my %blastHash=();
		while (<GETHITSR>) {
			my @tmp = split;
			my $query = $tmp[0];
			my $subject = $tmp[1];
			if (!(exists $blastHash{$query})) {
				$blastHash{$query} = $subject;
			}
		}
		close GETHITSR;
		
		#my %hits=();
		#print OUT "(BlastForAgreement-GetHits) here are the blast results\n";
		#my $cmd = "cat $blastResultsPath";
		#my $temp = qx($cmd);
		#print OUT $temp . "\n";
		
		#print OUT "(BlastForAgreement-GetHits) here is blastHash\n";
		#print OUT Dumper \%blastHash;
		my @hitsArr = values %blastHash;
		return \@hitsArr;
	}
}



sub GetSubset {		#extracts target fam from %taxon_gi_assigned, but at most 10 taxa, returns as a hash ref
	my $target = shift(@_);
	my %ans=();
	
	my @tempArr=();
	
	foreach my $taxon (@taxaArr) {
		if (exists $taxon_gi_assigned{$target}{$taxon}) {
			my $gi = GetTopHit($taxon_gi_assigned{$target}{$taxon});
			my @temp = ($taxon, $gi);
			push @tempArr, [@temp];
		}
	}
	
	if ($#tempArr >= 0) {
		fisher_yates_shuffle(\@tempArr);
	}
	
	
	for (my $i=0; $i<10; $i++) {
		if ($i <= $#tempArr) {
			my @row = @{$tempArr[$i]};
			$ans{$row[0]} = $row[1];
		}
	}
	
	return \%ans;
}
