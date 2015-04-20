#use lib '/nv/hp10/jstern7/perl5reinstall/lib';
#use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';
use strict; use warnings;
use Fcntl ':flock';
use Data::Dumper;
use List::MoreUtils qw / uniq /;
use Statistics::Descriptive;

my %indexedNames_reverse=();
my $STORIdir = "/home/ec2-user/STORI/";
my $arrayFile = $STORIdir . "job_data_STORI.txt";
my $blastdbcmdPath = "/home/ec2-user/STORI/blastdbcmd";
my $blastdbDir="/home/ec2-user/STORI/universal20150110/blast";
my $clustalw_loc="/home/ec2-user/clustalw21/clustalw2";
my $clustalo_loc="/home/ec2-user/clustalo/clustalo";
	
my %converged=();
my %unconverged=();
my %paused=();
my %fin_=();
my %bact_fin_=();
my %arch_fin_=();
my %euk_fin_=();
my $repeat=1;
my %index_master=();
my %groups=();
LoadGroups();
my %runHash = %{LoadRunData($arrayFile)};


#print "here is the converged hash after LoadRunData\n";
#print Dumper \%converged;
#print "here is groups after LoadRunData\n";
#print Dumper \%groups;

my %run_params = %{LoadParams($arrayFile)};

my %clipboard_reverse_hash=();
my %clipboard_hash=();
my %clipboard=();

ShowMenu();


while ($repeat==1) {
	
	
	my $input="";
	$SIG{"ALRM"} = sub {die "TimeOut"};  #keep looping even without input
	eval {alarm 300; print "STORI> "; $input = <>; alarm 0};
	if ($@) {print "\n";}
	else {
		#print "You entered $input";
		#print "\n";
	}
	
	
	
	if ($input =~ m/show\s(.+)/) {
		if ($1 =~ m/groups/) {
			print "showing the groups\n";
			print Dumper \%groups;
			#print Dumper \%converged;
			}
		elsif ($1 =~ m/runs/) {
			print "showing the runs\n"; 
			%index_master = %{ShowRuns()};
			#Print2Darr(\@runArr);
			#print Dumper \%runHash;
			#print Dumper \%run_params;
			}
		elsif ($1 =~ m/clipboard\s?(.*)/) {
			my $index=$1;
			print "index\= $index\n";
			if($index =~ m/^\d+$/) {
				ShowClipboard($index); }
			elsif ($index =~ m/^\-all\s?(.*)/) {				#allows output of the clipboard contents so that sequences are matched
															# with their parent taxon id, name and phylip abbreviation
				my $orgFile = $1;
				print "orgFile\= $orgFile\n";
				if ($orgFile =~ m/(\S+)\s\-f\s(.+)/) {		#this is a special case where the user did retrievals using subsampled
															#taxa lists, and wishes to only output a table containing those subsampled taxa,
															# and did not make separate 'subsample.txt' files, and instead prefers to
															#filter the main taxa list (eg 'archaea.txt') through the 'taxa-master[subsample].txt' file,
															# which contains only taxids but no nomenclature
															# e.g., see figure 2 of Stern & Gaucher, 2013
					my $org2 = $1;
					my $filter = $2;
					$filter = "taxa-master\[" . $filter . "\].txt";
					print "filtering orgFile\= $org2 through $filter\n";
					ShowOrganizedFilteredClipboard($org2, $filter);
				}
				else {
					print "showing entire clipboard using org file $orgFile\n";
					ShowOrganizedClipboard($orgFile);
				}
			}
			elsif ($index =~ m/\-fulltaxa\s(.+)/) {
				my $orgFile = $1;
				print "showing clipboard taxa for which all families are present, using org file $orgFile\n";
				ShowFullClipboard($orgFile);
			}
			else {
				ShowClipboard();
			}
		}
		else {
			my $name = $1;
			if (exists $groups{$name}) {
				ShowTheseRuns($groups{$name});
			}
			elsif (exists $runHash{$name}) {
				ShowRunResult($name);
			}
			elsif ($name =~ m/(.+)\s\-c/) {
				ShowRunResultConcise($1);
			}
			elsif ($name =~ m/(.+)\s\-org\s(.+)/) {
				ShowRunResultOrganized($1, $2);
			}
		}
	}
	elsif ($input =~ m/^group\s(\S+)\s\-f\s(\S+)\sthru\s(.+)/) {
		my $newName = $1;
		my $groupToFilter = $2;
		my $expr = $3;
		print "filter option.\n newName\= $newName\n groupToFilter\= $groupToFilter\n expr\= $expr\n";
		
		MakeFilteredGroup($newName, $groupToFilter, $expr);
	}
	
	elsif ($input =~ m/^group\s(.+)/) {
		my @input=split / /, $1;
		my $groupName=shift(@input);
		print "grouping groupName=$groupName groups=" . join(" ", @input) . "\n";
		
		Group($groupName,\@input);
	}
	elsif ($input =~ m/annotate\s(.+)/) {
		if (exists $clipboard_hash{$1}) {
			my $ref = $clipboard{$clipboard_hash{$1}};
			GetAnnotation($ref);
		}
	}
	elsif ($input =~ m/^name\s(.+)/) {
		my $ref = $clipboard{$clipboard_hash{$1}};
		DoWordCensus($ref);
	}
	elsif ($input =~ m/Name\s(.+)/) {
		AssignNames($1);
	}
	elsif ($input =~ m/rename\s(\d+)\s(\S+)/) {
		if (exists $clipboard_hash{$1}) {
			my $ref = $clipboard{$clipboard_hash{$1}};
			my $oldName = $clipboard_hash{$1};
			my $newName = $2;
			$clipboard_hash{$1} = $newName;
			$clipboard_reverse_hash{$newName} = $1;
			delete $clipboard_reverse_hash{$oldName};
			Rename($ref, $oldName, $newName);
		}
		else {
			print "Sorry, index not found.\n";
		}
	}
	elsif ($input =~ m/stats\s(.+)/) {
		my $name=$1;
		print "showing stats for run/group $name\n";
		if (exists $groups{$name}) {
			GroupStats($name);
		}
		else {
			my $score = GetRunScore($name);
			my $time = GetRunTime($name);
			my $count = GetNumFamilies($name);
			my @scoreTrend = @{GetAllRunScores($name)};
			my @scoreTrend_r = @{Round(@scoreTrend)};
			my @famTrend = @{GetNumFamsTrend($name)};
			my @famTrend_r = @{Round(@famTrend)};
			print "score: " . $score . "\ntime: " . $time . "\nnumFamilies: " . $count . "\n";
			print join(" \| ", @scoreTrend_r) . "\n";
			print join(" \| ", @famTrend_r) . "\n";
		}

	}
	elsif ($input =~ m/ungroup\s(.+)/) {
		my $groupName=$1;
		print "ungrouping $groupName\n";
		delete $groups{$groupName};
	}
	elsif ($input =~ m/clear\s(.+)/) {
		my $indices=$1;
		my @forClearing = @{IndexParser($indices)};
		Clear(\@forClearing);
	}
	elsif ($input =~ m/detail\s(.+)/) {
		my $name=$1;
		print "details $name\n";
		my $ref = $runHash{$name};
		print Dumper $ref;
	}
	elsif ($input =~ m/compare\s(.+)\s(.+)/) {
		my $run1=$1;
		my $run2=$2;
		PairwiseRunCompare($run1,$run2);
	}
	elsif ($input =~ m/getid\s(.+)/) {
		my $name = $1;
		if (exists $groups{$name}) {
			GetID($name);
		}
	}
	elsif ($input =~ m/parse\s(.+)/) {
		IndexParser($1);
	}
	elsif ($input =~ m/summarize\s(.+)/) {
		my $s_in = $1;
		my $orgFile = -1; my $run=-1;
		if ($s_in =~ m/(.+)\s\-org\s(.+)/) {
			$run = $1;
			$orgFile = $STORIdir . $2;
			SummarizeFiltered($run, $orgFile);
		}
		else {
			$run = $s_in;
			if ($run =~ m/(\s|\.\.|\-)/) {
				my @indices = @{IndexParser($run)};
				foreach my $index (@indices) {
					my $name = $index_master{$index};
					Summarize($name);
				}
			}
			elsif (exists $groups{$run}) {
				my @runList = @{$groups{$run}};
				foreach my $name (@runList) {
					Summarize($name);
				}
			}
			else {
				Summarize($run);
			}
		}
	}
	elsif ($input =~ m/distance\s(identity|bitscore|length|sd)\s(.+)/) {
		my $measure = $1;
		my $indices = $2;
		PrintDistanceMatrix($measure, $indices);
		#distance <identity \|\| bitscore> <clipboard indices>
	}
	elsif ($input =~ m/menu/) {
		ShowMenu();	}
	elsif ($input =~ m/overlap/) {
		FindOverlappingFamilies();
	}
	elsif ($input =~ m/export\sclipboard\s(\S+)\s(\S+)\s(\S+)/) {
		my $famSuffix = $1;
		my $fileName = $2;
		my $orgFile = $3;
		PrintToFileOrganizedClipboard($famSuffix, $fileName, $orgFile);
	}
	elsif ($input =~ m/exit/) {
		$repeat=0;
	}
	SaveGroups();
	LoadGroups();
}
print "bye!\n";



sub ShowMenu {
	print "Welcome to STORI Stats!
	\nMulti-node STORI ©2012 JG Stern & EA Gaucher
Single-node STORI ©2015 JG Stern
	\n\t\tCommands:
	\t\tshow <\'runs\' \|\| \'groups\' \|\| \'clipboard\' \[index \|\| \-all <org_file> \|\| \-fulltaxa <org_file>\] \|\| run-name [-c \|\| -org <file>]>
	\t\tgroup <new-group-name> \(<1 or more run-names> \|\| \-f <existing_group> thru <regexp>\)
	\t\tstats <run-name \|\| group-name>
	\t\tungroup <group-name>
	\t\tdetail <run-name>
	\t\tgetid <group>
	\t\tsummarize <run-name | indices | group> \[\-org <org_file> \(iff summarizing a single run\)\]
	\t\tclear <clipboard indices>
	\t\tdistance <identity \|\| bitscore \|\| sd> <clipboard indices>
	\t\tscramble <clipboard indices>
	\t\tannotate <clipboard index>
	\t\trename <clipboard index> <new_name>
	\t\tName <cutoff> (automatically renames families in clipboard\; moderately accurate)
	\t\toverlap (identifies overlaps between families on the clipboard)
	\t\texport clipboard <family-suffix> <file_name \(no path\)> <org_file \(no path\)>
	\t\texit\n";
}


sub Rename {
	#$ref, $oldName, $newName
	my $ref = shift(@_);
	my $oldName = shift(@_);
	my $newName = shift(@_);
	my %txgi = %{$ref};
	delete $clipboard{$oldName};
	$clipboard{$newName} = {%txgi};
}

sub PrintDistanceMatrix {
	my $measure = shift(@_);
	my $indices = shift(@_);
	my @indexArr = @{IndexParser($indices)};
	print "measure=$measure\n";
	print "indices= " . join(" ", @indexArr) . "\n";
	
	my %numberedFams=();
	my @famArr=();
	my $counter=0;
	foreach my $index (@indexArr) {
		if (exists $clipboard_hash{$index}) {
			my $fam = $clipboard_hash{$index};
			$numberedFams{$counter} = $fam;
			push @famArr, $fam;
			$counter++;
		}
	}
	
	my %allpairs_lo=(); my %allpairs_mean=(); my %allpairs_hi=();
	foreach my $numA (keys %numberedFams) {
		foreach my $numB (keys %numberedFams) {
			#if ($numA != $numB) {
				my @temp = ($numA, $numB);
				my @temp_s = sort { $a <=> $b } @temp;  #sort ascending
				my $famA = $numberedFams{$temp_s[0]};
				my $famB = $numberedFams{$temp_s[1]};
				my $famA_ref = $clipboard{$famA};
				my $famB_ref = $clipboard{$famB};
				if (!(exists $allpairs_lo{$famA}{$famB})) {
					print "getting dist of $famA vs $famB\n";
					my $ref = Distance($famA_ref, $famB_ref, $measure);
					$allpairs_lo{$famA}{$famB} = $ref->[0];
					$allpairs_mean{$famA}{$famB} = $ref->[1];
					$allpairs_hi{$famA}{$famB} = $ref->[2];
				}
			#}
		}
	}
	
	#print the lo matrix
	print "-2stdev\n";
	foreach my $famA (@famArr) {
		foreach my $famB (@famArr) {
			if (exists $allpairs_lo{$famA}{$famB}) {
				print $allpairs_lo{$famA}{$famB} . " "; }
			else {
				print "- ";
			}
		}
		print "\n";
	}
	
	#print the mean matrix
	print "mean\n ";
	foreach my $fam (@famArr) {
		print $fam . " "; }
	print "\n";
	foreach my $famA (@famArr) {
		print $famA . " ";
		foreach my $famB (@famArr) {
			if (exists $allpairs_mean{$famA}{$famB}) {
				print $allpairs_mean{$famA}{$famB} . " "; }
			else {
				print "- ";
			}
		}
		print "\n";
	}
	
	#print the hi matrix
	print "+2stdev\n";
	foreach my $famA (@famArr) {
		foreach my $famB (@famArr) {
			if (exists $allpairs_hi{$famA}{$famB}) {
				print $allpairs_hi{$famA}{$famB} . " "; }
			else {
				print "- ";
			}
		}
		print "\n";
	}
}



sub Distance {
	my $famA_ref = shift(@_);
	my $famB_ref = shift(@_);
	my $measure = shift(@_);
	
	my @ans=(-1,-1,-1);
	
	#print "finding distance between\n";
	#print Dumper $famA_ref;
	#print Dumper $famB_ref;
	
	my $output;
	
	if ($measure =~ m/sd/) {
		my $locA = $STORIdir . "tempA.fasta";
		MakeCleanFasta($famA_ref, $famB_ref, $locA);  #make a fasta with acc format gi_A or gi_B
		my $belvuIn = $STORIdir . "tempA.belvuIn";
		AlignForBelvu($locA, $belvuIn);  #align the fasta using clustalw and format the output for belvu
		$output = RunScoredist($belvuIn);  #run the scoredist algorithm in Belvu and compute the mean fam dist ± 2 SEM
		print $output . "\n";
	}
	else {
		my $locA = $STORIdir . "tempA.fasta";
		my $locB = $STORIdir . "tempB.fasta";
		MakeFasta($famA_ref, $locA);
		MakeFasta($famB_ref, $locB);
		
		my $cmd = "perl $STORIdir" . "smith-waterman-scripts/familyPairDistanceSW-2.pl " . $locA . " " . $locB;
		$output = qx($cmd);
		
		print $output . "\n";
	}
	
	if ($measure =~ m/identity/) {
		#report distance as %identity
		if ($output =~ m/identity\s(\S+)\s(\S+)\s(\S+)\n/) {
			my $lo = $1; my $mean = $2; my $hi = $3;
			$ans[0] = $lo;
			$ans[1] = $mean;
			$ans[2] = $hi;
		}
	}
	elsif ($measure =~ m/bitscore/) {
		#report distance as bitscore
		if ($output =~ m/bitscore\s(\S+)\s(\S+)\s(\S+)\n/) {
			my $lo = $1; my $mean = $2; my $hi = $3;
			$ans[0] = $lo;
			$ans[1] = $mean;
			$ans[2] = $hi;
		}
	}
	elsif ($measure =~ m/length/) {
		#report alignment length
		if ($output =~ m/length\s(\S+)\s(\S+)\s(\S+)\n/) {
			my $lo = $1; my $mean = $2; my $hi = $3;
			$ans[0] = $lo;
			$ans[1] = $mean;
			$ans[2] = $hi;
		}
	}
	elsif ($measure =~ m/sd/) {
		#report alignment length
		if ($output =~ m/scoredist\s(\S+)\s(\S+)\s(\S+)\n/) {
			my $lo = $1; my $mean = $2; my $hi = $3;
			$ans[0] = $lo;
			$ans[1] = $mean;
			$ans[2] = $hi;
		}
	}
	return \@ans;
}


sub RunScoredist {
	my $belvuIn = shift(@_);
	my $matrixFile = $STORIdir . "matrix.out";
	my $cmd = $STORIdir . "belvu.LIN_64bit -T b -T p $belvuIn -o tree > $matrixFile";
	print "cmd is $cmd\n";
	system($cmd);
	
	open MATRIX, "$matrixFile";
	my %distances=();
	my @seqNames=();
	my $row=-1;
	while (<MATRIX>) {
		my $line=$_;
		chomp($line);
		my @lineArr = split /\t/, $line;
		#pop(@lineArr);
		
		my $col=0;
		if ($lineArr[0] =~ m/\//) {
			my $index=0;
			foreach my $elt (@lineArr) {
				$lineArr[$index] =~ s/(.+)\/.+/$1/;
				$index++;
			}
			@seqNames = @lineArr;
			#print join("\n", @seqNames);
		}
		elsif ($line =~ m/\d+/) {
			foreach my $dist (@lineArr) {
				$dist =~ m/\s*(\S+\.\S+)/;
				$dist = $1;
				#print "$row $col $dist\n";
				$distances{$seqNames[$row]}{$seqNames[$col]}=$dist;
				$col++;
			}
		}
	
		$row++;
	}
	close MATRIX;
	
	#print Dumper \%distances;
	
	my @AvsB=(); my @BvsA=();
	
	foreach my $gi1 (keys %distances) {
		my %temp = %{$distances{$gi1}};
		#print "$gi1\n";
		#print Dumper \%temp;
		foreach my $gi2 (keys %temp) {
			#print "$gi1 $gi2\n";
			if (($gi1 =~ m/A*/) && ($gi2 =~ m/B*/)) {
				push(@AvsB, $temp{$gi2});
			}
			elsif (($gi1 =~ m/B*/) && ($gi2 =~ m/A*/)) {
				push(@BvsA, $temp{$gi2});
			}
		}
	}
	
	#print Dumper \@AvsB;
	#print Dumper \@BvsA;
	
	my @alldist = (@AvsB, @BvsA);
	my @ans = @{DoStats(\@alldist)};
	my $final = "scoredist $ans[0] $ans[1] $ans[2]\n";
	return $final;
}


sub DoStats {
	my $arr_ref = shift(@_);
	
	my @temp = @{$arr_ref};
	my $n = $#temp;
	my $sqrtn = sqrt($n);
	
	my $mean = Mean(\@temp);
	my $sd = Sd(\@temp);
	my $stderr = ($sd/$sqrtn);
	my $loBar = RoundScalar(($mean - (2*$stderr)));
	my $hiBar = RoundScalar(($mean + (2*$stderr)));
	$mean = RoundScalar($mean);
	
	#print "doStats  $loBar  $mean  $hiBar\n";
	my @ans = ($loBar, $mean, $hiBar);
	return \@ans;
}


sub AlignForBelvu {
#($locA, $belvuIn)
	my $locA = shift(@_);
	my $belvuIn = shift(@_);
	
	my $temp = $STORIdir . "temp.aln";
	my $cmd = "$clustalw_loc -INFILE=$locA -ALIGN -TYPE=PROTEIN -OUTFILE=$temp -OUTPUT=CLUSTAL -OUTORDER=ALIGNED -MATRIX=GONNET -GAPOPEN=8 -TYPE=PROTEIN -QUIET";
	print "cmd is $cmd\n";
	system($cmd);

	$cmd = "$clustalo_loc -i $temp -o $belvuIn --outfmt=selex --force";
	print "cmd is $cmd\n";
	system($cmd);
}

sub MakeCleanFasta {
#($famA_ref, $famB_ref, $locA)
	my $famA_ref = shift(@_);
	my %txgiA = %{$famA_ref};
	my $famB_ref = shift(@_);
	my %txgiB = %{$famB_ref};
	my $filepath = shift(@_);
	
	my $cmd = "touch $filepath";
	system($cmd);
	$cmd = "rm $filepath";
	system($cmd);
	
	my %lookupA=(); my %lookupB=();
	
	foreach my $taxon (keys %txgiA) {
		if ($txgiA{$taxon} > -1) {
			$lookupA{$txgiA{$taxon}} = $taxon;
		}
	}
	foreach my $taxon (keys %txgiB) {
		if ($txgiB{$taxon} > -1) {
			$lookupB{$txgiB{$taxon}} = $taxon;
		}
	}
	
	open (fastaFile, ">>$filepath");
	
	foreach my $gi (keys %lookupA) {
		my $taxon = $lookupA{$gi};
		if ($gi > 0) {
			my $cmd="$blastdbcmdPath -entry $gi -db " . $blastdbDir . "/$taxon -outfmt \"\%f\"";
			#print out "cmd is: " . $cmd . "\n";
			my $seq = qx($cmd);
			#print out $seq;
			my $replacementAcc = "\>A\_" . $gi . "\n";
			$seq =~ s/\>.+\n/$replacementAcc/;
			print fastaFile $seq;
		}
	}
	foreach my $gi (keys %lookupB) {
		my $taxon = $lookupB{$gi};
		if ($gi > 0) {
			my $cmd="$blastdbcmdPath -entry $gi -db " . $blastdbDir . "/$taxon -outfmt \"\%f\"";
			#print out "cmd is: " . $cmd . "\n";
			my $seq = qx($cmd);
			#print out $seq;
			my $replacementAcc = "\>B\_" . $gi . "\n";
			$seq =~ s/\>.+\n/$replacementAcc/;
			print fastaFile $seq;
		}
	}
	
	close fastaFile;
	
}

sub MakeFasta {
	
	my $fam_ref = shift(@_);
	my %txgi = %{$fam_ref};
	
	#print out "(begin MakeFasta) here is txgi\:\n";
	#print out Dumper \%txgi;
	
	my $filepath = shift(@_);
	
	my $cmd = "touch $filepath";
	system($cmd);
	$cmd = "rm $filepath";
	system($cmd);
	
	my %lookup=();
	
	foreach my $taxon (keys %txgi) {
		if ($txgi{$taxon} > -1) {
			$lookup{$txgi{$taxon}} = $taxon;
		}
	}
	
	open (fastaFile, ">>$filepath");
	
	foreach my $gi (keys %lookup) {
		my $taxon = $lookup{$gi};
		if ($gi > 0) {
			my $cmd="$blastdbcmdPath -entry $gi -db " . $blastdbDir . "/$taxon -outfmt \"\%f\"";
			#print out "cmd is: " . $cmd . "\n";
			my $seq = qx($cmd);
			#print out $seq;
			print fastaFile $seq;
		}
	}
	
	close fastaFile;
	#print out "(MakeFasta)\n";
}


sub Clear {
	my $arr_ref = shift(@_);
	my @indices = @{$arr_ref};
	print "Clearing from clipboard\: " . join(" ", @indices) . "\n";
	foreach my $index (@indices) {
		delete $clipboard{$clipboard_hash{$index}};
		delete $clipboard_hash{$index};
	}
}


sub FindOverlappingFamilies {
	#Identify families on the clipboard that share members
	my @pairScores=();
	my %fortest=();
	
	foreach my $famA (keys %clipboard) {
		foreach my $famB (keys %clipboard) {
			if (!($famA eq $famB)) {
				my @temp = ($famA, $famB);
				my @temp_s = sort {$a cmp $b} @temp;
				my $concat = $temp_s[0] . $temp_s[1];
				$fortest{$concat} = [@temp_s];
			}
		}
	}
	
	foreach my $arr_ref (values %fortest) {
		my $famA = $arr_ref->[0];
		my $famB = $arr_ref->[1];
		my $agreementScore = GetAgreementScore($clipboard{$famA},$clipboard{$famB});
		my $sizeA = GetFamSize($clipboard{$famA});
		my $sizeB = GetFamSize($clipboard{$famB});
		my @temp = ($famA, $famB, $agreementScore, $sizeA, $sizeB);
		push @pairScores, [@temp];
	}
	
	my @pairScores_sorted = sort { $b->[2] <=> $a->[2] } @pairScores;  #sort descending
	
	#print Dumper \%clipboard_reverse_hash;
	
	my @forClearing = ();
	foreach my $line_ref (@pairScores_sorted) {
		if ($line_ref->[2] > 0.5) {
			my $famA = $line_ref->[0];
			my $famB = $line_ref->[1];
			my $sizeA = $line_ref->[3];
			my $sizeB = $line_ref->[4];
			my $famAid = $clipboard_reverse_hash{$famA};
			my $famBid = $clipboard_reverse_hash{$famB};
			print "$famA \(i\=$famAid\ | s\=$sizeA) $famB \(i\=$famBid\ | s\=$sizeB) $line_ref->[2]\n";
			if ($sizeA > $sizeB) {
				push @forClearing, $famBid; }
			if ($sizeB > $sizeA) {
				push @forClearing, $famAid; }
		}
	}
	Clear(\@forClearing);
}


sub GetAgreementScore {	#returns the #agreements normalized to the #taxa with a member in both famimlies
	my %fam1 = %{shift(@_)};
	my %fam2 = %{shift(@_)};
	my $score=0;
	my $size=0;
	my @temp = keys %fam1;
	my $size1=($#temp + 1);
	@temp = keys %fam2;
	my $size2=($#temp + 1);
	if ($size1 < $size2) {
		$size=$size1;}
	else {
		$size=$size2;}
	
		
	foreach my $taxon (keys %fam1) {
		if (exists $fam2{$taxon}) {
			$fam1{$taxon} =~ m/\{*(\d+)\}*/;
			my $gi1 = $1;
			$fam2{$taxon} =~ m/\{*(\d+)\}*/;
			my $gi2 = $1;
			if ($gi1 == $gi2) {
				$score++;
			}
		}
	}
	my $final=0;
	if ($size>0) {
		$final=($score/$size);
	}
	return $final;
}


sub PrintToFileOrganizedClipboard {
	my $suffix = shift(@_);
	my $outFile = shift(@_);
	my $orgFile = shift(@_);
	
	$orgFile = $STORIdir . $orgFile;
	my $filePath = $STORIdir . $outFile;
	open clipFile, ">$filePath";
	#print Dumper \%clipboard_hash;
	#my $choice = shift(@_);

	my @famNames=();
	my @indices = keys %clipboard_hash;
	my @indices_s = sort { $a <=> $b } @indices;  #sort ascending
	foreach my $index (@indices_s) {
		push @famNames, $clipboard_hash{$index};
	}
	
	#print "\n";
	print clipFile "taxID name phylum ";
	foreach my $family (@famNames) {
		print clipFile $family . $suffix . " ";
	}
	print clipFile "\n";
	
	print "orgFile\= $orgFile\n";
	
	if (!($orgFile eq "-1")) {
		open Org, $orgFile;
	}
	
	my @taxids_ordered=();
	my @names_ordered=();
	my @phyla_ordered=();
	
	while (<Org>) {
		my $line = $_;
		chomp($line);
		if ($line =~ m/(\d+)\t(.+)\t(.+)/) {
			push @taxids_ordered, $1;
			push @names_ordered, $2;
			push @phyla_ordered, $3;
			#print "file\: $1 $2 $3\n";
		}
	}
	close Org;
	
	my $r=0;
	foreach my $taxon (@taxids_ordered) {
		print clipFile $taxon . " " . $names_ordered[$r] . " " . $phyla_ordered[$r] . " ";
		foreach my $family (@famNames) {	
			if (exists $clipboard{$family}{$taxon}) {
				print clipFile $clipboard{$family}{$taxon} . " ";
			}
			else {
				print clipFile "-1 ";
			}
		}
		print clipFile "\n";
		$r++;
	}
	#print "\n\n";
	close clipFile;
}


sub ShowOrganizedFilteredClipboard {
	my $orgFile = shift(@_);
	my $filter = shift(@_);
	$orgFile = $STORIdir . "/" . $orgFile;
	$filter = $STORIdir . "/" . $filter;
	#print Dumper \%clipboard_hash;
	#my $choice = shift(@_);

	my @famNames=();
	my @indices = keys %clipboard_hash;
	my @indices_s = sort { $a <=> $b } @indices;  #sort ascending
	foreach my $index (@indices_s) {
		push @famNames, $clipboard_hash{$index};
	}
	

	print "\n";
	print "taxID name phylum ";
	foreach my $family (@famNames) {
		print $clipboard{$family}{"score"} . " ";
	}
	print "\n";
	
	print "\n";
	print "taxID name phylum ";
	foreach my $family (@famNames) {
		print $family . " ";
	}
	print "\n";
	
	#print "orgFile\= $orgFile\n";
	
	if (!($orgFile eq "-1")) {
		open Org, $orgFile;
	}
	
	my @taxids_ordered=();
	#my @names_ordered=();
	#my @phyla_ordered=();
	
	my %taxaNomenclature=();
	
	while (<Org>) {
		my $line = $_;
		chomp($line);
		if ($line =~ m/(\d+)\t(.+)\t(.+)/) {
			#push @taxids_ordered, $1;
			#push @names_ordered, $2;
			#push @phyla_ordered, $3;
			
			$taxaNomenclature{$1}{"name"} = $2;
			$taxaNomenclature{$1}{"phylum"} = $3;
			
			#print "file\: $1 $2 $3\n";
		}
	}
	close Org;
	
	if (!($filter eq "-1")) {
		open Filter, $filter;
	}
	
	while (<Filter>) {
		my $line = $_;
		chomp($line);
		if ($line =~ m/(\d+)/) {
			push @taxids_ordered, $1;
		}
	}
	close Filter;
	
	my $r=0;
	foreach my $taxon (@taxids_ordered) {
		#print $taxon . " " . $names_ordered[$r] . " " . $phyla_ordered[$r] . " ";
		print $taxon . " " . $taxaNomenclature{$taxon}{"name"} . " " . $taxaNomenclature{$taxon}{"phylum"} . " ";
		foreach my $family (@famNames) {	
			if (exists $clipboard{$family}{$taxon}) {
				print $clipboard{$family}{$taxon} . " ";
			}
			else {
				print "-1 ";
			}
		}
		print "\n";
		$r++;
	}
	print "\n\n";
}


sub ShowOrganizedClipboard {
	my $orgFile = shift(@_);
	$orgFile = $STORIdir . "/" . $orgFile;
	#print Dumper \%clipboard_hash;
	#my $choice = shift(@_);

	my @famNames=();
	my @indices = keys %clipboard_hash;
	my @indices_s = sort { $a <=> $b } @indices;  #sort ascending
	foreach my $index (@indices_s) {
		push @famNames, $clipboard_hash{$index};
	}
	

	print "\n";
	print "taxID name phylum ";
	foreach my $family (@famNames) {
		print $clipboard{$family}{"score"} . " ";
	}
	print "\n";
	
	print "\n";
	print "taxID name phylum ";
	foreach my $family (@famNames) {
		print $family . " ";
	}
	print "\n";
	
	#print "orgFile\= $orgFile\n";
	
	if (!($orgFile eq "-1")) {
		open Org, $orgFile;
	}
	
	my @taxids_ordered=();
	my @names_ordered=();
	my @phyla_ordered=();
	
	while (<Org>) {
		my $line = $_;
		chomp($line);
		if ($line =~ m/(\d+)\t(.+)\t(.+)/) {
			push @taxids_ordered, $1;
			push @names_ordered, $2;
			push @phyla_ordered, $3;
			#print "file\: $1 $2 $3\n";
		}
	}
	close Org;
	
	my $r=0;
	foreach my $taxon (@taxids_ordered) {
		print $taxon . " " . $names_ordered[$r] . " " . $phyla_ordered[$r] . " ";
		foreach my $family (@famNames) {	
			if (exists $clipboard{$family}{$taxon}) {
				print $clipboard{$family}{$taxon} . " ";
			}
			else {
				print "-1 ";
			}
		}
		print "\n";
		$r++;
	}
	print "\n\n";
}


sub ShowFullClipboard {
	my $orgFile = shift(@_);
	$orgFile = $STORIdir . "/" . $orgFile;
	#print Dumper \%clipboard_hash;
	#my $choice = shift(@_);

	my @famNames=();
	my @indices = keys %clipboard_hash;
	my @indices_s = sort { $a <=> $b } @indices;  #sort ascending
	foreach my $index (@indices_s) {
		push @famNames, $clipboard_hash{$index};
	}
	
	print "\n";
	print "taxID name phylum ";
	foreach my $family (@famNames) {
		print $family . " ";
	}
	print "\n";
	
	#print "orgFile\= $orgFile\n";
	
	if (!($orgFile eq "-1")) {
		open Org, $orgFile;
	}
	
	my @taxids_ordered=();
	my @names_ordered=();
	my @phyla_ordered=();
	
	while (<Org>) {
		my $line = $_;
		chomp($line);
		if ($line =~ m/(\d+)\t(.+)\t(.+)/) {
			push @taxids_ordered, $1;
			push @names_ordered, $2;
			push @phyla_ordered, $3;
			#print "file\: $1 $2 $3\n";
		}
	}
	close Org;
	
	my $r=0;
	foreach my $taxon (@taxids_ordered) {
		my $fullFlag = 1;
		my $line = $taxon . " " . $names_ordered[$r] . " " . $phyla_ordered[$r] . " ";
		foreach my $family (@famNames) {	
			if (exists $clipboard{$family}{$taxon}) {
				$line = $line . $clipboard{$family}{$taxon} . " ";
			}
			else {
				$fullFlag = 0;
			}
		}
		if ($fullFlag == 1) {
			print $line . "\n";
		}
		$r++;
	}
	print "\n\n";
}


sub ShowClipboard {
	#print Dumper \%clipboard_hash;
	my $choice = shift(@_);
	if (defined $choice) {
		my $fam = $clipboard_hash{$choice};
		my %famHash = %{$clipboard{$fam}};
		print Dumper \%famHash;
	}
	else {
		my @indices = keys %clipboard_hash;
		my @indices_s = sort { $a <=> $b } @indices;  #sort ascending
		foreach my $index (@indices_s) {
			print $index . ": " . $clipboard_hash{$index} . "\n";
		}
	}
}


sub ShowRunResultOrganized {
	my $name = shift(@_);
	my $partialFile = shift(@_);
	my $orgFile = $STORIdir . $partialFile;
	
	my @s_counts = keys %{$runHash{$name}};
	my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
	my $sc = $sorted_counts[0];
		
	my $runNameA = $name . "a";
	my $runNameB = $name . "b";
	
	my $sourceFilesDirA = $run_params{$name}{"parentDir"} . "/" . $runNameA;
	my $sourceFilesDirB = $run_params{$name}{"parentDir"} . "/" . $runNameB;
	
	my $fileA = $sourceFilesDirA . "/STORI_out_sc" . $sc . "_" . $runNameA . ".txt";
	my $fileB = $sourceFilesDirB . "/STORI_out_sc" . $sc . "_" . $runNameB . ".txt";
	
	my $cmd = "perl " . $STORIdir . "checkSTORI.pl " . $fileA . " " . $fileB . " -org " . $orgFile;
	print "cmd is: $cmd\n";
	system($cmd);
}


sub SummarizeFiltered {
	my $name = shift(@_);
	my $orgFile = shift(@_);
	my @s_counts = keys %{$runHash{$name}};
	my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
	my $sc = $sorted_counts[0];
		my $runNameA = $name . "a";
		my $runNameB = $name . "b";
		
		my $sourceFilesDirA = $run_params{$name}{"parentDir"} . "/" . $runNameA;
		my $sourceFilesDirB = $run_params{$name}{"parentDir"} . "/" . $runNameB;
		
		my $fileA = $sourceFilesDirA . "/STORI_out_sc" . $sc . "_" . $runNameA . ".txt";
		my $fileB = $sourceFilesDirB . "/STORI_out_sc" . $sc . "_" . $runNameB . ".txt";
	
	my $cmd = "perl " . $STORIdir . "checkSTORI.pl " . $fileA . " " . $fileB . " -org " . $orgFile;
	my $output = qx($cmd);
	my @scores=();
	my @families=();
	my @lines = split /\n/, $output;
	my %clipboard_temp=();
	my %famScores=();
	foreach my $line (@lines) {
		if ($line =~ m/score\s(.+)/) {
			@scores = split / /, $1;
			print Dumper @scores;
			print "\n";
		}
		elsif ($line =~ m/phylum\s(.+)/) {
			@families = split / /, $1;
			my $r=0;
			foreach my $family (@families) { #family names can get long because of iterative name concatenations. in this block we retitle the family with a shorter unique name
				my $abbrev = $name . "_" . (substr $family, 0, 4) . (substr $family, -4) . "_" . $r;
				$families[$r] = $abbrev;
				my $score = shift(@scores);
				$famScores{$abbrev} = $score;
				$r++;
			}
			#print "families_abbrev: " . join("\n", @families) . "\n";
			print Dumper \%famScores;
		}
		elsif ($line =~ m/(\d+)\s\S+\s\S+\s(.+)/) {
			my $taxon = $1;
			my @gis = split / /, $2;
			my $x=0;
			foreach my $fam (@families) {
				my $temp = $gis[$x];
				if ((!($temp =~ m/-1/)) && (!($temp =~ m/\?/))) {
					$clipboard{$fam}{$taxon} = $temp;
					$clipboard_temp{$fam}{$taxon} = $temp; }
				$x++;
			}
		}
	}
	foreach my $fam (@families) {
		$clipboard{$fam}{"score"} = $famScores{$fam};
	}
	
	my @temp = keys %clipboard_hash;
	my $counter=0;
	if ($#temp != -1) {
		$counter = (Max(\@temp) + 1);
	}
	my $numAdded = 0;
	foreach my $fam (keys %clipboard_temp) {
		$clipboard_hash{$counter} = $fam;
		$clipboard_reverse_hash{$fam} = $counter;
		$counter++;
		$numAdded++;
	}
	
	print $numAdded . " families added to clipboard.\n";
}

sub Summarize {
	my $name = shift(@_);
	my @s_counts = keys %{$runHash{$name}};
	my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
	my $sc = $sorted_counts[0];
		my $runNameA = $name . "a";
		my $runNameB = $name . "b";
		
		my $sourceFilesDirA = $run_params{$name}{"parentDir"} . "/" . $runNameA;
		my $sourceFilesDirB = $run_params{$name}{"parentDir"} . "/" . $runNameB;
		
		my $fileA = $sourceFilesDirA . "/STORI_out_sc" . $sc . "_" . $runNameA . ".txt";
		my $fileB = $sourceFilesDirB . "/STORI_out_sc" . $sc . "_" . $runNameB . ".txt";
	
	my $cmd = "perl " . $STORIdir . "checkSTORI.pl " . $fileA . " " . $fileB;
	my $output = qx($cmd);
	my @families=();
	my @lines = split /\n/, $output;
	my %clipboard_temp=();
	my @scores=();
	my %famScores=();
	foreach my $line (@lines) {
		if ($line =~ m/score\s(.+)/) {
			@scores = split / /, $1;
			#print Dumper @scores;
			print "\n";
		}
		elsif ($line =~ m/family\s(.+)/) {
			@families = split / /, $1;
			my $r=0;
			foreach my $family (@families) { #family names can get long because of iterative name concatenations. in this block we retitle the family with a shorter unique name
				my $abbrev = $name . "_" . (substr $family, 0, 4) . (substr $family, -4) . "_" . $r;
				$families[$r] = $abbrev;
				my $score = shift(@scores);
				$famScores{$abbrev} = $score;
				$r++;
			}
			#print "families_abbrev: " . join(" ", @families) . "\n";
			#print Dumper \%famScores;
		}
		elsif ($line =~ m/(\d+)\s(.+)/) {
			my $taxon = $1;
			my @gis = split / /, $2;
			my $x=0;
			foreach my $fam (@families) {
				my $temp = $gis[$x];
				if ((!($temp =~ m/-1/)) && (!($temp =~ m/\?/))) {
					$clipboard_temp{$fam}{$taxon} = $temp;
					$clipboard{$fam}{$taxon} = $temp; }
				$x++;
			}
		}
	}
	foreach my $fam (@families) {
		$clipboard{$fam}{"score"} = $famScores{$fam};
	}
	
	my @temp = keys %clipboard_hash;
	my $counter=0;
	if ($#temp != -1) {
		$counter = (Max(\@temp) + 1);
	}
	my $numAdded=0;
	foreach my $fam (keys %clipboard_temp) {
		$clipboard_hash{$counter} = $fam;
		$clipboard_reverse_hash{$fam} = $counter;
		$counter++;
		$numAdded++;
	}
	
	print $numAdded . " families added to clipboard.\n";	
}

sub Scrunch {
	#meh this is just playing around
	my $input = shift(@_);
	my $target = shift(@_);
	my $numBlocks = ($target/3);
	if ($numBlocks =~ m/(\d+)\.\d+/) {
		$numBlocks = $1; }
		
	my @chars = split / */, $input;
	my @windows = ();
	for (my $x=0; $x<($#chars-2); $x++) {
		my $temp = $chars[$x] . $chars[($x+1)] . $chars[($x+2)];
		push @windows, $temp;
	}
	print Dumper \@windows;
	
	return "derp";
}



sub AssignNames {			#name all families containing more sequences than the cutoff
	my $cutoff=shift(@_);
	my %assignments=();
	my %conflictCheck=();
	
	foreach my $clipboardName (keys %clipboard) {
		my $ref = $clipboard{$clipboardName};
		if (GetFamSize($ref) > $cutoff) {
			my $name_ref = DoWordCensus($ref);
			my @temp = @{$name_ref};
			my $choice = $name_ref->[0]->[0];
			my $r=1;
			while ((exists $conflictCheck{$choice}) && ($r < $#temp)) { #moves on to the second-best name if a conflict exists
				$choice = $name_ref->[$r]->[0];
				$r++;
			}
			
			$choice =~ s/[^\w\d\s]//g;		#remove characters that may cause trouble later, eg ()/
			$assignments{$clipboardName} = $choice;
			if (!(exists $conflictCheck{$choice})) {
				my @temp = ($clipboardName);
				$conflictCheck{$choice} = [@temp]; }
			else {
				my @temp = @{$conflictCheck{$choice}};
				push @temp, $clipboardName;
				$conflictCheck{$choice} = [@temp];
			}
		}
		else {
			my @temp = ($clipboard_reverse_hash{$clipboardName});
			Clear(\@temp);
		}
	}
	
	foreach my $name (keys %conflictCheck) {
		my @temp = @{$conflictCheck{$name}};
		if ($#temp > 0) {
			print "naming conflict over $name\:\n";
			foreach my $code (@temp) {
				print "\(" . $clipboard_reverse_hash{$code} . "\) " . $code . "\n";
			}
			delete $conflictCheck{$name};
			#print "attempting to conflict resolution\n";
			#ResolveConflicts()
		}
	}
		
	foreach my $name (keys %conflictCheck) {      #name change block
		my @temp = @{$conflictCheck{$name}};
		my $oldName = $temp[0];
		my $newName = $name;
		my $oldIndex = $clipboard_reverse_hash{$oldName};
		
		$clipboard_hash{$oldIndex} = $newName;
		$clipboard_reverse_hash{$newName} = $oldIndex;
		delete $clipboard_reverse_hash{$oldName};
		
		my $ref = $clipboard{$oldName};
		Rename($ref, $oldName, $newName);
	}
	
	print Dumper \%assignments;
	print Dumper \%conflictCheck;
}

sub GetFamSize {
	my $ref=shift(@_);
	my %fam = %{$ref};
	my @temp = keys %fam;
	#my $ans = $#temp + 1;
	my $ans = $#temp; #extra taxon is the 'score'
	return $ans;
}


sub DoWordCensus {
	my %wordHash=();
	my $hash_ref = shift(@_);
	my @names = @{GetAnnotationBare($hash_ref)};
	foreach my $name (@names) {
		my @words = split / /, $name;
		for (my $x=0; $x <= ($#words - 2); $x++) {
			my $subword = $words[$x] . "_" . $words[$x+1] . "_" . $words[$x+2];
			if (!(exists $wordHash{$subword})) {
				$wordHash{$subword} = 1;
			}
			else {
				$wordHash{$subword}++;
			}
		}
	}
	
	my @nameArr=();
	foreach my $word (keys %wordHash) {
		my @temp = ($word, $wordHash{$word});
		push @nameArr, [@temp];
	}
	
	my @nameArr_s = sort { $b->[1] <=> $a->[1] } @nameArr;  #sort descending
	return \@nameArr_s;
	#print Dumper @nameArr_s;
}


sub GetAnnotationBare {
	my $hash_ref = shift(@_);
	my %txgi = %{$hash_ref};
	my @ansArr=();
	
	foreach my $taxon (keys %txgi) {
		my $gi = $txgi{$taxon};
		
		if (!($taxon eq "score")) {
			if ($taxon > 0) {
				my $annot = -1;
				if ($gi > 0) {
					my $cmd="$blastdbcmdPath -entry $gi -db $blastdbDir/$taxon -outfmt \"\%t\"";
					$annot = qx($cmd);
					chomp($annot);
				}
				push @ansArr, $annot;
			}
		}
	}
	return \@ansArr;
}


sub GetAnnotation {
	my $hash_ref = shift(@_);
	my %txgi = %{$hash_ref};

	foreach my $taxon (keys %txgi) {
		my $gi = $txgi{$taxon};
		
		if ($taxon > 0) {
			my $annot = -1;
			if ($gi > 0) {
				my $cmd="$blastdbcmdPath -entry $gi -db $blastdbDir/$taxon -outfmt \"\%t\"";
				$annot = qx($cmd);
				chomp($annot);
			}
			print $taxon . "\t" . $gi . "\t" . $annot . "\n";
		}
	}
}


sub ShowRunResult {
	my $name = shift(@_);
	my @s_counts = keys %{$runHash{$name}};
	my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
	my $sc = $sorted_counts[0];
	
	my $runNameA = $name . "a";
	my $runNameB = $name . "b";
	
	my $sourceFilesDirA = $run_params{$name}{"parentDir"} . "/" . $runNameA;
	my $sourceFilesDirB = $run_params{$name}{"parentDir"} . "/" . $runNameB;
	
	my $fileA = $sourceFilesDirA . "/STORI_out_sc" . $sc . "_" . $runNameA . ".txt";
	my $fileB = $sourceFilesDirB . "/STORI_out_sc" . $sc . "_" . $runNameB . ".txt";
	
	my $cmd = "perl " . $STORIdir . "checkSTORI.pl " . $fileA . " " . $fileB;
	#print "cmd is: $cmd\n";
	system($cmd);

}

sub ShowRunResultConcise {
	my $name = shift(@_);
	my @s_counts = keys %{$runHash{$name}};
	my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
	my $sc = $sorted_counts[0];
	
		my $runNameA = $name . "a";
		my $runNameB = $name . "b";
		
		my $sourceFilesDirA = $run_params{$name}{"parentDir"} . "/" . $runNameA;
		my $sourceFilesDirB = $run_params{$name}{"parentDir"} . "/" . $runNameB;
		
		my $fileA = $sourceFilesDirA . "/STORI_out_sc" . $sc . "_" . $runNameA . ".txt";
		my $fileB = $sourceFilesDirB . "/STORI_out_sc" . $sc . "_" . $runNameB . ".txt";
	
	my $cmd = "perl " . $STORIdir . "checkSTORI.pl " . $fileA . " " . $fileB . " -c";
	system($cmd);

}

sub GetID {
	my $name = shift(@_);
	my $ref = $groups{$name};
	my @runs = @{$ref};
	
	my @idlist=();
	
	foreach my $rname (@runs) {
		my @s_counts = keys %{$runHash{$rname}};
		my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
		my $sc = $sorted_counts[0];
		my $temp = $runHash{$rname}{$sc}{"ids"}->[0];
		push @idlist, $temp;
		$temp = $runHash{$rname}{$sc}{"ids"}->[1];
		push @idlist, $temp;
	}
	
	foreach my $id (@idlist) {
		print $id . "\n";
	}
	
}

sub ShowTheseRuns {
	my @runs = @{shift(@_)};
	
	foreach my $run (@runs) {
		my @avgScores =();
		my @s_counts = keys %{$runHash{$run}};
		my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
		
		foreach my $scnt (@sorted_counts) {
			if ($runHash{$run}{$scnt}{"scores"}->[0] != -1) {
				my $avgScore = Mean($runHash{$run}{$scnt}{"scores"});
				$avgScore = sprintf("%.2f", $avgScore);
				push @avgScores, $avgScore;
			}
			else {
				push @avgScores, -1;
			}
		}
		
		
		print $indexedNames_reverse{$run} . ": " . $run . " ";

			
		foreach my $avgScore (@avgScores) {
				if (defined $avgScore) {
					print $avgScore . " ";}
			}
		print "\n";
		
	}
}


sub SaveGroups {
	#print "saving groups\n";
	my $tmp = $STORIdir . "group_backup.txt";
	open grpFile, ">$tmp";
	my @grpNames=keys %groups;
	
	foreach my $name (@grpNames) {
		my @members=@{$groups{$name}};
		print grpFile $name . ": " . join(" ", @members) . "\n";
	}
	close grpFile;
}

sub LoadGroups {
	#print "loading groups\n";
	my $tmp = $STORIdir . "group_backup.txt";
	open grpFile, $tmp;
	while (my $line=<grpFile>) {
		chomp($line);
		$line=~ m/(.+):\s(.+)/;
		my $grpName = $1;
		my @members = split / /, $2;
		$groups{$grpName} = [@members];
	}
	close grpFile;
}


sub PairwiseRunCompare {
	my $run1 = shift(@_);
	my $run2 = shift(@_);
	
	my @STORIcounts = (keys %{$runHash{$run1}});
	my @STORIcounts_sorted = sort { $a <=> $b } @STORIcounts;  #sort ascending
	my $STORIcount = $STORIcounts_sorted[$#STORIcounts_sorted];
	my $id_1a = $runHash{$run1}{$STORIcount}{"ids"}->[0];
	my $id_1b = $runHash{$run1}{$STORIcount}{"ids"}->[1];

	@STORIcounts = (keys %{$runHash{$run2}});
	@STORIcounts_sorted = sort { $a <=> $b } @STORIcounts;  #sort ascending
	$STORIcount = $STORIcounts_sorted[$#STORIcounts_sorted];
	my $id_2a = $runHash{$run2}{$STORIcount}{"ids"}->[0];
	my $id_2b = $runHash{$run2}{$STORIcount}{"ids"}->[1];
		
	my $cmd = "perl " . $STORIdir . "checkSTORI.pl STORI-pbs.o" . $id_1a . " STORI-pbs.o" . $id_2a . " -q";
	#print "cmd is: $cmd\n";
	my $output = qx($cmd);
	$output =~ m/score\s(.+)\s$/;
	my @scoresA = split / /, $1;

	$cmd = "perl " . $STORIdir . "checkSTORI.pl STORI-pbs.o" . $id_1a . " STORI-pbs.o" . $id_2b . " -q";
	#print "cmd is: $cmd\n";
	$output = qx($cmd);
	$output =~ m/score\s(.+)\s$/;
	my @scoresB = split / /, $1;
	
	$cmd = "perl " . $STORIdir . "checkSTORI.pl STORI-pbs.o" . $id_1b . " STORI-pbs.o" . $id_2b . " -q";
	#print "cmd is: $cmd\n";
	$output = qx($cmd);
	$output =~ m/score\s(.+)\s$/;
	my @scoresC = split / /, $1;

	$cmd = "perl " . $STORIdir . "checkSTORI.pl STORI-pbs.o" . $id_1b . " STORI-pbs.o" . $id_2a . " -q";
	#print "cmd is: $cmd\n";
	$output = qx($cmd);
	$output =~ m/score\s(.+)\s$/;
	my @scoresD = split / /, $1;
	
	my @scoresCombined = (@scoresA, @scoresB, @scoresC, @scoresD);
	my $avg=Mean(\@scoresCombined);
	#print "avg agreement between $run1 and $run2: " . $avg . "\n";
	#print "all scorings: " . join(" ", @scoresCombined) . "\n";
	return $avg;
}

sub GroupStats {
	my $name = shift(@_);
	my $ref = $groups{$name};
	my @runs = @{$ref};
	
	my @times=(); my @numFams=(); my @scores=();
	
	foreach my $rname (@runs) {
		push @times, GetRunTime($rname);
		push @numFams, GetNumFamilies($rname);
		push @scores, GetRunScore($rname);
	}
	
	print "Stats for group   $name   [" . join("\, ", @runs) . "]\n";
	my $times_mean = Mean(\@times);
	my $times_sd = Sd(\@times);
	my $numFams_mean = Mean(\@numFams);
	my $numFams_sd = Sd(\@numFams);
	my $scores_mean = Mean(\@scores);
	my $scores_sd = Sd(\@scores);
	
	my @allscores=(); my %pairs=();
	foreach my $runA (@runs) {
		foreach my $runB (@runs) {
			my @tempArr= ($runA, $runB);
			my @tempSorted = sort { $a cmp $b } @tempArr;  #sort ascending
			if (!(exists $pairs{$tempSorted[0]}{$tempSorted[1]})) {
				$pairs{$tempSorted[0]}{$tempSorted[1]} =1;
				if (!($runA eq $runB)) {
					push (@allscores, PairwiseRunCompare($runA, $runB));
				}
			}		
		}
	}
	my $betweenRunScore_mean = Mean(\@allscores);
	my $betweenRunScore_sd = Sd(\@allscores);
	
	($times_mean, $times_sd, $numFams_mean, $numFams_sd, $scores_mean, $scores_sd, $betweenRunScore_mean, $betweenRunScore_sd) = @{Round($times_mean, $times_sd, $numFams_mean, $numFams_sd, $scores_mean, $scores_sd, $betweenRunScore_mean, $betweenRunScore_sd)};
	
	print "     Avg   Stdev\n";
	print "time  $times_mean  $times_sd\n";
	print "#fams  $numFams_mean  $numFams_sd\n";
	print "score  $scores_mean  $scores_sd\n";
	print "betweenSc   $betweenRunScore_mean  $betweenRunScore_sd\n";
}


sub Round {
	my @ans=();
	foreach my $num (@_) {
		push @ans, sprintf("%.2f", $num);
	}
	return \@ans;
}

sub RoundScalar {
	my $num=shift(@_);
	if ($num =~ m/\./) {
		$num = sprintf("%.2f", $num);
	}
	return $num;
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


sub MakeFilteredGroup {
	my $newName = shift(@_);
	my $groupToFilter = shift(@_);
	my $expr = shift(@_);
	
	my @preFiltered = @{$groups{$groupToFilter}};
	my @postFilter = ();
	
	foreach my $name (@preFiltered) {
		if ($name =~ m/($expr)/) {
			push @postFilter, $name;
		}
	}
	
	$groups{$newName} = [@postFilter];
}

sub Group {
	my $name = shift(@_);
	my $members_ref = shift(@_);
	my @members = @{$members_ref};
	
	my @final=();
	foreach my $mem (@members) {
		if (exists $index_master{$mem}) {
			push @final, $index_master{$mem};
		}
		else {
			push @final, $mem;
		}
	}
	
	$groups{$name} = [@final];
	#print "groups{$name} = [" . join(" ", @final) . "]\n";
}

sub GetRunTime {
	my $run = shift(@_);
	my @s_counts = keys %{$runHash{$run}};
	my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
	my $sc = $sorted_counts[0];
	my $total_rt=0;
	if ($runHash{$run}{$sc}{"runtime"} == -1) {
		shift(@sorted_counts);
	}
	foreach my $STORI (@sorted_counts) {
		$total_rt += $runHash{$run}{$STORI}{"runtime"};
	}
	return $total_rt;
}


sub IndexParser {
	my $input = shift(@_);
	my @noSpaces = @{SpaceParse($input)};
	my @comb1=(); my @comb2=(); my @comb3=();
	foreach my $noSpace (@noSpaces) {
		my @noCommas = @{CommaParse($noSpace)};
		foreach my $noComma (@noCommas) {
			my @noDashes = (); my @noDotDots=(); my @unproc=();
			if ($noComma =~ m/\-/) {
				@noDashes = @{DashParse($noComma)};
			}
			elsif ($noComma =~ m/\.\./) {
				@noDotDots = @{DotDotParse($noComma)};
			}
			else {
				@unproc = ($noComma);
			}
			@comb3 = (@noDashes, @noDotDots, @unproc);
			@comb2 = (@comb2, @comb3);
		}
		@comb1 = (@comb1, @comb2);
		@comb2=();
	}
	
	#print join(" ", @comb1) . "\n";
	return \@comb1;
	
	sub SpaceParse {
		my $input = shift(@_);
		my @spaceDelim = split / /, $input;
		return \@spaceDelim;
	}
	
	sub CommaParse {
		my $input = shift(@_);
		my @commaDelim = split /\,/, $input;
		#print "commaDelim\: " . join(" ", @commaDelim) . "\n";
		return \@commaDelim;
	}
	
	sub DashParse {
		my $input = shift(@_);
		my @list=();
		if ($input =~ m/(\d+)\-(\d+)/) {
			foreach my $x ($1..$2) {
				push @list, $x;
			}
		}
		return \@list;
	}
	
	sub DotDotParse {
		my $input = shift(@_);
		my @list=();
		if ($input =~ m/(\d+)\.\.(\d+)/) {
			foreach my $x ($1..$2) {
				push @list, $x;
			}
		}
		return \@list;
	}
}

sub GetNumFamilies {
	my $run = shift(@_);
	my @s_counts = keys %{$runHash{$run}};
	my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
	
	my $fams=0;
	
	my $sc = shift(@sorted_counts);
	my $scores_ref = $runHash{$run}{$sc}{"scores"};
	if ($scores_ref->[0] == -1) {
		$sc = shift(@sorted_counts);
		$scores_ref = $runHash{$run}{$sc}{"scores"};
	}
	if (defined $scores_ref) {
		$fams=Count($scores_ref);
	}
	
	sub Count {
		my $temp_ref = shift(@_);
		my @tempArr = @{$temp_ref};
		my $count=0;
		foreach my $c (@tempArr) {
			$count++;
		}
		return $count;
	}
	
	return $fams;
}

sub GetRunScore {
	my $run = shift(@_);
	my @s_counts = keys %{$runHash{$run}};
	my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
	
	my $score=0;
	
	my $sc = shift(@sorted_counts);
	my $scores_ref = $runHash{$run}{$sc}{"scores"};
	if ($scores_ref->[0] == -1) {
		$sc = shift(@sorted_counts);
		$scores_ref = $runHash{$run}{$sc}{"scores"};
	}
	if (defined $scores_ref) {
		$score=Mean($scores_ref);
	}
	return $score;
}


sub GetNumFamsTrend {
	my $run = shift(@_);
	my @s_counts = keys %{$runHash{$run}};
	my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
	
	my $num=0;
	my @famTrend=();
	
	my $sc = shift(@sorted_counts);
	my $scores_ref = $runHash{$run}{$sc}{"scores"};
	if ($scores_ref->[0] == -1) {
		$sc = shift(@sorted_counts);
		$scores_ref = $runHash{$run}{$sc}{"scores"};
	}
	if (defined $scores_ref) {
		$num=Count($scores_ref);
		push @famTrend, $num;
	}
	
	foreach $sc (@sorted_counts) {
		my $scores_ref = $runHash{$run}{$sc}{"scores"};
		if (defined $scores_ref) {
			$num=Count_t($scores_ref);
			push @famTrend, $num;
		}
	}
	
	sub Count_t { #exclude -1
		my($ref_arr) = @_;
		my @tempArr = @{$ref_arr};
		my @array;
		foreach my $elt (@tempArr) {
			if ($elt != -1) {
				push @array, $elt;
			}
		}
		my $answer = ($#array + 1);
		return $answer;
	}
	
	return \@famTrend;
}

sub GetAllRunScores {
	my $run = shift(@_);
	my @s_counts = keys %{$runHash{$run}};
	my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
	
	my $score=0;
	my @scores=();
	
	my $sc = shift(@sorted_counts);
	my $scores_ref = $runHash{$run}{$sc}{"scores"};
	if ($scores_ref->[0] == -1) {
		$sc = shift(@sorted_counts);
		$scores_ref = $runHash{$run}{$sc}{"scores"};
	}
	if (defined $scores_ref) {
		$score=Mean($scores_ref);
		push @scores, $score;
	}
	
	foreach $sc (@sorted_counts) {
		my $scores_ref = $runHash{$run}{$sc}{"scores"};
		if (defined $scores_ref) {
			$score=Mean($scores_ref);
			push @scores, $score;
		}
	}
	return \@scores;
}

sub ShowRuns {
	#print Dumper \%runHash;
	my @runNames = keys %runHash;
	my $numConvergences=0;
	my $numPaused=0;
	my %indexedNames=();
	my $count=0;
	foreach my $run (@runNames) {
		my @avgScores =();
		my @s_counts = keys %{$runHash{$run}};
		my @sorted_counts = sort { $b <=> $a } @s_counts;  #sort descending
		
		foreach my $scnt (@sorted_counts) {
			if ($runHash{$run}{$scnt}{"scores"}->[0] != -1) {
				my $avgScore = Mean($runHash{$run}{$scnt}{"scores"});
				$avgScore = sprintf("%.2f", $avgScore);
				push @avgScores, $avgScore;
			}
			else {
				push @avgScores, -1;
			}
		}		
		
		if ((!(exists $converged{$run})) && (!(exists $paused{$run}))) {
			print $count . ": " . $run . " ";
			
			foreach my $avgScore (@avgScores) {
				if (defined $avgScore) {
					print $avgScore . " ";}
			}
			print "\n";
		}
		$indexedNames{$count} = $run;
		$indexedNames_reverse{$run} = $count;
		$count++;
	}
	
	#print Dumper \%converged;
	my @temp = keys %converged;
	$numConvergences = ($#temp + 1);
	print "($numConvergences converged runs)\n";
	@temp = keys %paused;
	$numPaused = ($#temp + 1);
	print "($numPaused paused runs)\n";
	return \%indexedNames;
}



sub LoadRunData {		#loads a specified CSV file into a hash and returns the hash's memory reference
	my $file = shift(@_);
	my %runHash_pr=();
	my @runArr=();
	
	if (-e $file) {
		open (arrFile, $file);
		unless (flock(arrFile, 1)) {
			warn "File $file already locked; waiting...\n";
			flock(arrFile, 1) or die;
		}
		
		my @tempArr = <arrFile>;
		close arrFile;
		
		#print "begin load\n" . Dumper %runHash_pr;
		
		SET_LOOP: foreach my $line (@tempArr) {
			#print "parsing line: $line \n";
			if ($line =~ m/BEGIN\sPARAMETERS/) {
				last SET_LOOP;
			}
			chomp($line);
			my @row = split /\,/, $line;
			push @runArr, [@row];
		}
			
		
		if ($#runArr >= 0) {
			foreach my $row_ref (@runArr) {			#load the runArr, which is from the file, into a hash
				my @tempRow = @{$row_ref};
				my $runName = shift(@tempRow);
				my $runid_A = shift(@tempRow);
				my $runid_B = shift(@tempRow);
				my $score=""; 
				my $STORIcount=0;
				
				#print "tempRow= " . join(" ", @tempRow) . "\n";
				
				if ($runName =~ m/\_fin\_/) {
					my @fin_Arr = (); $fin_{$runName} = 1; @fin_Arr = keys %fin_; $groups{"finished"} = [@fin_Arr];
				}
				if ($runName =~ m/arch.*\_fin\_/) {
					my @arch_fin_Arr = (); $arch_fin_{$runName} = 1; @arch_fin_Arr = keys %arch_fin_; $groups{"arch_finished"} = [@arch_fin_Arr];
				}
				if ($runName =~ m/bact.*\_fin\_/) {
					my @bact_fin_Arr = (); $bact_fin_{$runName} = 1; @bact_fin_Arr = keys %bact_fin_; $groups{"bact_finished"} = [@bact_fin_Arr];
				}
				if ($runName =~ m/euk.*\_fin\_/) {
					my @euk_fin_Arr = (); $euk_fin_{$runName} = 1; @euk_fin_Arr = keys %euk_fin_; $groups{"euk_finished"} = [@euk_fin_Arr];
				}
				
				while (@tempRow) {
					#print join(" ", @tempRow) . "\n";
					my @scores=();
					while (!($score =~ m/r/)) {
						if ($score =~ m/\d+/) {
							push @scores, $score; }
						$score = shift(@tempRow);
					}
					my $temp="scores";
					$runHash_pr{$runName}{$STORIcount}{$temp} = [@scores];
					
					chomp($score);
					$score =~ m/r(.+)/;
					my $time=$1;
					$temp="runtime";
					$runHash_pr{$runName}{$STORIcount}{$temp} = $time;
					
					$score = shift(@tempRow);
					if (defined $score) {
						if ($score =~ m/converged/) {
							my @convergedArr = ();							
							$converged{$runName} = 1;	
							@convergedArr = keys %converged;
							$groups{"converged"} = [@convergedArr];
						}
						elsif ($score =~ m/paused/) {
							my @pausedArr = ();							
							$paused{$runName} = 1;	
							@pausedArr = keys %paused;
							$groups{"paused"} = [@pausedArr];
						}
						#print Dumper \%converged;
					}
					$STORIcount++;		
				}	
				$STORIcount--;
				my $temp = "ids"; my @tempArr = ($runid_A, $runid_B);
				$runHash_pr{$runName}{$STORIcount}{$temp} = [@tempArr];
			}
			foreach my $runName (keys %runHash_pr) {
				if ((!(exists $converged{$runName})) && (!(exists $paused{$runName}))) {
					my @unconvergedArr = ();
					$unconverged{$runName} = 1;
					@unconvergedArr = keys %unconverged;
					$groups{"unconverged"} = [@unconvergedArr];
				}
			}
			my @allRuns = keys %runHash_pr;
			$groups{"all"} = [@allRuns];
		}
	}
	#print "here is groups in LoadRunData\n";
	#print Dumper \%groups;
	#print "end load\n" . Dumper %runHash_pr;
	return \%runHash_pr;
}



sub LoadParams {

	#print "(LoadParams runHash)\n";
	#print Dumper \%runHash;
	
	my %run_params=();
	my $file = shift(@_);
	my $readFlag=0;
	if (-e $file) {
		open (arrFile, $file);
		unless (flock(arrFile, 1)) {
			warn "File $file already locked; waiting...\n";
			flock(arrFile, 1) or die;
		}
		
		my @tempArr = <arrFile>;
		my @tempParamsArr=();
		close arrFile;
		my $r =0;
		foreach my $line (@tempArr) {
			chomp($line);
			if ($line =~ m/BEGIN\sPARAMETERS/) {
				$readFlag=1;
			}
			if ($readFlag==1) {
				my @row = split /\,/, $line;
				push @tempParamsArr, [@row];
				$r++;
			}
		}
		
		foreach my $ref (@tempParamsArr) {
			my @tempRow = @{$ref};
			my $runName = shift(@tempRow);
			foreach my $elt (@tempRow) {
				$run_params{$runName}{"parentDir"} = shift(@tempRow);
				$run_params{$runName}{"taxaFile"} = shift(@tempRow);
				$run_params{$runName}{"windowSize"} = shift(@tempRow);
				$run_params{$runName}{"finalMaxFams"} = shift(@tempRow);
			}
		}
	}
	return \%run_params;
}




sub Max {
	my($ref_arr) = @_;
	my @array = @{$ref_arr};
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	my $answer = $stat->max();
	return $answer;
}

sub Print2Darr {	
	my($ref_arr) = @_;
	
	#print "entered Print2Darr. \n ref_arr rows = $rows and cols = $cols \n";
	
	my @array = @{$ref_arr};
	
	print "\n";
	for (my $r=0; $r<= $#array; $r++) {
		for (my $c=0; $c<= $#{$array[$r]}; $c++) {
			print "$array[$r][$c] ";
		}
		print "\n";
	}
	print "\n";
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