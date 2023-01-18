#STORI begins by accepting from the
#user a set of seed protein sequences and an upper limit to the number of
#different families (sets of orthologs) it may retrieve. STORI randomly
#scatters the seeds into two parallel, independent iterator processes. 

#command: perl /path/to/beginSTORI.pl runNumber /path/to/dir/containing/source/dir taxaFile windowSize finalMaxFams

use Data::Dumper;

my $home = "/home/ec2-user/STORI";
my $blastdb_dir = $home . "/universal20150110/blast";
my $blastdbcmdPath = $home . "/blastdbcmd";
my $runNumber=shift(@ARGV);
my $parentDirPath=shift(@ARGV);
my $searchOrderFile=shift(@ARGV);
my $windowSize=shift(@ARGV);
my $finalMaxFams=shift(@ARGV);

if ((!(defined $runNumber)) || (!(defined $parentDirPath)) || (!(defined $searchOrderFile)) || (!(defined $windowSize)) || (!(defined $finalMaxFams))) {
	print "usage: perl /path/to/beginSTORI.pl runNumber /path/to/dir/containing/source/dir taxaFile windowSize finalMaxFams\n";
	exit;
}

my $sourceDirPathA = $parentDirPath . "/" . $runNumber . "a";
my $sourceDirPathB = $parentDirPath . "/" . $runNumber . "b";


my $queryGiFileA = "query-gis\[" . $runNumber . "a-iter0\].txt";
my $queryGiFileB = "query-gis\[" . $runNumber . "b-iter0\].txt";

my $runA_id;
my $runB_id;

my $runIDfile = $parentDirPath . "/tempRunIDfile.txt";
my $cmd = "touch $runIDfile";
system($cmd);

open beginTemp, ">$runIDfile";
print beginTemp "runA -1\;\nrunB -1\;\n";
close beginTemp;



my $go="yes";
if ((-d $sourceDirPathA) || (-d $sourceDirPathB)) {
	$go="no";
	print "a sourceDirPath already exists. enter yes if you want to continue & overwrite.\n";
	$go=<>;
	
	if ($go =~ m/yes/) {
		my $cmd = "rm -r $sourceDirPathA";
		print "cmd is: " . $cmd . "\n";
		system($cmd);
		my $cmd = "rm -r $sourceDirPathB";
		print "cmd is: " . $cmd . "\n";
		system($cmd);
	}
}
if ($go =~ m/yes/) {
	my $cmd = "mkdir " . $sourceDirPathA;
	print "cmd is: " . $cmd . "\n";
	system($cmd);
	my $cmd = "mkdir " . $sourceDirPathB;
	print "cmd is: " . $cmd . "\n";
	system($cmd);
	
	my $orderAfile = $sourceDirPathA . "/taxa-master[". $runNumber . "a].txt";
	my $orderBfile = $sourceDirPathB . "/taxa-master[". $runNumber . "b].txt";
	
	my $cmd = "cp $home/taxa-master[" . $searchOrderFile . "].txt " .  $orderAfile;
	print "cmd is: " . $cmd . "\n";
	system($cmd);
	my $cmd = "cp $home/taxa-master[" . $searchOrderFile . "].txt " .  $orderBfile;
	print "cmd is: " . $cmd . "\n";
	system($cmd);
	
	my $queryGiPathA = $sourceDirPathA . "/". $queryGiFileA;
	my $queryGiPathB = $sourceDirPathB . "/". $queryGiFileB;
	
	my %annot_hash = ();
	my @taxaArr = @{MakeTaxaArr($orderAfile)};
	
	my $satisFlag="no"; my $annotQuery="derp";
	while (!($satisFlag =~ m/yes/)) {
		my $cmd = "rm " . $queryGiPathA;
		print "cmd is: " . $cmd . "\n";
		system($cmd);
		my $cmd = "rm " . $queryGiPathB;
		print "cmd is: " . $cmd . "\n";
		system($cmd);
		$annotQuery="derp";
		print "\nPlease enter an expression to match with protein names: ";
		$annotQuery=<>;
		print "you entered " . $annotQuery . "\n";
		
		print "what offset factor\? (usually 3)\n";
		my $offset_factor = <>;
		
		my $seed1 = MakeQueryGis_new(\@taxaArr, $queryGiPathB, $annotQuery, $offset_factor);
		my $seed2 = MakeQueryGis_new(\@taxaArr, $queryGiPathA, $annotQuery, $offset_factor);
		
		my $r = $seed1 / $seed2;
		
		if ( ($r > 3) || ($r < 0.33) ) {
			print "warning- larger than expected discrepency between number of seed gis. consider repeating the search.\n\n";
		}
		if (($seed1 < ($finalMaxFams/2)) || ($seed2 < ($finalMaxFams/2))) {
			print "warning- one of the seed files is smaller than expected given the number of max families. consider reducing the number of families you expect\n\n";
		}
		
		print "satisfied?\n";
		$satisFlag=<>;
	}
	
	print "\nI have created suggestions for the query-gi files. Please adjust them now, and type \"blastoff\" when you have finished and are ready for the STORI to begin.\nThe files are located here:\n$sourceDirPathA\n$sourceDirPathB\n";
	my $input;
	while (!($input eq "blastoff")) {
		print "3 2 1\>";
		$input=<>;
		chomp($input);
	}
	
	if ($input eq "blastoff") {
		
		my $command = "exec perl $home/STORI.pl $runNumber" . "a $sourceDirPathA $windowSize $finalMaxFams 0 > $sourceDirPathA" . "/STORI-log-sc0.txt 2>&1";
		
		$runA_id = fork();

		die "unable to fork: $!" unless defined($runA_id);
		if (!$runA_id) {  # child
			#print "i am the child proc my id is " . $runA_id . "\n";
			#print "cmd is: " . $command . "\n";
			exec($command);
			die "unable to exec: $!";
			exit;
		}
		print "runA_id is " . $runA_id . "\n";
		
		$command = "exec perl $home/STORI.pl $runNumber" . "b $sourceDirPathB $windowSize $finalMaxFams 0 > $sourceDirPathB" . "/STORI-log-sc0.txt 2>&1";
		
		$runB_id = fork();
		
		die "unable to fork: $!" unless defined($runB_id);
		if (!$runB_id) {  # child
			#print "i am the child proc my id is " . $runB_id . "\n";
			#print "cmd is: " . $command . "\n";
			exec($command);
			die "unable to exec: $!";
			exit;
		}
		print "runB_id is " . $runB_id . "\n";
		print "liftoff!\n";
	}
}

if ( ($runA_id > 0) && ($runB_id > 0) ) {
	open beginTemp, ">$runIDfile";
	print beginTemp "runA $runA_id\;\nrunB $runB_id\;\n";
	close beginTemp;
}

sub MakeRandomSearchOrder {
	my $taxaFile = shift(@_);
	open taxa, $taxaFile;
	my @taxaArr=<taxa>;
	close taxa;
	
	my $count=0;
	foreach my $taxon (@taxaArr)
	{
		chomp($taxon);  #removes any newline characters at the end of the string
		$taxaArr[$count] = $taxon;
		$count++;
	}
	
	fisher_yates_shuffle(\@taxaArr);
	
	open taxa, ">$taxaFile";
	foreach my $taxon (@taxaArr)
	{
		print taxa $taxon . "\n";
	}
	close taxa;
	
	return \@taxaArr;
}

sub MakeTaxaArr {
	my $taxaFile = shift(@_);
	open taxa, $taxaFile;
	my @taxaArr=<taxa>;
	close taxa;
	
	my $count=0;
	foreach my $taxon (@taxaArr)
	{
		chomp($taxon);  #removes any newline characters at the end of the string
		$taxaArr[$count] = $taxon;
		$count++;
	}
	fisher_yates_shuffle(\@taxaArr);
	return \@taxaArr;
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


sub MakeQueryGis_new {
	%annot_hash = ();
	my $arr_ref = shift(@_);
	my $queryGiPath = shift(@_);
	my $queryText = shift(@_);
	my $offset_factor = shift(@_);
	
	chomp($queryText);
	my @smallTaxaArr=@{$arr_ref};
	print "searching taxID ";
	foreach my $taxon (@smallTaxaArr) {			
		print $taxon . " ";
		my $cmd="$blastdbcmdPath -entry all -db $blastdb_dir\/$taxon -outfmt \"\%a \%t\"";
		my $annot = qx($cmd);
		my @annotArr = split /\n/, $annot;
		
		while (@annotArr) {
			my $subject = shift(@annotArr);
			if ($subject =~ m/($queryText)/) {
				$subject =~ m/^(.+?)\s(.+)$/;
				my $gi = $1;					#not actually the GI number anymore, but leaving variable name as is
				my $annotation = $2;
				$annot_hash{$taxon}{$gi} = $annotation;
			}
		}
		
	}
	
	@smallTaxaArr = keys %annot_hash;  #we only want to choose from taxa that contain hits
	my $taxaSize = ($#smallTaxaArr + 1);
	
	print "\nsmallTaxaArr\: $#smallTaxaArr $taxaSize\n";
	fisher_yates_shuffle(\@smallTaxaArr);

	my @taxaSample = ();
	
	my @seedGis=();
	while (($#seedGis < $finalMaxFams) && ($#smallTaxaArr > -1)) {
		my $randomTaxon = shift(@smallTaxaArr);
		push @taxaSample, $randomTaxon;
		my %temp = %{$annot_hash{$randomTaxon}};
		foreach my $gi (keys %temp) {
			push @seedGis, $gi;
		}
	}
	print "taxa sampled\: " . join(" ", @taxaSample) . "\n";
	
	open giFile, ">>$queryGiPath";

#	The below foreach loop ensures that when the STORI.pl iterator begins, 
#	each seed has its own unique
#“quasi”-family. The quasi-families are not biologically meaningful families.
#However, as iteration progresses in STORI.pl, these quasi-families will become plausible
#predictions of families (orthologous proteins). 
#	
#		To prevent abundant families from outcompeting sparse families, beginSTORI.pl can boost
#		the seed sequence scores. Using a
#		single parameter $offset_factor, the user defines a range for all initial seed scores,
#		specifying the extent to which the seeds persist for multiple taxa list
#		traversals within the first pair of STORI.pl instances.

	my $x=0;
	foreach my $gi (@seedGis) {
		my $randScore = GenerateRandomScore($offset_factor);
		print giFile $gi . " " . $randScore . "\n";
		print giFile "h fam" . $x . "\n\n";
		$x++;
	}
	
	close giFile;
	
	my $ans = $#seedGis;
	print "Number of seed GIs\: $ans\n\n";
	return $ans;
}

sub GenerateRandomScore {
	my $ans;
	
	my $offset_factor = shift(@_);
	my $offset_raw = $offset_factor * ($windowSize * $windowSize);
	
	my $min = (($windowSize * $windowSize) + $offset_raw);
	my $max = (1.1 * $min);
	my $range = ($max - $min);
	my $ans = (int rand ($range)) + $min;
	
	return $ans;
}
