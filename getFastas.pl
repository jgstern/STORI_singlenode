# Retrieves all protein sequences for the genomes specified in $taxaFile
# by querying the NCBI nuccore database and its links to the protein database.
# The sequences are saved as FASTA files, one file per taxon, in $projectDir
# (Actually, the accessions in $genomeFile need not be genomes - 
# all that is required is for the Nucleotide GI to be a valid record, and to
# have links to the Protein database.)
# the format for $taxaFile:
# <taxon id> <nucleotide GI 1> <nucleotide GI 2> [...] [<nucleotide GI n>]
# <taxon id> na
# If the genome for a particular taxon is
# incomplete, you can type "na" next to the taxID, and the program will retrieve
# protein sequences directly from the nr protein database. The disadvantage
# of retrieval directly from nr is that there may be organellar sequences.
# Using the links to nr from RefSeq 
# chromosome nucleotide sequences is the preferred method of retrieval,
# if the genome is complete.


#use lib '/nv/hp10/jstern7/perl5reinstall/lib';
#use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';
use LWP::Simple;
use Data::Dumper;

$projectDir = "/home/ec2-user/STORI/universal20150110";  
#get taxa list
$taxaFile = "/home/ec2-user/STORI/taxids_GIs.txt";   #a list of taxids and Genome GIs


if (!(-d $projectDir)) {
	my $cmd = "mkdir $projectDir";
	system($cmd);
}

$blastDirPath = $projectDir . "/blast";
if (!(-d $blastDirPath)) {
	my $cmd = "mkdir $blastDirPath";
	system($cmd);
}

$fastaDirPath = $projectDir . "/fasta";
if (!(-d $fastaDirPath)) {
	my $cmd = "mkdir $fastaDirPath";
	system($cmd);
}

my %accessionHash=();
open (TAXA, $taxaFile);
@fileArr=<TAXA>;
close TAXA;
$count=0;
foreach my $line (@fileArr) {
	my @lineArr=split / /, $line;
	my $taxID = shift(@lineArr);
	my @temp = ();
	SETLOOP: foreach my $elt (@lineArr) {
		if ($elt =~ m/^(\d+)$/) {
			push @temp, $1;
		}
		elsif ($elt =~ m/^na$/) {
			$accessionHash{$taxID} = "na";
			last SETLOOP;
		}
	}
	$accessionHash{$taxID} = [@temp] if ($#temp >= 0);
 }

#print Dumper \%accessionHash;

my %retrievalReport=();

foreach my $taxon (keys %accessionHash) {
	my $filename = $fastaDirPath . "/" . $taxon . ".fasta";
	
	my $outcome;
	
	if ($accessionHash{$taxon} eq "na") {
		$outcome = RetrieveFasta_nr_nofilter($taxon, $filename);
	}
	else {
		$outcome = RetrieveFasta_linked($taxon, $filename);
	}
	
	$retrievalReport{$taxon} = $outcome;
}

print Dumper \%retrievalReport;


sub RetrieveFasta_linked {
	my $taxon = shift(@_);
	my $filename = shift(@_);
	
	my @giArr = @{$accessionHash{$taxon}};
	
	my $total=0;
	my $successful=0;
	my $efetch_out="";
	foreach my $id (@giArr) {
	
		# Download protein records linked to nucleotide records.
		my $db1 = 'nuccore';  # &dbfrom
		my $db2 = 'protein';     # &db
		my $linkname = 'nuccore_protein'; # desired link &linkname
		#input UIDs in $db1 (protein GIs)
		
		#assemble the elink URL
		$base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
		$url = $base . "elink.fcgi?dbfrom=$db1&db=$db2&id=$id";
		$url .= "&linkname=$linkname&cmd=neighbor_history";
	
		#post the elink URL
		print "$url\n";
		$output = get($url);
		
		#parse WebEnv and QueryKey
		$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
		$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
		
		### include this code for ELink-EFetch
		#assemble the efetch URL
		$url = $base . "efetch.fcgi?db=$db2&query_key=$key&WebEnv=$web";
		$url .= "&rettype=fasta&retmode=text";
		
		my $fetched = get($url);
		
		if (!($fetched =~ m/\<ERROR\>/)) {
			#open output file for writing
			$efetch_out .= $fetched;
			$efetch_out .= "\n";
			print "$url\n";
			$successful++;
		}
		$total++;
		sleep 5;
	}
	
	open(OUT, ">$filename") || die "cant open file\n";
	
	print OUT "$efetch_out";
	
	close OUT;
	
	$ans = $successful . "/" . $total;
	return $ans;
}


sub RetrieveFasta_nr {
	my $taxon = shift(@_);
	my $filename = shift(@_);
	
	my $successes = 0;
	my $total = 0;
	
	$query = "(txid" . $taxon . "[orgn])+NOT+((chloroplast*[Protein+Name])+OR+(mitochondrial*[Protein+Name]))";

	#assemble the esearch url
	$base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
	$url = $base . "esearch.fcgi?db=protein&term=$query&usehistory=y";

	#post the esearch url
	$output = get($url);
	print "$url\n";

	#parse WebEnv, QueryKey and Count (# records retreived)
	$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
	$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
	$count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);
		#open output file for writing
	open(OUT, ">$filename") || die "cant open file\n";
	
	#retrieve data in batches of 500
	$retmax = 500;
	for ($retstart = 0; $retstart < $count; $retstart += $retmax) {
		$efetch_url = $base . "efetch.fcgi?db=protein&WebEnv=$web";
		$efetch_url .= "&query_key=$key&retstart=$retstart";
		$efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
		$efetch_out = get($efetch_url);
		if (!($efetch_out =~ m/.*ERROR.*/)) {
			print OUT "$efetch_out";
			#print "done\n";
			$successes++;
		}
		else {
			#print "error retrieving $url\n $efetch_url\n file not written\n\n";
		}
		$total++;
		sleep 15;
	}
	#print "done\n";
	close OUT;
	my $ans = $successes . "/" . $total;
	return $ans;
}

sub RetrieveFasta_nr_nofilter {
	my $taxon = shift(@_);
	my $filename = shift(@_);
	
	my $successes = 0;
	my $total = 0;
	
	$query = "(txid" . $taxon . "[orgn])";

	#assemble the esearch url
	$base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
	$url = $base . "esearch.fcgi?db=protein&term=$query&usehistory=y";

	#post the esearch url
	$output = get($url);
	print "$url\n";

	#parse WebEnv, QueryKey and Count (# records retreived)
	$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
	$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
	$count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);
		#open output file for writing
	open(OUT, ">$filename") || die "cant open file\n";
	
	#retrieve data in batches of 1000
	$retmax = 1000;
	for ($retstart = 0; $retstart < $count; $retstart += $retmax) {
		$efetch_url = $base . "efetch.fcgi?db=protein&WebEnv=$web";
		$efetch_url .= "&query_key=$key&retstart=$retstart";
		$efetch_url .= "&retmax=$retmax&rettype=fasta&retmode=text";
		$efetch_out = get($efetch_url);
		if (!($efetch_out =~ m/.*ERROR.*/)) {
			print OUT "$efetch_out";
			#print "done\n";
			$successes++;
		}
		else {
			#print "error retrieving $url\n $efetch_url\n file not written\n\n";
		}
		$total++;
		sleep 10;
	}
	#print "done\n";
	close OUT;
	my $ans = $successes . "/" . $total;
	return $ans;
}