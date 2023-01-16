use lib '/nv/hp10/jstern7/perl5reinstall/lib';
use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';

use Statistics::Descriptive;

#takes fasta files as input, makes them nonredundant, and then makes them into blast databases

$projectDir = "/home/josh/scratch/viruses_2017"; 
$hitDir = $projectDir . "/hits";

#get taxa list
$taxaFile = "/home/josh/STORI/taxids_GIs.txt";   #a list of taxids

$makeblastdbPath="/home/josh/STORI/makeblastdb";  #path to the blast+ program makeblastdb
$blastDirPath = $projectDir . "/blast";

$nrdbLoc = "/home/josh/STORI/bp_nrdb_SHA.pl";


open (taxa, $taxaFile);
@fileArr=<taxa>;
$count=0;
my @taxaArr=();
foreach my $line (@fileArr) {
	chomp($line); 
	if ($line =~ m/^(\d+)\s.+/) {
		push @taxaArr, $1;
	}
}
	
foreach $taxon (@taxaArr) {
	$filename = $projectDir . "/fasta/" . $taxon . ".fasta";
	$nrfilename = $projectDir . "/fasta/nr" . $taxon . ".fasta";
	$nrGIfilename = $projectDir . "/fasta/nrgi" . $taxon . ".fasta";
	
	$command="perl $nrdbLoc -o $nrfilename $filename";
	print "cmd is: $command\n";
	system("$command");
	RemoveExtraGIs($nrfilename, $nrGIfilename, $taxon);
	$command="$makeblastdbPath -parse_seqids -in $nrGIfilename -title $taxon -out $blastDirPath/$taxon";
	print "cmd is: $command\n";
	system("$command");
        #print "ok lets take 5";
        #sleep 5;
        #print "\n";
        #sleep 300;
}


sub RemoveExtraGIs {  #doesn't actually remove anything. NCBI no longer uses GIs in a way that is compatible with STORI
		      #furthermore, even if a sequence were to have multiple accessions in our custom BLAST databases,
		      #I have a hunch it won't matter and the sequence will just be retrievable via any of its accessions
		      #Keeping the name out of laziness.
		      #All this does is create the special-case hitfiles for what the result should be when we blast a 
		      #taxon against itself.
		      #The only scenario in which having multiple accessions per sequence might be a problem is if
		      #the BLAST results return all of the accessions rather than just 1 of them,
		      #and in a subsequent search STORI interprets all of the accessions as a single accession
		      #this will need testing to clear up.
	my $nrFileName = shift(@_);
	my $nrGIfilename = shift(@_);
	my $taxID = shift(@_);
	
	open (nrfile, $nrFileName);
	open (nrgifile, ">$nrGIfilename");
	
	my $pathHitfile = $hitDir . "/" . $taxID . "_" . $taxID;
	open (hitfile, ">$pathHitfile");
	
	@fastaArr=<nrfile>;

	my $r=0;
	foreach $line (@fastaArr)
	{
		if ($line =~ m/^>(.+?)\s(.+)$/) { 		#match before space nongreedily
			my $accession = $1;
			print hitfile $accession . "\t" . $accession . "\t100\n";
		}
		chomp($line);
		$fastaArr[$r] = $line;
	
		print nrgifile "$fastaArr[$r]\n";
		
		$r++;
	}

	close nrfile;
	close nrgifile;
	close hitfile;
}
