use lib '/home/josh/perl5/lib';
use lib '/home/josh/perl5/lib/perl5';

use Statistics::Descriptive;

#takes fasta files as input, makes them nonredundant, and then makes them into blast databases
#also outputs simulated results for what would happen if each taxon were BLASTed against itself

$homeDir = "/home/josh/";
$projectDir = $homeDir . "scratch/archaea"; 
$hitDir = $projectDir . "/hits";
#get taxa list
$taxaFile = $homeDir . "STORI/taxids_GIs_archS.txt";   #a list of taxids

$blastdbcmdPath = $homeDir . "STORI/blastdbcmd";
$makeblastdbPath= $homeDir . "STORI/makeblastdb";  #path to the blast+ program makeblastdb
$blastDirPath = $projectDir . "/blast";
$nrdbLoc = $homeDir . "STORI/bp_nrdb_SHA.pl";


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
	
	$command="perl $nrdbLoc -o $nrfilename $filename";
	print "cmd is: $command\n";
	system("$command");
	$command="$makeblastdbPath -dbtype prot -parse_seqids -in $nrfilename -title $taxon -out $blastDirPath/$taxon";
	print "cmd is: $command\n";
	system("$command");
	$command = "$blastdbcmdPath -entry all -db $blastDirPath\/$taxon -outfmt \"\%a \%t\"";
	my $annot = qx($command);
	my @annotArr = split /\n/, $annot;

	my $pathHitfile = $hitDir . "/" . $taxon . "_" . $taxon;
        open (hitfile, ">$pathHitfile");

	while (@annotArr) {
                my $subject = shift(@annotArr);
                $subject =~ m/^(.+?)\s(.+)$/;
                my $accession = $1;
		print hitfile $accession . "\t" . $accession . "\t100\n";
        }
	close hitfile;
}
