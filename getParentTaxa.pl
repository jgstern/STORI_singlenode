#!/usr/bin/perl -w
use lib '/nv/hp10/jstern7/perl5reinstall/lib';
use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';
use strict;
use Fcntl ':flock';

my $qgiFile = shift(@ARGV);
my $hitDir = shift(@ARGV);
my $taxaFile = shift(@ARGV);
my $sourceFilesDir = shift(@ARGV);
my $parentTaxaDataPath = $sourceFilesDir . "/parentTaxaData.txt";


#1: build gi_lookup_hash{$gi} = $taxon; load @taxaArr
#get taxa list
my $cmd="touch $taxaFile";
system($cmd);
open (taxa, $taxaFile);
my @taxaArr=<taxa>;

my $count=0;
foreach my $taxon (@taxaArr)
{
	chomp($taxon);  #removes any newline characters at the end of the string
	$taxaArr[$count] = $taxon;
	$count++;
}

#use the hitfiles for a taxon blasted against itself to make a hash of every gi with its parent taxon
my %gi_lookup_hash=();
my @giArr=();
foreach my $taxon (@taxaArr)
{
	my $targetName = $taxon . "_" . $taxon;
	open (hitFile, "$hitDir/$targetName");
	unless (flock(hitFile, 1)) {
		warn "File $hitDir/$targetName already locked; waiting...\n";
		flock(hitFile, 1) or die;
	}
	@giArr = <hitFile>;
	close hitFile;
	$count=0;
	foreach my $gi (@giArr) #we just want the query gis, not the hits
	{
		chomp($gi);
		$gi =~ s/^(.+?)\s+.+?\s+\d+.*$/$1/;
		$giArr[$count] = $gi;
		#print $giArr[$count] . "\n"; 
		$count++;
	}
	
	my @nrGiArr=(); my $tempGi=-1;  #filter out consecutive duplicate GIs
	while ($#giArr>=0)
	{
		$tempGi=shift(@giArr);
		#print "comparing $giArr[0] to $tempGi\n";
		if (defined $giArr[0]) { 
			if ($giArr[0] ne $tempGi)
			{
				push @nrGiArr, $tempGi;
				#print "added $tempGi to nrGiArr\n";
			}
		}
		else {
			push @nrGiArr, $tempGi;
			#print "added $tempGi to nrGiArr\n";
		}
	}
	
	foreach my $gi (@nrGiArr)  #add gi & its taxon to gi_lookup_hash
	{
		$gi_lookup_hash{$gi} = $taxon;  # this should add (key, value) to gi_lookup_hash where key is gi & value is parent taxon
		#print "adding $gi of taxon $taxon\n";
	}
}


#2: assign query gis to their parent taxa
open (BATCH_QUERY, "$qgiFile");
my @queryArr_unparsed = <BATCH_QUERY>; my %tempHash=();
close BATCH_QUERY;

foreach my $line (@queryArr_unparsed) {		
	chomp($line);
	if (!($line =~ m/^h\s.+$/)) {
		if ($line =~ m/^(.+?)\s(\d+)$/) {
			#tempHash{<gi>} = <parent taxon>
			$tempHash{$1} = $gi_lookup_hash{$1};
			#print "tempHash\{" . $1 . "\} = " . $gi_lookup_hash{$1} . "\n";
		}
	}
}

#3: output tempHash to file
open out, ">$parentTaxaDataPath";

foreach my $gi (keys %tempHash) {
	print out $gi . " " . $tempHash{$gi} . "\n";
}

close out;
