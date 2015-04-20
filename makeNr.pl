#use lib '/nv/hp10/jstern7/perl5reinstall/lib';
#use lib '/nv/hp10/jstern7/perl5reinstall/lib/perl5';

use Statistics::Descriptive;

#takes fasta files as input, makes them nonredundant, and then makes them into blast databases

$projectDir = "/home/ec2-user/STORI/universal20150110"; 
$hitDir = $projectDir . "/hits";

#get taxa list
$taxaFile = "/home/ec2-user/STORI/taxids_GIs.txt";   #a list of taxids

$makeblastdbPath="/home/ec2-user/STORI/makeblastdb";  #path to the blast+ program makeblastdb
$blastDirPath = $projectDir . "/blast";

$nrdbLoc = "/home/ec2-user/STORI/bp_nrdb_SHA.pl";


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
}

sub MaxDex {
	my($ref_arr) = @_;
	my @array = @{$ref_arr};
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	$answer = $stat->maxdex();
	return $answer;
}

sub RemoveExtraGIs {
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
		if ($line =~ m/^(>.+?;)(\s:gi.+\n)/) { 		#match before semicolon nongreedily
			@deflines = split /gi\|/,$line; my $index=0;
			shift @deflines; #first elt is >
			
			foreach $def (@deflines) {
				my @gis=();
				if ($def =~ m/(\d+?)\|.*/) {
					push @gis, $1;	
				}
				$index = MaxDex(\@gis);
			}
			
			my $final;
			#$fastaArr[$r] = $1;
			if ($deflines[$index] =~ m/(.+)\;\s:$/) {
				$final = $1; }
			else {
				$final = $deflines[$index];
			}
			
			my $gi;
			if ($final =~ m/^(\d+)\D+/) {
				$gi = $1;
				print hitfile $gi . "\t" . $gi . "\t100\n";
			}
			
			$fastaArr[$r] = ">gi|" . $final;
			#print "match!! num1 is: $1 \n\n and num2 is: $2 \n\n";
		}
		else {
			chomp($line);
			my $gi;
			if ($line =~ m/^gi|(\d+)\D+/) {
				$gi = $1;
				print hitfile $gi . "\t" . $gi . "\t100\n";
			}
			$fastaArr[$r] = $line; 
		}
	
		print nrgifile "$fastaArr[$r]\n";
		
		$r++;
	}

	close nrfile;
	close nrgifile;
	close hitfile;
}


