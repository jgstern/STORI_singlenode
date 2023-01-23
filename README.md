STORI
=====

### Selectable Taxon Ortholog Retrieval Iteratively

#### For the motivation behind STORI, check out our PrePrint:
#### [Accessing and applying molecular history](https://dx.doi.org/10.7287/peerj.preprints.1293v1)

User Guide: [https://tinyurl.com/mr395m6z](https://tinyurl.com/mr395m6z)
Thesis: [http://linkd.in/1fZO63l](http://linkd.in/1fZO63l)

### 4/20/2015 - Introducing single-node STORI

https://github.com/jgstern/STORI_singlenode

Now STORI should run on any box with bash, Perl, and Python.
The initial release of STORI requires a cluster using the job
scheduler Moab, but this latest release runs on a single node.
I tested it on CentOS 7 using NCBI's modern accession format (rather than GIs).

The user guide for multi-node STORI should be sufficient for
setup and usage of single-node STORI, with a few adjustments.

1. STORIcontrol.pl is now STORIcontrol.py  
	I refactored STORIcontrol from Perl to Python, to help me learn
	Python. This means you'll need to have a working install of
	Python 2.7, and invoke STORIcontrol with  
		  >python STORIcontrol.py  
	STORIcontrol.py still needs to run in background to monitor convergence
	and continue runs as necessary, but the overhead is low.
2. The more runs you start, the slower everything goes [assuming only 1 core]  
	I changed beginSTORI.pl and continueSTORIfast_t.pl to kick off
	new STORI.pl PIDs rather than submit new jobs to a job scheduler.
3. Wall-clock limit for STORI.pl is uniform and now controlled
by maxRuntime in STORIcontrol.py  
	In the future I would like to make this limit adaptive and a
	unique property of each run.
4. Comparing CPU hours until convergence between different runs is now
complicated, since more running STORI.pl processes mean slower performance.
The runtime attribute should instead be thought of as "total process wall time for this run".
5. Fixed divide by zero bug in sub DecideReduction of STORI.pl.

Eventually I would like to integrate STORIstats with STORIcontrol,
and rewrite the methods for saving and reading the job_data_STORI
file to use JSON. 



### Method for preparing a VM for STORI:
1. Install CentOS 7

2. Installed perlbrew  
	\curl -L http://install.perlbrew.pl | bash  
	echo "source ~/perl5/perlbrew/etc/bashrc" >> .bash_profile

3. Installed GNUscreen, gcc, "Development Tools", et al.:  
	sudo su  
	yum install screen  
	yum install gcc  
	yum groupinstall "Development Tools"  
	yum install bzip2  
	yum install perl-core  
	yum install wget  
	exit

4. Installed Perl for ec2-user:  
	screen  
	perlbrew install perl-5.16.0

5. Installed cpanm for (I think) all perlbrew perls  
	perlbrew install-cpanm

6. Install openssl  
	sudo su  
	cd /usr/src  
	wget https://www.openssl.org/source/openssl-3.0.7.tar.gz  
	tar -zxf openssl-3.0.7.tar.gz  
	rm openssl-3.0.7.tar.gz  
	cd /usr/src/openssl-3.0.7  
	./config  
	make  
	make test  
	make install  
	ln -s /usr/local/lib64/libssl.so.3 /usr/lib64/libssl.so.3  
	ln -s /usr/local/lib64/libcrypto.so.3 /usr/lib64/libcrypto.so.3  
	yum install zlib-devel  
	exit  
	[exit terminal then reopen]  
	openssl version  

7. Install cpan modules  
	cpanm -i Statistics::Descriptive Data::Dumper List::MoreUtils Time::Elapse Bio::SeqIO Getopt::Long LWP::Simple LWP::UserAgent HTTP::CookieJar::LWP LWP::Protocol::https  

8. Install Python 2.7.9 (python 2.6 lacks a necessary module)  
	sudo wget https://www.python.org/ftp/python/2.7.9/Python-2.7.9.tgz  
	sudo chmod 755 Python-2.7.9.tgz  
	tar -xvf Python-2.7.9.tgz  
	cd Python-2.7.9  
	./configure  
	sudo make altinstall  
	cd /home/ec2-user  
	mkdir bin  
	cd bin  
	ln /home/ec2-user/Python-2.7.9/python python279
	
9. The original STORI repo contains the original BLAST+ executables. However, now it is 10 years later, and those binaries will crash if your FASTA deflines use a pdb accession format with more than 1 letter for the chain. Hence you should grab the latest version: 
	wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz  
	
	Next, do `tar -xvf ncbi-blast-2.13.0+-x64-linux.tar.gz` and grab these 3 binaries: makeblastdb, blastp, and blastdbcmd  


### End-to-end example
1. Set up a CentOS VM as above  
2. Make a directory in your home dir for STORI and put these scripts in it  
3. Decide which taxa are of interest and use their txids to populate ``taxids_GIs.txt``  
4. It is also best to populate ``taxids_GIs.txt`` with the Nucleotide accession(s).version for each chromosome of a completed genome for each taxon ID  
5. If the genome is not complete or the Nucleotide accession is otherwise unavailable, you can simply type ``na``  
6. In the latter scenario you might get organellar/plastid protein sequences; that is why a Nucleotide finished chromosome(s) accession.version is preferable  
7. ``perl getFastas.pl``  
8. ``perl makeNr.pl``  
9. ``python STORIcontrol.py``  

### For steps 3-4 above, how do I get the taxon IDs and Nucleotide accessions?
1. Go to https://www.ncbi.nlm.nih.gov/taxonomy  
2. In the search bar, enter 1[uid] and click Search  
3. You should see a single result that says "root" - the root of the tree of life - click it  
4. Now you should see "Archaea, Bacteria, Eukaryota, Viruses..." click any of them  
5. 

# TODO
- make the end-to-end example more detailed & complete
- better parameterize the various file paths at the beginning of each of the scripts, so that setup doesn't require updating the beginning of each script  
