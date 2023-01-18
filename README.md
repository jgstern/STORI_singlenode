STORI
=====

###Selectable Taxon Ortholog Retrieval Iteratively

####For the motivation behind STORI, check out our PrePrint:
####[Accessing and applying molecular history](https://dx.doi.org/10.7287/peerj.preprints.1293v1)

User Guide: [https://tinyurl.com/mr395m6z](https://tinyurl.com/mr395m6z)
Thesis: [http://linkd.in/1fZO63l](http://linkd.in/1fZO63l)

###4/20/2015 - Introducing single-node STORI

https://github.com/jgstern/STORI_singlenode

Now STORI should run on any box with bash, Perl, and Python.
The initial release of STORI requires a cluster using the job
scheduler Moab, but this latest release runs on a single node.
I tested it on RHEL 6.5.

The user guide for multi-node STORI should be sufficient for
setup and usage of single-node STORI, with a few minor adjustments.

1. STORIcontrol.pl is now STORIcontrol.py  
	I refactored STORIcontrol from Perl to Python, to help me learn
	Python. This means you'll need to have a working install of
	Python 2.7, and invoke STORIcontrol with  
		  >python STORIcontrol.py  
	STORIcontrol.py still needs to run in background to monitor convergence
	and continue runs as necessary, but the overhead is low.
2. The more runs you start, the slower everything goes  
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



###Method for preparing a VM for STORI:
1. Install CentOS 7

2. Installed perlbrew  
	\curl -L http://install.perlbrew.pl | bash  
	echo "source ~/perl5/perlbrew/etc/bashrc" >> .bash_profile

3. Installed GNUscreen:  
	sudo su  
	yum install screen  
	exit

4. Installed gcc  
	sudo su  
	yum install gcc  
	exit

5. Installed "Development Tools"  
	sudo su  
	yum groupinstall "Development Tools"  
	yum install bzip2  
	yum install perl-core  
	yum install wget  
	exit  

6. Installed Perl for ec2-user:  
	screen  
	perlbrew install perl-5.16.0

7. Installed cpanm for (I think) all perlbrew perls  
	perlbrew install-cpanm

8. Install openssl  
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

9. Install cpan modules  
	cpanm -i Statistics::Descriptive Data::Dumper List::MoreUtils Time::Elapse Bio::SeqIO Getopt::Long LWP::Simple LWP::UserAgent HTTP::CookieJar::LWP LWP::Protocol::https  

10. Install Python 2.7.9 (python 2.6 lacks a necessary module)  
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
	
11. To save GitHub space I omitted the BLAST executables from the STORI_singlenode repo, however you will find them in the STORI repo. so these commands should retrieve them:  
	wget https://github.com/jgstern/STORI/raw/master/blastdbcmd  
	wget https://github.com/jgstern/STORI/raw/master/blastp  
	wget https://github.com/jgstern/STORI/raw/master/makeblastdb  


# TODO
- make sure STORIstats.pl works with accessions rather than GIs  
- remove hard-coded file paths from continueSTORIfast_t.pl  
