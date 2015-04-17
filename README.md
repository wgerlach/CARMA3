
  Copyright (C) 2010 CeBiTec, Bielefeld University
  Written by Wolfgang Gerlach.

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  version 2 as published by the Free Software Foundation.

  This file is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this file; see the file LICENSE.  If not, write to
  the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
  Boston, MA 02111-1307, USA.


CARMA3
======

This software is freely available at
https://github.com/wgerlach/CARMA3



If you use CARMA please cite:
Wolfgang Gerlach, Jens Stoye
Taxonomic classification of unassembled metagenomic reads with CARMA3
Nucleic acids research 39 (14), e91-e91


Content
=======

* Installation
* Instructions for use




Installation
============

a)
Requires parts of the Boost C++ Libraries

under ubuntu:
```bash
sudo apt-get install libboost-all-dev libbz2-dev zlib1g-dev
```

manually:
www.boost.org/
Download and install the boost library if not already available on your system:
```bash
wget http://downloads.sourceforge.net/project/boost/boost/1.51.0/boost_1_51_0.tar.gz
tar xvfz boost_1_51_0.tar.gz
cd boost_1_51_0
```

Just in case you want to use boost version 1.46: sorry, you'll have to manually modify one line in libs/filesystem/v3/src/path.cpp
See details https://svn.boost.org/trac/boost/attachment/ticket/4688/boost_filesystem.2.patch

```bash
./bootstrap.sh --with-libraries=iostreams,program_options,filesystem,system --prefix=/yourpath/CARMA3/
./bjam toolset=gcc --disable-filesystem2 address-model=64 install
```

b)
Compile CARMA:

under ubuntu you might need to install automake first:
```bash
sudo apt-get install automake
```

```bash
./autogen.sh
```
If autogen.sh fails, maybe you need new macro files (m4_ax_boost*.m4 in CARMA directory) for the BOOST library, they can be found here: http://www.gnu.org/software/autoconf-archive/The-Macros.html

If the boost library is installed in system default directories (e.g. via apt-get), then you may omit "--with-boost" in the following:

```bash
./configure  --with-boost=/yourpath/CARMA3/ LDFLAGS="-m64 -Wl,-R /yourpath/CARMA3/lib"
```
if your boost install prefix was /yourpath/CARMA3/

```bash
make
```

c)
Install BLAST (with the NR database) and/or HMMER (with Pfam database, see next step for details)
CARMA3 has been developed to parse BLAST output produced with blastall-2.2.21 and blastall-2.2.24. We did not test other versions yet. 
Unfortunately, CARMA3 does not yet work with BLAST+. Please use the legacy BLAST, version up to 2.2.26.
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/

Let us know if you encounter any problems.

d)
If you want to use HMMER, CARMA3 needs:
```bash
wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam24.0/database_files/pfamA.txt.gz
wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam24.0/database_files/gene_ontology.txt.gz
wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam24.0/Pfam-A.hmm.gz
wget ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62
```

These files will be needed only once, afterwards they can be deleted:
```bash
wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam24.0/database_files/pfamseq.txt.gz
wget ftp://ftp.sanger.ac.uk/pub/databases/Pfam/releases/Pfam24.0/Pfam-A.fasta.gz
```

Create pfamId2TaxId-mapping file:
```bash
cd src
g++ -o getPfamNCBITaxiID getPfamNCBITaxiID.cpp
cd ..
./src/getPfamNCBITaxiID pfamseq.txt.gz | gzip > pfamid2taxid.txt.gz
```

Create pfam_fasta_dir, a directory that contains the multiple alignments of each Pfam family:
```bash
mkdir pfam_fasta_dir
zcat Pfam-A.fasta.gz | tools/splitPfamFasta.pl ./pfam_fasta_dir/ 
```

e) (not recommended, experimental)
If you want to use RDP for 16S classification:
download and modify RDP files (replace "10_24" by the current version):
```bash
wget http://rdp.cme.msu.edu/download/release10_24_arch_aligned.fa.gz
wget http://rdp.cme.msu.edu/download/release10_24_bact_aligned.fa.gz
wget http://rdp.cme.msu.edu/download/release10_24_unaligned.fa.gz
./src/carma --create_specific_RDP_ali --input release10_24_arch_aligned.fa.gz | gzip > release10_24_arch_aligned_specific.fa.gz
./src/carma --create_specific_RDP_ali --input release10_24_bact_aligned.fa.gz | gzip > release10_24_bact_aligned_specific.fa.gz
./src/carma --create_specific_RDP_seq --input release10_24_unaligned.fa.gz | gzip > release10_24_unaligned_specific.fa.gz
```

Create RDP BLAST database
```bash
zcat release10_24_unaligned.fa.gz | formatdb -i stdin -n rdp -p F -o T
```

f) 
download NCBI taxonomy files:
```bash
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
```
extract the files and place them somewhere...

g)
Configure carma.cfg 

h)
Several perl scripts are located in the tools subdirectory. For an explanation see WebCARMA manual section online.
These scripts also require configuration, like the paths to the NCBI taxonomy and GO files. You can modiy corresponding entries
in the headers of these scripts.

i)
for SGE submission, adapt the following line in carma3-SGE.pl to your SGE configuration:
```perl
my $native_options = "-cwd -l arch=sol-amd64,";
```

In some situations the DRMAA causes "Text file busy"-errors, that are related to the submission of shell scripts.
In this case add the following parameters to $native_options:
```text
"-shell yes -b n -S /bin/shell"
```
(thanks to Koen Illeghems for this hint)



Instructions for use
====================

CARMA3 can provide taxonomic classifcation with BLASTx(NR) and HMMER3(Pfam). Our evaluations have shown
that taxonomic classification accuracy of the BLASTx-based variant is better than the HMMER-variant. Therefore 
we currently recommend to use the BLASTx-variant of CARMA3 for taxonomic classification.

For the functional classification we recommend to use our HMMER-variant of CARMA3. If you use only the HMMER-variant
of CARMA3, note that detection of frameshifts is not supported. But if you use the BLAST-variant prior to the 
HMMER-variant, the HMMER-variant can also use the read translation of BLASTx which includes frameshifts.


CARMA3 main binary:
```bash
./src/carma
```
You can also work with this binary on one single machine, if you have not too much data to analyse. Also, when 
you already have BLAST results from somewhere else, you can use the binary for the taxonomic classification only as this 
is usually much faster.

```bash
./carma3-SGE.pl
```
If you have a compute cluster with Sun Grind Engine as submission system, you can use our pipeline script, to perform
the blast search, classification of blast results, HMMER3-search, HMMER classification and creation of different output
files (e.g. visualisation of the taxonomic and functional profiles as histograms) in one step. If you have some other
submission system you may be able to adapt our script for your system. In this case we would happy to include your script
in the source code of CARMA3.



If you find bugs, have any problems or questions, please feel free to contact me via
mailAwolfgang-gerlach.com
(Important, please replace the "A" with "@"!)


Wolfgang Gerlach


