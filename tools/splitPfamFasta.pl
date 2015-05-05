#!/usr/bin/perl -w



# pipe Pfam-A.fasta in this script, will create one FASTA file for each family in your /data/emptydirectory/


if (@ARGV != 1) {
	print "Usage: cat Pfam-A.fasta | splitPfamFasta.pl /data/emptydirectory/";
	print "   or: zcat Pfam-A.fasta.gz | splitPfamFasta.pl /data/emptydirectory/";
	exit(1);
}

my $output_dir=$ARGV[0];

unless (-d $output_dir) {
	print STDERR "error $output_dir is not a directory...\n";
	die;
}

my $currentPfam="";

my $FAMILY;


while (my $line=<STDIN>)
{

	if ($line =~ /^>\S*\s+(PF\d+\.\d+)/) {
		if ($1 ne $currentPfam) {
			$currentPfam = $1;
			if (defined $FAMILY) {
				close($FAMILY);
				}
			
			open $FAMILY, ">>$output_dir/$currentPfam.fas" or die $!;
		}
		
	}
	
	unless (defined $FAMILY) {
		die;
	}
	print $FAMILY $line;

}
