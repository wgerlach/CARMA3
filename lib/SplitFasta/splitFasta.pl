#!/usr/bin/perl -w

use FindBin;
use lib "$FindBin::Bin";
use SplitFasta;


1;


if (@ARGV != 4) {
	print "Usage: splitFasta NumberOfChunks FASTA-file output-directory prefix\n";
	print "e.g.:  splitFasta 100 input.fas /tmp/fastadir/ chunk\n";
	print "Numbering starts with 1 !\n";
	exit(1);
}

my $numberOfChunks=$ARGV[0];
my $fasta_file=$ARGV[1];
my $output_dir=$ARGV[2];
my $prefix = $ARGV[3];

SplitFasta::splitFasta($numberOfChunks, $fasta_file, $output_dir, $prefix);


