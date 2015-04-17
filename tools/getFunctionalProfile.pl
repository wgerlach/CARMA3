#!/usr/bin/env perl
# Written by Wolfgang Gerlach, University of Bielefeld (Germany)

########################################################################################
#### Configuration


# Download Ontology File ("OBO v1.2 format") from http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo
# Specify where you have saved the file (better use absolut path name):
use constant ONTOLOGY_FILE => '/vol/carma/data/gene_ontology.1_2.obo';

use constant GNUPLOT_TERMINAL => 'postscript enhanced color ",12"';


########################################################################################
use strict;
use Getopt::Std;
use File::Temp qw(tempfile tempdir);


sub init_ONTOLOGY_FILE(){

my $i=0;
	my $gene_ontology = {};
	my $alt_mapping = {};
	
	
	open (STREAM, ONTOLOGY_FILE);
	my ($go_id, $term, $category);
	my @alt_ids=();
	
	while (my $line = <STREAM>) {
		#print $line;

		if ($line =~ /^\[Term\]/) {
			if ((defined $go_id)&&(defined $term)&&(defined $category)) {
				#print "GOT: [$go_id, $term, $category]\n";
				$gene_ontology->{$go_id}=[$term, $category, 0];
				
				while (my $alt_id =shift(@alt_ids)) {
					$alt_mapping->{$alt_id}=$go_id;
				}
				
				@alt_ids = ();
				
				$i++;
			}
		
			$go_id=undef;
			$term=undef;
			$category=undef;
		} elsif ($line =~ /^id:\s*(.+)$/) {
			$go_id=$1;
		} elsif ($line =~ /^name:\s*(.+)$/) {
			$term=$1;
		} elsif ($line =~ /^namespace:\s*(.+)$/) {
			$category=$1;
		} elsif ($line =~ /^alt_id:\s*(.+)$/) {
			push(@alt_ids, $1);
		}

		
		
		#print "GOT: [$go_id, $term, $category]\n";
		
	}

	close(STREAM);
	
	#print "i: $i\n";
	
	return [$gene_ontology, $alt_mapping];
}


sub parseCARMAFile(){
	my ($input_file, $gene_ontology, $alt_mapping) = @_;

	my $sequence_count = 0;

	open (STREAM, $input_file);
	
	while (my $line = <STREAM>) {
	
		if ($line =~ /^\#/) { # comment lines
			next;
		}
	#print $line."\n";
		chomp($line);
		if ( (length($line) > 0) && ($line =~ /\>/) ) {
			$sequence_count++;
			my ($evalue, $go_list);
			unless ((undef, undef, $evalue, $go_list) = split(/\=\+\=/, $line) ) {
				print STDERR "error: could not parse line: $line\n";
				die;
			}
			
			#print "org: ".$evalue."\n";
			unless(defined $main::opt_e && $evalue > $main::opt_e) {
				# remove opening an closing brackets
				$go_list =~ s/\{//;
				$go_list =~ s/\}//;
				# TODO use evalue cutoff
			
				my @go_array = split(/\,/, $go_list);
			
				foreach my $go_id (@go_array) {
					#print "go_id: ($go_id)\n";
	
					unless ( defined $gene_ontology->{$go_id} ) {
				
						unless (defined $alt_mapping->{$go_id}) {
							print STDERR "undefined: $go_id\n";
							die;
							}
					
						$go_id=$alt_mapping->{$go_id};
						}
					# increase counter by one:
					$gene_ontology->{$go_id}[2]++;
				}
			}
		
		}
	
	}

	close(STREAM);

	if ($sequence_count == 0) {
		print STDERR "Error: no sequences found, maybe this is no EGT file ? (in FASTA-format, .egt !)\n";
		exit(1);
	}

	return;
}

sub usage{
	print STDERR <<USAGE;
	
Copyright (C) 2009 CeBiTec, Bielefeld University.
Written by Wolfgang Gerlach.
This software comes with ABSOLUTELY NO WARRANTY; This is free
software, and you are welcome to redistribute it under
certain conditions.
Read the COPYRIGHT for details.


Create a Gene Ontology based profile from CARMA output:

Usage: getFunctionalProfile.pl <options> carma_egt_file	
	

	Options
	
	-o <c>  : output file (tab separated values)
	-g <c>  : visualize profile with gnuplot/postscript
	            e.g.: "-g output.ps"
	-e <n>	: e-value cut-off for EGTs	            
	-s 	: sort data (recommended)
	-l <n>  : limit data to the n most abundant GO's
                    (recommended in conjunction with -g)
	-h      : this help

USAGE
}


##################
# Get options
##################
our($opt_h, $opt_e, $opt_g, $opt_s, $opt_l, $opt_o);

getopts('he:g:sl:o:');


if($opt_h){
	usage;
	exit;
}


my $output_file=$opt_o;

if (defined $output_file) {
	if (-e $output_file) {
		print STDERR $output_file . " already exists.\n";
		exit(1);
	}
}

if (defined $opt_e) {
	# simple check
	my $test = $opt_e+1.5;
	if ($test == 1.5 && not $opt_e =~ /0+\.?0*/) {
		print STDERR "ERROR: Value \"".$opt_e."\" for -e not a valid number.\n";
		exit(1);
		}		
}



unless (-e ONTOLOGY_FILE) {
	print STDERR "Could not find ONTOLOGY_FILE=\"".ONTOLOGY_FILE."\"\n";
	print STDERR "Configuration is done in this perl script.\n";
	die;
}

my $gnuplot_path;
if($opt_g) {

	$gnuplot_path = `which gnuplot`;
	
	if ($gnuplot_path) {
		print STDERR "gnuplot: $gnuplot_path\n";
	} else {
		print STDERR "ERROR: \"which\" could not find gnuplot.\n";
		exit(1);
	}
	
	if (-e $opt_g) {
		print STDERR "ERROR: File $opt_g already exists!\n";
		exit(1);
	}
	
}



unless (@ARGV) {
	usage();
	exit(0);
}


my $input_file = shift(@ARGV);




my ($gene_ontology, $alt_mapping) = @{init_ONTOLOGY_FILE()};
unless ( defined $gene_ontology ) {
	die;
}



&parseCARMAFile($input_file, $gene_ontology, $alt_mapping);


my ($fh, $tmp_file);
if ($opt_g) {
	($fh, $tmp_file) = tempfile(UNLINK => 1, SUFFIX => ".dat");
	#print "tmp_file: $tmp_file\n";
}

# write results into output array
my @output_array;
while ( my ($go_id, $array_ref) = each(%$gene_ontology) ) {
	my ($term, $category, $count) = @$array_ref;
	if ($count > 0) {
		#my $output_line = "$go_id\t$term\t$category\t$count\n";
		push(@output_array, [$go_id, $term, $category, $count])
        }
   }


if ($opt_s) {
	@output_array = sort {$b->[3]  <=>  $a->[3]} @output_array;
	
	if ($opt_l) {
		while (@output_array > $opt_l) {
			pop(@output_array);
		}
	}
}


if ($opt_g) {
	foreach my $ref (@output_array) {
		#my $output_line = join("\t", @$ref)."\n";
		#$output_line =~ s/ /-/;
		my ($go_id, $term, $category, $count) = @$ref;
		
		if ($category =~ /function/) {
			$category = "function";
		} elsif ($category =~ /process/) { 
			$category = "process";
			}
		elsif ($category =~ /component/) {
			$category = "component";
		}
		
       		print $fh "$go_id\t\"$term ($category)\"\t$category\t$count\n";
       	}
} 

if (defined $output_file) {
	open (STREAM, ">".$output_file);
	foreach my $ref (@output_array) {
		#print join("\t", @$ref)."\n";
		my ($go_id, $term, $category, $count) = @$ref;
		print STREAM "$go_id\t\"$term ($category)\"\t$category\t$count\n";
	}
	close(STREAM);
	print STDERR "File \"$output_file\" written.\n";
}


if ($opt_g)   {
	my $terminal = GNUPLOT_TERMINAL;

	my $gnuplot_commands= <<GNUPLOT;

#set size 1,0.5
set terminal $terminal
#set auto x
set output "$opt_g"
set boxwidth 0.9 absolute
set style fill   solid 1.00 border -1
set rmargin 5
set style histogram clustered gap 1 title  offset character 0, 0, 0
set style data histograms
set xtics border in scale 0,0.5 nomirror rotate by -45  offset character 0, 0, 0 font ",10"
set ylabel "Count of supporting EGTs"
set title "Functional Profile"
plot '$tmp_file' using 4:xticlabels(2) with histogram title ""
GNUPLOT


	open(PIPE, "| $gnuplot_path ") || die "gnuplot failed: $!\n";
	
	print PIPE $gnuplot_commands;
	close(PIPE);

	print STDERR "Gnuplot succeeded, so file \"$opt_g\" should have been written.\n";
		
}    

exit(0);


