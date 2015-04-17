#!/usr/bin/env perl
# Written by Wolfgang Gerlach, University of Bielefeld (Germany)

use FindBin;
use lib "$FindBin::Bin/../lib/";

use PerlTools;

########################################################################################
#### Configuration


#Download and unpack
#ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
use constant NAMES_FILE => '/vol/biodb/ncbi_taxonomy/names.dmp'; 
use constant NODES_FILE => '/vol/biodb/ncbi_taxonomy/nodes.dmp';
use constant MERGED_FILE => '/vol/biodb/ncbi_taxonomy/merged.dmp';

use constant DEFAULT_CONSIDERED_RANKS => ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species');

# other ranks you might want to consider:
#'tribe', 'forma', 'kingdom', 'subclass', 'varietas', 'subkingdom', 'suborder', 'infraorder'



use constant UNKNOWN_TAXA_RANK => 'superkingdom';						
use constant UNKNOWN_LABEL => 'Unknown';


#according ncbi-taxonomy these taxa(directly below root) are assigned no rank:					
#28384	other sequences (incl. Plasmids and Transposons)
#12908	unclassified sequences 
#12884	Viroids
#10239	Viruses
use constant NCBI_CELLULAR_ORGANISMS => 131567;
use constant NON_CELLULAR_ORGANISMS_NAME => 'Non-cellular';
use constant NON_CELLULAR_ORGANISMS_RANK => 'superkingdom';

use constant GNUPLOT_TERMINAL => 'postscript enhanced color "Arial"';



########################################################################################
use strict;
use Getopt::Std;
use File::Temp qw(tempfile tempdir);
use File::Basename;

my $NCBI_NAMES;
my $NCBI_NODES;
my $CONSIDERED_RANKS_HASH;
my @CONSIDERED_RANKS;

#sub init_CONSIDERED_RANKS_HASH() {
#	my ($optr) = @_;
#	
#	my %hash=();
#	
#	
#	if (defined $optr) {
#		@CONSIDERED_RANKS = split(/,/,$optr);
#	} else {
#		@CONSIDERED_RANKS = DEFAULT_CONSIDERED_RANKS;
#	}
#	
#	foreach my $rank (@CONSIDERED_RANKS) {
#			%hash->{lc($rank)}=1;
#		}
#	
#	unless (%hash) {
#		die;
#		}

#	$CONSIDERED_RANKS_HASH=\%hash;

#}





sub parseCARMAFile(){
	my ($input_file) = @_;

	my $counts_per_rank;
	
	my $ncbi_code;
	my $rank;
	my $name;
	my $parent;

	unless (-e $input_file) {
		print STDERR "error: file \"".$input_file."\" not found\n";
		exit(1);
	}

	open (STREAM, $input_file);
	
	while (my $line = <STREAM>) {
		chomp($line);
		my $predicted_ncbi_code=undef;
		(undef, undef, undef, $predicted_ncbi_code) = split(/\t/, $line);
		unless ( defined $predicted_ncbi_code) {
			print STDERR "Error parsing file.\n";
			print STDERR "LINE: $line\n";
			die;
		}
		$ncbi_code=$predicted_ncbi_code;
		
		if ($ncbi_code==0) {
			my $label = UNKNOWN_LABEL;
			$counts_per_rank->{lc( UNKNOWN_TAXA_RANK )}->{$label}++;
			next ;
		}
		
		my $is_cellular=0;
		# check if it is a cellular organism
		while ($ncbi_code!=1) {
		
			if ($ncbi_code == NCBI_CELLULAR_ORGANISMS) {
				$is_cellular=1;
				last;
			}
			unless (defined @{$NCBI_NODES->{$ncbi_code}}) {
				print STDERR "Line: ".$line."\n";
				print STDERR "Got nothing for ncbi_code  \"$ncbi_code\" in ncbi_nodes, skipping...\n";
				#exit(1);
				last;
			}
				
			unless ( ($parent, $rank) = @{$NCBI_NODES->{$ncbi_code}} ) {
				die;
			}
			
			$ncbi_code= $parent;
	
		
		}
		
		if ($is_cellular==0) {
			my $nc_name = NON_CELLULAR_ORGANISMS_NAME;
			my $nc_rank = NON_CELLULAR_ORGANISMS_RANK;
			$counts_per_rank->{lc($nc_rank)}->{$nc_name}++;
			next;
		}
		
		
		$ncbi_code=$predicted_ncbi_code;
		
		
		while ($ncbi_code!=1) {
		
			unless (defined @{$NCBI_NODES->{$ncbi_code}}) {
				print STDERR "Got nothing for ncbi_code $ncbi_code in ncbi_nodes, skipping...\n";
				last;
			}
				
			unless ( ($parent, $rank) = @{$NCBI_NODES->{$ncbi_code}} ) {
				die;
			}
		

			if ( (defined $CONSIDERED_RANKS_HASH->{lc($rank)} && $CONSIDERED_RANKS_HASH->{lc($rank)}==1 )) {
				$name = $NCBI_NAMES->{$ncbi_code};
				
				unless ($name) {
					print STDERR "error: nothing found for ncbi_code $ncbi_code in NCBI_NAMES.\n";
					die;
				}
				
				$counts_per_rank->{lc($rank)}->{$name}++;
				
			}
			
			#next ncbi_code:
			$ncbi_code= $parent;
	
		
		}
		
	
	}

	close(STREAM);

	return $counts_per_rank;
}


sub call_gnuplot {
	my ($useheader, $numberOfDataSets, $data_file, $rank, $gnuplot_path, $out_dir, $fontsize) = @_;

	#my $numberOfDataSets = @$files;
	my $terminal = GNUPLOT_TERMINAL." ".$fontsize;
	
	my $output_file = $rank.".ps";

	
	
	
	if (defined $out_dir) {
		$output_file = File::Spec->catfile( ($out_dir), $output_file );
	}

	my $include_header = "set key autotitle columnhead";
	if ($useheader == 0) {
		#deactivate:
		$include_header= "\#".$include_header;
	}

	my $use_grey = 0;
	my $color = "";
	if ($use_grey == 1) {
		$color = "lc rgb \"black\"";
	}
	
	my $plotoptions="using 3:xticlabels(2) $color with histogram ";
	if ($useheader == 0) {
		$plotoptions.="title \"\"";
	}
	for (my $i = 4; $i < $numberOfDataSets+3; $i++) {

		#$plotoptions .= " , \'\' u $i with histogram ";
		
		if ($use_grey == 1) {
			if ($i == 4) {
				$color = "lc rgb \"dark-grey\"";
			} elsif ($i == 5) {
				$color = "lc rgb \"light-grey\"";
			}
		}
		
		$plotoptions .= " , \'\' u $i $color with histogram ";
		if ($useheader == 0) {
			$plotoptions.="title \"\"";
		}
	}
	$rank = ucfirst($rank); #Takes a string and retruns it with the first character in upper case. 
	
	my $gnuplot_commands= <<GNUPLOT; # (font ",10")


#set size 1,0.5
set terminal $terminal
#set auto x
#set palette gray
#unset colorbox
set yrange [0:]
set output "$output_file"
set boxwidth 0.9 absolute
set style fill solid 1.00 border -1
set rmargin 5
set style histogram clustered gap 1 title  offset character 0, 0, 0
set style data histograms
set xtics border in scale 0,0.5 nomirror rotate by -45  offset character 0, 0, 0 
set ylabel "Relative Abundance"
set title "Taxonomic Profile - $rank"
$include_header
#plot '$data_file' using 3:xticlabels(2) with histogram title "" , '' u 4   with histogram title ""
plot '$data_file' $plotoptions
GNUPLOT


	open(PIPE, "| $gnuplot_path ") || die "gnuplot failed: $!\n";
	
	print PIPE $gnuplot_commands;
	close(PIPE);

	print STDERR "Gnuplot succeeded, so file \"$output_file\" should have been written.\n";

}

sub usage{

	my @default_ranks_array=();
	foreach my $rank (keys %$CONSIDERED_RANKS_HASH) {
		if ($CONSIDERED_RANKS_HASH->{$rank} == 1) {
			push(@default_ranks_array, $rank);
		}
    	}
	my $default_ranks = join(",",@default_ranks_array);
	
	print STDERR <<USAGE;
	
Copyright (C) 2009 CeBiTec, Bielefeld University.
Written by Wolfgang Gerlach.
This software comes with ABSOLUTELY NO WARRANTY; This is free
software, and you are welcome to redistribute it under
certain conditions.
Read the COPYRIGHT for details.


Create a taxonomic profile from CARMA output:

Usage: getComparativeTaxonomicProfile.pl <options> carma_classifications_file
	

	Options
	
	-r <c>  : specify which ranks you want to consider
	            default:
	            $default_ranks
        -c <n>  : cut-off for species abundance
                    (recommended, e.g. "-c 0.01")
	-i <n>  : sort data accoding to i'th data set (default: 1)
	-g      : visualize profile with gnuplot/postscript
	            e.g. produces superkingdom.ps, family.ps ...
	-a      : use absolut values instead of relative values	            
	            
	Options only together with -g:
	
	-o <c>  : directory where you want to store the files
	            (default: current directory)
	-n <n>  : specify gnuplot legend. 
	            Argument: comma-separated list
	            E.g. "-n readsA,readsB,readsC"
	            Use "-n no" to hide legend.
	            (default legend: filenames)
	-f <n>  : font size (default: 10)
	-p <c>  : "-p 1" - convert ps to pdf
                  "-p 2" - convert ps to pdf & pdfcrop


USAGE
}

########################################################################################

##################
# Get options
##################
our($opt_h, $opt_r, $opt_g, $opt_a, $opt_o, $opt_i, $opt_c, $opt_n, $opt_f, $opt_p,);

getopts('hr:gao:i:c:n:f:p:');

#init_CONSIDERED_RANKS_HASH();
$CONSIDERED_RANKS_HASH= &PerlTools::init_CONSIDERED_RANKS_HASH(\@CONSIDERED_RANKS, DEFAULT_CONSIDERED_RANKS);

if($opt_h){
	usage;
	exit;
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
	
	
	
}



unless (@ARGV) {
	usage();
	exit(1);
}

unless (defined $opt_f) {
	$opt_f = 10;
}

if (defined $opt_o) {
	unless (-e $opt_o) {
		print STDERR "ERROR Directory \"$opt_o\" does not exist.\n";
		exit(1);
	}
	
	unless (-d $opt_o) {
		print STDERR "ERROR \"$opt_o\" is not a directory.\n";
		exit(1);
	}

	unless ( File::Spec->file_name_is_absolute($opt_o) ) {
		print STDERR "ERROR \"$opt_o\" is not an absolute path.\n";
		exit(1);
	
	}
	
}


my $ps2pdf;
my $pdfcrop;
if (defined $opt_p) {
	#system("ps2pdf $output_file");
	$ps2pdf = `which ps2pdf`;	
	chomp($ps2pdf);
	unless ($ps2pdf) {
		print STDERR "ERROR: Could not find ps2pdf.\n";
	}
	if ($opt_p == 2) {
		$pdfcrop = `which pdfcrop`;
		chomp($pdfcrop);
		unless ($pdfcrop) {
			print STDERR "ERROR: Could not find pdfcrop.\n";
		}

	}
	
}

if (defined $opt_r) {
	#&init_CONSIDERED_RANKS_HASH($opt_r);
	#my @temp_array=split(/,/,$opt_r);
	$CONSIDERED_RANKS_HASH= &PerlTools::init_CONSIDERED_RANKS_HASH(\@CONSIDERED_RANKS, split(/,/,$opt_r));
}


#my $input_file = shift(@ARGV);

##################################################################################



#init_NCBI_NAMES();
$NCBI_NAMES = &PerlTools::init_NCBI_NAMES(NAMES_FILE, MERGED_FILE);
$NCBI_NODES = &PerlTools::init_NCBI_NODES(NODES_FILE, MERGED_FILE);
#init_NCBI_NODES();

my ($fh, $data_file);
my @absoulte_profiles;
my @files;


# read absolute numbers from carma output
my $input_file;
while ( $input_file = shift(@ARGV) ) {
	push(@files, $input_file);

	my $counts_per_rank = &parseCARMAFile($input_file);
	push(@absoulte_profiles, $counts_per_rank);
}


my $numberOfDataSets=@files;
my @dataset_names;


if (defined $opt_n) {
	unless ($opt_n eq "no") {
		@dataset_names = split(/,/, $opt_n);
		my $num=@dataset_names;
		if ($num != $numberOfDataSets) {
			print STDERR "ERROR: $numberOfDataSets data sets (input files), but $num data set names: \"$opt_n\".\n";
			exit(1);
		}
	} 
} else {
		@dataset_names=@files;
}

# search for every possible taxon
my %bighash;
foreach my $file (@absoulte_profiles) {
	foreach my $rank (keys %$file) {
		my $hash = $file->{$rank};
		foreach my $taxon (keys %$hash) {
			$bighash{$rank}->{$taxon}=0;
		}
	}
}


my $sort_i=1;
if (defined $opt_i) {
	$sort_i = $opt_i;
}
my $master_file = $absoulte_profiles[$sort_i-1];



# put counts of master_file into bighash
foreach my $rank (keys %$master_file) {
	my $taxa_hash = $master_file->{$rank};

	foreach my $taxon (keys %$taxa_hash) {
		$bighash{$rank}->{$taxon}=$master_file->{$rank}->{$taxon};
	}
}


my %bighash_sorted;
#convert "hash of hashes" into "hash of sorted arrays":
foreach my $rank (keys %bighash) {
	my $rankhash = $bighash{$rank};
	foreach my $taxon (sort {$rankhash->{$b} <=> $rankhash->{$a} } keys %$rankhash) {
		push ( @{$bighash_sorted{$rank}} , $taxon);
	}

}

print "\n";

my $cutoff = 0;
if (defined $opt_c) {
	$cutoff = $opt_c;
}


my %total_counts;
#compute relative abundance BEFORE cut-off
#foreach my $file (@absoulte_profiles) {
for (my $i = 0 ; $i < @absoulte_profiles; $i++) {
	my $file = $absoulte_profiles[$i];
	my $filename = $files[$i];
	my $dataset_name = $dataset_names[$i];
	foreach my $rank (keys %bighash) {
		my $hash = $file->{$rank};
		my $totalcount=0;
		foreach my $taxon (keys %$hash) {
			$totalcount += $hash->{$taxon};
		}
		$total_counts{$rank}->{$dataset_name}=[$totalcount, -1];
		
		foreach my $taxon (keys %$hash) {
			if ($totalcount > 0) {
				my $absolute_abundance=$hash->{$taxon};
			
				my $relative_abundance = $absolute_abundance / $totalcount;
				
				$hash->{$taxon} = [$absolute_abundance, $relative_abundance];
			} else {
				$hash->{$taxon} = [0,0];
			}
		}
	}
}

my %filtered_ranks;
#remove taxa which are below cut-off (in every data set)
foreach my $rank (keys %bighash_sorted) {
	my $rankarray = $bighash_sorted{$rank};

	foreach my $taxon (@$rankarray) {
		
		
		#check if taxon is in any dataset above cut-off:
		my $isAboveCutOff=0;
		foreach my $file (@absoulte_profiles) {
			my $r_abundance = 0;
			
			if (defined $file->{$rank}->{$taxon}) {
				$r_abundance = @{$file->{$rank}->{$taxon}}[1];
			}

			if ($r_abundance > $cutoff) {
				$isAboveCutOff = 1;
			}
		}
		
		#if taxon is always below cut-off, delete it:
		if ($isAboveCutOff == 0) {
		
			foreach my $file (@absoulte_profiles) {
				if (defined $file->{$rank}->{$taxon}) {
					$filtered_ranks{$file}->{$rank} += @{$file->{$rank}->{$taxon}}[0];
					delete $file->{$rank}->{$taxon};
				}
			}
		
		}
		
		
	}

}



#recompute relative abundance AFTER cut-off
#foreach my $file (@absoulte_profiles) {
for (my $i = 0 ; $i < @absoulte_profiles; $i++) {
	my $file = $absoulte_profiles[$i];
	my $filename = $files[$i];
	my $dataset_name = $dataset_names[$i];
		
	foreach my $rank (keys %bighash) {

		my $hash = $file->{$rank};
		my $totalcount=0;
		foreach my $taxon (keys %$hash) {
			$totalcount += @{$hash->{$taxon}}[0];
		}
		@{$total_counts{$rank}->{$dataset_name}}[1]= $totalcount;
		
		foreach my $taxon (keys %$hash) {
			if ($totalcount > 0) {
				my ($absolute_abundance)= @{$hash->{$taxon}};
			
				my $relative_abundance = $absolute_abundance / $totalcount;
				@{$hash->{$taxon}}[1] = $relative_abundance;
			} 
		}
	}
}

# Output data (total counts):

print "Other:\tabsolut\tabsolut after cut\tOther relative\n";

print "rank";
foreach my $dataset_name (@dataset_names) {
	print "\t$dataset_name\t--\t--";
}
print "\n";

foreach my $rank (@CONSIDERED_RANKS) {
	my $rankhash = $total_counts{$rank};
	
	unless (defined $rankhash) {
		print STDERR "skip $rank...\n";
		next;
	}
	print "$rank";
	foreach my $dataset_name (@dataset_names) {
		unless (defined $rankhash->{$dataset_name}) {
			print STDERR "error: ". $dataset_name ." ($rank) not found in hash\n";
			die;
		}
		my ($tot, $totc) = @{$rankhash->{$dataset_name}};
		if ($tot==0) {
			print "\t$tot\t$totc\t0";
		} else {
			print "\t$tot\t$totc\t".sprintf("%.4f", (1-$totc/$tot));
		}
	}
	print "\n";	
}


# Output data:
#foreach my $rank (keys %bighash_sorted) {
foreach my $rank (@CONSIDERED_RANKS) {
	unless (defined $bighash_sorted{$rank} ) {
		next;
	}

	my $i=0;
	if ($opt_g) {
		($fh, $data_file) = tempfile(UNLINK => 0, SUFFIX => ".dat");
		 
	}
	
	
	
	my $header = "\nrank\t\"taxon\"";
	my $useheader=0;
	if (defined $opt_n) {
		my $descr = $opt_n;
		$descr =~ s/,/\t/g;
		$header.= "\t$descr\n";
	} else {
		foreach my $file (@files) {
			$header.="\t".basename($file);
		}
		$header.="\n";
	}

	unless ( (defined $opt_n) && ($opt_n eq "no") ) {
		$useheader=1;
		
		if ($opt_g) {
			print $fh $header;
		} else {
			print $header;
		}
	}
	
	my $rankarray = $bighash_sorted{$rank};
	foreach my $taxon (@$rankarray) {
		
		my @data_line;

		my $taxon_exists=0;
		foreach my $file (@absoulte_profiles) {
			my $abundance = 0;
			
			if (defined $file->{$rank}->{$taxon}) {
				$taxon_exists = 1;
				if (defined $opt_a) {
					$abundance = @{$file->{$rank}->{$taxon}}[0];
				} else {
					$abundance = sprintf("%.4f",@{$file->{$rank}->{$taxon}}[1]);
				}
			}
			push (@data_line, $abundance);
		}
		
		if ($taxon_exists == 1) {
			my $output_line = "$rank\t\"$taxon\"\t" . join("\t",@data_line) . "\n";
		
			if ($opt_g) {
				$i++;
				print $fh $output_line;
			} else {
				print $output_line;
			}
		}
		
	}
	if ($opt_g) {
		close($fh);
		if ($i > 0) {
			&call_gnuplot($useheader, $numberOfDataSets, $data_file, $rank, $gnuplot_path, $opt_o, $opt_f);
			}
		
	}

}


if (defined $opt_p && defined $opt_g) {
	foreach my $rank (keys %bighash_sorted) {
		my $ps_file = $rank.".ps";
		if (defined $opt_o) {
			$ps_file = File::Spec->catfile( ($opt_o), $ps_file );
		}

		my $pdf_file = $ps_file;
		$pdf_file=~ s/\.ps$/\.pdf/;

		print "call: $ps2pdf $ps_file $pdf_file\n";
		system("$ps2pdf $ps_file $pdf_file");
		
		if ($opt_p == 2) {
			
			print "call: $pdfcrop $pdf_file $pdf_file\n";
			system("$pdfcrop $pdf_file $pdf_file");
		}

	}
	
}

print "\nreads that have been filtered away by cut-off:\n";
foreach my $rank (@CONSIDERED_RANKS) {
	print "$rank";
	foreach my $file (@absoulte_profiles) {
		
		if (defined $filtered_ranks{$file}->{$rank}) {
			print "\t".$filtered_ranks{$file}->{$rank};
	
		} else {
			print "\t0";
		}
		

	}
	print "\n";
}


print "\nrelative filtered away:\n";
foreach my $rank (@CONSIDERED_RANKS) {
	print "$rank";
	my $i = 0;
	foreach my $file (@absoulte_profiles) {

		my $dataset_name = $dataset_names[$i];
		my $total = @{$total_counts{$rank}->{$dataset_name}}[1];
		if ($total > 0) {
		print "\t".sprintf("%.4f",$filtered_ranks{$file}->{$rank}/$total);
		} else {
			print "\t--";
		}

		
		
		$i++;
	}
	print "\n";
}

exit(0);



