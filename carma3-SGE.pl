#!/usr/bin/env perl

#
#  Copyright (C) 2007 CeBiTec, Bielefeld University
#  Written by Wolfgang Gerlach.
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  version 2 as published by the Free Software Foundation.
#
#  This file is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this file; see the file LICENSE.  If not, write to
#  the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
#  Boston, MA 02111-1307, USA.
#

use strict;

require Scheduler::DRMAA;
use FindBin;
use lib "$FindBin::Bin/lib/SplitFasta";
use lib "$FindBin::Bin/lib/Schedule";
use lib "$FindBin::Bin/lib";

use SplitFasta;
use FileUtils;
use File::Temp qw(tempfile tempdir);
use File::Glob ':glob';
use File::Copy;
use Getopt::Std;
use Cwd;

use File::Path;

1;

our($opt_p, $opt_P,$opt_M, $opt_h, $opt_c, $opt_l, $opt_m, $opt_x, $opt_n, $opt_d, $opt_a, $opt_e, $opt_I, $opt_C);

my $config_file = "carma.cfg";
my $prefix = "chunk_";


#defined by configuration file:
my $carma_binary; 
#my $blastall_script;
#my $blastx_database;
#my $blastx_evalue; 
#my $blastn_database;
#my $blastn_evalue; 
my $cluster_tmp_dir;
my $maxNumberOfChunks; 
my $message; 
my $LD_LIBRARY_PATH;
my $clusteremail;
my $notificationemail;
############################


sub usage{
	print STDERR <<USAGE;
	
Copyright (C) 2010 CeBiTec, Bielefeld University.
Written by Wolfgang Gerlach.
This software comes with ABSOLUTELY NO WARRANTY; This is free
software, and you are welcome to redistribute it under certain conditions.
Read the COPYRIGHT for details.


carma3-SGE.pl: taxonomic classification of short
                  environmental DNA fragments
	
Usage: carma3-SGE.pl <options> <fasta_input_file>
  where: fasta_file contains metagenomic DNA-fragments

	Options
	-p <c>  : project output directory (required)
	-m      : use Pfam HMM (for DNA)
	-M      : use Pfam HMM (for amino acids)	
	-x      : use BLASTX
	-P      : use BLASTP
	-c <n>	: submit <n> jobs
	-l <n>  : at most <n> CPUs are used on the cluster
	-d      : delete intermediate results (e.g. BLAST result)
	-a      : configuration file
	-C      : pass to --config_overlay option (experimental)
	-e <c>  : email notification via UNIX "mail"
	-h	: show this help message


USAGE
}


sub DIE_handler {
	my($signal) = @_;

	my $mail_message = "CARMA3 pipeline failed";

	if (defined $signal) {
		$mail_message .= ": ".$signal;
	}
	mailme($mail_message);
}
$SIG{__DIE__}  = 'DIE_handler';



sub mailme {
	my ($mail_message) = @_;

	if (defined $opt_e) {

	
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);

		my $timestamp = sprintf "%4d-%02d-%02d %02d:%02d:%02d", $year+1900,$mon+1,$mday,$hour,$min,$sec;

		$mail_message = "[".$timestamp."] ".$mail_message;
	
		my $mail_command = "echo \"Subject:CARMA3 job status notification\\n".$mail_message."\" \| mail ".$opt_e;
	
		print STDERR "call: ".$mail_command."\n";
	
		system( $mail_command );
	}
}


my %config_file_hash;

sub readConfiguration {

	my ($_config_file) = @_;
	
	print $carma_binary."\n";
	
	

	open(FILE,$_config_file)  or die $!;
	while (<FILE>) { 
		
		#print $_; 
		chomp($_);
		if ($_ =~ /^#/) {
			next;
		}
		if ($_ =~ /^\s*$/) {
			next;
		}
		if (my ($key, $value) = $_ =~ /(\S+)\s*\=\s*(\S+)/) {
			if (my ($inner_value) = $value =~ /^\"(\S*)\"$/ ) {
				$value = $inner_value;
			}
			#print $key.": ".$value."\n";
			$config_file_hash{$key} = $value;
		} else {

			die("Error parsing Configuration file: \"".$_."\"\n");
		}
	}

	close (FILE);

	# take values from configuration hash

	$carma_binary = getConfigValue("carma_binary");
	
	$clusteremail = getConfigValue("clusteremail");
	#$blastall_script = getConfigValue("blastall_script");
	#$blastx_database = getConfigValue("blast_nr_database");
	#$blastn_database = getConfigValue("blast_nt_database");
	
	#$blastx_evalue = getConfigValue("blastx_evalue");
	#$blastn_evalue = getConfigValue("blastn_evalue");
	
	$maxNumberOfChunks = getConfigValue("maxNumberOfChunks");
	$message = getConfigValue("message");
	$cluster_tmp_dir = getConfigValue("cluster_tmp_dir");
	$LD_LIBRARY_PATH = getConfigValue("LD_LIBRARY_PATH");


}

sub getConfigValue() {
	my ($key) = @_;
	
	unless (defined $config_file_hash{$key}) {
		die( "Error: Key \"".$key."\" not found in configuration file !\n");
		
	}
	return $config_file_hash{$key};
}

sub submit_SGE_jobs_meta{
	my ($actual_command_line, $this_project_dir, $num_jobs, @l_options) = @_;
	
	my $max_repeats = 5;
	my $success_jobs = 0;
	
	while ($success_jobs < $num_jobs && $max_repeats > 0) {
	
		$success_jobs = submit_SGE_jobs($actual_command_line, $this_project_dir, $num_jobs, @l_options);
		if ($success_jobs == $num_jobs) {
			return 0;
		}
		print STDERR "Warning: Only $success_jobs of $num_jobs succeeded. Submit again... ($max_repeats)\n";
		$max_repeats--;
		sleep 10;
	}

	die ("Error: Cluster submission failed. Only $success_jobs of $num_jobs succeeded.\n");

	return 1;
}

sub submit_SGE_jobs{
	my ($actual_command_line, $this_project_dir, $num_jobs, %l_options) = @_;
	

	my $stop_file = $this_project_dir."/../stop";
	my $last_stop_test = time();
	if (-e $stop_file) {
		print STDERR "warning: stop file detected\n";
		exit(0);
	}
	
	require Scheduler::DRMAA;
		 
	my @job_cmd_array;
	
	my $template = "CARMA3_XXXXXX";
	my ($fh, $script_file) = tempfile ($template, UNLINK => 1, DIR => $cluster_tmp_dir, SUFFIX => ".sh");
	
	print STDERR "script_file: $script_file\n";
	

	push(@job_cmd_array, $script_file );
	

 	
 	open (MYFILE, '>'.$script_file) or die;
	print MYFILE "#!/bin/bash\n";
	print MYFILE "ulimit -c 0\n"; # to avoid huge dump cores
	#print MYFILE "export LD_LIBRARY_PATH=".$LD_LIBRARY_PATH."\n";
	print MYFILE $actual_command_line."\n";
	close (MYFILE); 
	
 	chmod(0775, $script_file); 
 	
 	# open file 
 	#open (MYFILE, '<'.$script_file) or die;
 	#unlink $script_file;
 		 

	sleep 10;
	
	unless (-d $this_project_dir) {
		print STDERR $message."project directory not found.\n";
		exit(1);
	}
	
	#my $tmp_dir = FileUtils::get_tmp_dir_object($cluster_tmp_dir);
	
	#my $errorfile = $this_project_dir."\$SGE\_TASK\_ID.stderr";
	
	#delete old stderr-files:
	my @old_stderr_files = glob $this_project_dir."*.stderr";
	foreach my $file (@old_stderr_files) {
		unlink $file or die "Could not unlink $file: $!";
	}
	
	

	my $scheduler = Scheduler::DRMAA->new(temp_dir => $this_project_dir, commandline => [@job_cmd_array]);
	
	#my $native_options = "-cwd -l arch=sol-amd64,vf=20G";
	my $native_options = "-cwd -l arch=sol-amd64";
	
	while ( my ($key, $value) = each(%l_options)) {
		$native_options.=",$key=$value";
	}
	unless (defined $opt_I) {
		$native_options.=",idle=1";
	}
	
	if ( $clusteremail =~ /\@/ ) {
		$native_options.=" -m a -M ".$clusteremail;
		#$native_options.=" -m a -M wgerlach\@cebitec.uni-bielefeld.de";
	}
	
	if (defined $opt_l) { #max_running_tasks -- limit number of CPUs, not jobs
		$native_options.=" -tc ".$opt_l;
	}
	
	#if ($LD_LIBRARY_PATH ne "") { #does not work!
	#	$native_options .= " -v LD_LIBRARY_PATH\=".$LD_LIBRARY_PATH;
	#}
	print STDERR "use native_options: $native_options\n";
	$scheduler->set_native_option($native_options);
	
	
	sleep 60; # 60 #work-around such that filesystem has time to sync..
	my $jobids = $scheduler->submit2($num_jobs, 0);
	
	my $arrayid = @$jobids[0];
	($arrayid) = $arrayid =~ /^(\d+)\./;
	print STDERR "Array job ID: ".$arrayid."\n";	

	print STDERR "Done.\nSubmitted $num_jobs jobs!\n";
	system("echo $arrayid > ".$this_project_dir."/jobid.txt");
	

	
	
	my @error_jobs=();
	my $success_jobs = 0;
	
	print STDERR "waiting and parsing results...\n";
	$scheduler->iterate2(
		sub {
			my ($name) = @_;

			
			$success_jobs++;
			print STDERR "Job \"" . $name . "\" finished successfully. ";
			if ($num_jobs != (@error_jobs+$success_jobs)) {
				my $left = ($num_jobs-(@error_jobs+$success_jobs));
				if ($left > 1) {
					print STDERR  "$left jobs are left.\n";
				} else {
					print STDERR  "1 job is left.\n";
				}
			} else {
				print STDERR "\n";
			}
			if (($last_stop_test + 600) >= time()) {
				$last_stop_test = time();
				if (-e $stop_file) {
					print STDERR "warning: stop file detected\n";
					exit(0);
				}
			}
		},
		sub {
		 	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);

			printf STDERR "%4d-%02d-%02d %02d:%02d:%02d ", $year+1900,$mon+1,$mday,$hour,$min,$sec;

			printf STDERR "Job $arrayid\.%s failed: %s\n",$_[0],$_[1];

			push (@error_jobs, $_[0]);
			
			if (($last_stop_test + 600) >= time()) {
				$last_stop_test = time();
				if (-e $stop_file) {
					print STDERR "warning: stop file detected\n";
					exit(0);
				}
			}
		});
		
	if (@error_jobs > 0) {
		print STDERR "Sorry, list of jobs that failed: " . join(',', @error_jobs ) . "\n";
		print STDERR "Command was: ".$actual_command_line."\n";
		print STDERR $message."Cluster job crashed.\n";
		#print STDERR "Simply rerun CARMA with same parameters, only jobs that failed will be resubmitted!\n";
		#print STDERR "(Will not work if called via carma.pl using option \"-o\"!)\n";
		#return 1;
	} else {
		print STDERR "All submitted jobs finished successfully!\n";
	}
	
	unlink $this_project_dir."/jobid.txt"  or warn "Could not unlink ".$this_project_dir."/jobid.txt: $!";
	
	unlink $script_file;
	$scheduler->DESTROY();
	
	# delete stdout-files, since they should be empty...
	my @old_stdout_files = glob $this_project_dir."*.stdout";
	foreach my $file (@old_stdout_files) {
		unlink $file or die "Could not unlink $file: $!";
	}
	
	return $success_jobs;
}

sub makeDirectory{
	my ($dir) = @_;
	
	unless(-d $dir) {
		print "create $dir\n";
		mkdir $dir or die "Couldn't create dir: [$dir] ($!)";
	}

}

sub getArrayMinPosition{

	my ($array_ref) = @_;
	
	my $size = @{$array_ref};
	
	my $min_value=$$array_ref[0];
	my $min_position=0;
	
	for (my $i = 0; $i < $size; $i++) {
		if ($$array_ref[$i] == 0) {
			return $i;
		}
	
		if ($$array_ref[$i] < $min_value) {
			$min_value = $$array_ref[$i];
			$min_position=$i;
		}
	}
	
	
	return $min_position;
}

sub mergeTaxFiles{
	my ($output_dir, $result_file) = @_;

	unless (-e $result_file) {
		my $files_wildcard = "$output_dir/*.tax";

		my @files = glob($files_wildcard);

		unless (@files) {
			print STDERR "ERROR: No (tax) result files found !??\n";
			return 1;
		}

		unless (open(OUTPUT, ">".$result_file)) {
			die $!;
		}

		foreach my $file (@files) {
			#print $file."\n";
			copy($file,\*OUTPUT);

		} 

		close OUTPUT;
	} else {
		print STDERR "File \"$result_file\" already exists, skipping... \n";
		
	}

	return 0
}


sub combineTaxResultsByEvalue {
	my $outputfile = shift(@_);
	my @files = @_;
	my @hashes;
	my $superhash;
	my $i=0;
	
	# read files
	foreach my $file (@files) {
		print $file."\n";
		unless (open(FILE, $file)) {
			die $!;
		}
		
		while (<FILE>) {
			#print $_;
			chomp($_);
			my ($desc, $prop, $goterms, $taxid, $tax_string, $evalue) = split(/\t/, $_);
			
			my ($read_id) = $desc =~ /^(.*)\_-?\d\_-?\d$/;
			unless ($read_id) {
				$read_id = $desc;
			}
			#print $desc."\n";
			$superhash->{$read_id} = 1;
			if (defined $hashes[$i]->{$read_id}) {
				die ( "already defined???\n");
				
			}
			$hashes[$i]->{$read_id} = [$desc, $prop, $goterms, $taxid, $tax_string, $evalue];
		}
		close(FILE);
		$i++;
	}
	
	
	open(OUTFILE, ">".$outputfile)  or die $!;
	
	# iterate trough all reads
	for my $key ( keys %$superhash ) {
        	#my $value = $superhash->{$key};
	        #print "$key => $value\n";
	        my $min_classi;
	        my $min_evalue=100000;
	       # my $count = 0;
	        foreach my $hash (@hashes) {
	        	if (defined $hash->{$key}) {
		        	#$count++;
	        		#print $hash->{$key}[4]."\n";
	        		
	        		if (($hash->{$key}[5] < $min_evalue) && ($hash->{$key}[5] > 0)) {
	        			$min_evalue = $hash->{$key}[5];
	        			$min_classi = $hash->{$key};
	        		}
	        		
	        	}
	        	
	        }
	        
	        #print "max: $min_evalue\n";
	        if (defined $min_classi) {
	        	print OUTFILE join("\t", @$min_classi)."\n";
	        }
	        #if ($count > 1) {
	        #	print "c: $count\n";
	        #	exit(0);
	        #}
	        
    	}
	close(OUTFILE);
}

sub combineTaxResults {
	my $outputfile = shift(@_);
	my @files = @_;
	my @hashes;
	my $superhash;
	my $i=0;
	
	open(OUTFILE, ">".$outputfile)  or die $!;
	
	# read files
	foreach my $file (@files) {
		print "read file: ".$file."\n";
		unless (open(FILE, $file)) {
			die $!;
		}
		
		while (<FILE>) {
			#print $_;
			chomp($_);
			my ($desc, $prop, $goterms, $taxid, $tax_string, $evalue) = split(/\t/, $_);
			
			my ($read_id) = $desc =~ /^(.*)\_-?\d\_-?\d$/;
			unless ($read_id) {
				$read_id = $desc;
			}
			#print $desc."\n";
			
			unless (defined $superhash->{$read_id}) {
				print OUTFILE $_."\n";
				$superhash->{$read_id} = 1;
			}
			
			#if (defined $hashes[$i]->{$read_id}) {
			#	print STDERR "already defined???\n";
			#	abort();
			#}
			#$hashes[$i]->{$read_id} = [$desc, $prop, $goterms, $taxid, $tax_string, $evalue];
		}
		close(FILE);
		$i++;
	}
	
	
	close(OUTFILE);
}

sub get_egt_by_family_hash {
	my ($egt_result_file) = @_;	

	my $_egt_by_family_hash = {};	
	
	open(FILE, "<".$egt_result_file)  or die $!;
		my $pfam_family;
		my $descr;
		my $sequence ="";
		while (<FILE>) { 
			#print $_; 
			chomp($_);
	
			if ($_ =~ /^#/) {
				next;
			}
	
			if ($_ =~ /^>/) {
				 ($pfam_family) = $_ =~/^>(PF\d+\.\d+)\=\+\=/;
				 unless (defined $pfam_family) {

				 	die("error parsing line:\n"."\"".$_."\"\n");
				 }
				 $descr = $_;
			} else {
				$sequence = $_;
				push ( @{$_egt_by_family_hash->{$pfam_family}} , [$descr, $sequence]);
			}
		}

	close (FILE);	
	
	return $_egt_by_family_hash;
}



sub collectEGTs() {
	my ($_egt_result_file, $search_hmm_output_dir)= @_;
	
	unless (-e $_egt_result_file) {
		my $files_wildcard = "$search_hmm_output_dir/*.egt";

		my @files = glob($files_wildcard);

		unless (@files) {
			print STDERR "ERROR: No EGT files found !??\n";
			print STDERR $message."No homology could be detected.\n";
			exit(0);
		}

		unless (open(OUTPUT, ">".$_egt_result_file)) {
			die $!;
		}

		foreach my $file (@files) {
			copy($file,\*OUTPUT);
		} 

		close OUTPUT;
	} else {
		print STDERR "File \"$_egt_result_file\" already exists, skipping... \n";
	}
	
	my $egt_result_file_filesize = -s $_egt_result_file; //;;;
	return $egt_result_file_filesize;
	
}




sub createNewEGTchunks {
	my ($egt_by_family_hash, $classify_hmm_input_dir) = @_;

	my $_realNumberOfEGTChunks;


	my $numberOfFamilies = keys( %$egt_by_family_hash );
	
	if ($numberOfFamilies < $maxNumberOfChunks) {
		$_realNumberOfEGTChunks = $numberOfFamilies;
	} else {
		$_realNumberOfEGTChunks = $maxNumberOfChunks;
	}


	my $first_chunkfilename = $classify_hmm_input_dir.$prefix."1.egt";

	if (-e $first_chunkfilename) {

		print STDERR "EGT chunks already exist, skip creation...\n";
	} else {

		my @num_seq_per_file=();
		my @egt_by_chunk=();

		for (my $i=0; $i < $maxNumberOfChunks ;$i++) {
			$num_seq_per_file[$i]=0;
		}

		foreach my $family (keys %{$egt_by_family_hash}){

			my $array_ref = $egt_by_family_hash->{$family};
			my $pos = getArrayMinPosition(\@num_seq_per_file);
			my $num_egt_in_family = @{$array_ref};
	
	
			$num_seq_per_file[$pos] += $num_egt_in_family;
			push( @{$egt_by_chunk[$pos]}, @{$array_ref} );
	
		}

		

		for (my $i=0; $i < $_realNumberOfEGTChunks ;$i++) {

			my $chunkfilename = $classify_hmm_input_dir.$prefix.($i+1).".egt";
			#print "chunkfilename: $chunkfilename\n";
			my $CHUNK_FILE;
			open $CHUNK_FILE, ">>".$chunkfilename or die $!;
	
			foreach my $tupel (@{$egt_by_chunk[$i]} ) {
				my ($descr, $seq) = @{$tupel};
				print $CHUNK_FILE $descr."\n";
				print $CHUNK_FILE $seq."\n";
			}
	
			close($CHUNK_FILE);
	
		}

	}

	return $_realNumberOfEGTChunks;
}



sub all_files_exist{
	my ($path, $prefix, $suffix, $number) = @_;
	
	
	my $all_exist = 1;
	for (my $i = 1; $i <= $number; $i++) {
		unless (-e $path.$prefix.$i.$suffix) {
			$all_exist = 0;
			last;
		}
	}

	return $all_exist;
}




######################################################################################################################################
# START

require Scheduler::DRMAA;
# example:
#./carma3-SGE.pl -p /vol/cluster-data/wgerlach/test/ /vol/metagenomics/metagenomes/25genomes/metasim/25metagenome_454sim_default_formatted.fna



##################
# Get options
##################


getopts('p:PMc:l:hmxda:e:IC:');


####################
# check parameters
####################
unless(@ARGV==1){
	usage;
	print STDERR "FASTA input file is missing.\n";
	exit 1;
}

unless(defined($opt_p)){
	usage;
	print STDERR "Argument \"-p\" is missing.\n";
	exit 1;
}

my $project_dir = $opt_p;
my $first_char = substr($project_dir,0,1);
if ($first_char ne "/") {
	print STDERR "error: please provide absolut path for project directory\n";
	exit(1);
}

my $fasta_file=$ARGV[0];
unless(-e $fasta_file) {
	print STDERR "ERROR: Could not read file \"$fasta_file\"!\n";
	exit 1;
}


unless(defined($opt_m) || defined($opt_x) || defined($opt_n) || defined($opt_P)|| defined($opt_M)){
	usage;
	print STDERR "ERROR: give at leat one of the options -m, -x, -M or -P\n";
	exit 1;	
}

if ( (defined($opt_x) || defined($opt_m)) && defined($opt_P)) {
	print STDERR "ERROR: combination of -P with others not supported\n";
	exit 1;	
}

if ( (defined($opt_x) || defined($opt_m)) && defined($opt_M)) {
	print STDERR "ERROR: combination of -M with others not supported\n";
	exit 1;	
}



unless (-d $project_dir) {
	die("Directory $project_dir not found.\n");
}

if (defined $opt_a) {
	if (-e $opt_a) {
		$config_file = $opt_a;
	}
}

my $config_overlay_arg="";
if (defined $opt_C) {
	$config_overlay_arg = " --config_overlay ".$opt_C." ";
}

################################	
# config file

#make absolue path name:
unless (substr($config_file, 0, 1) eq "/") {

	#my $thisdir = `pwd`;
	my $thisdir = &Cwd::cwd();
	$config_file = $thisdir."/". $config_file;
	
}


unless (-e $config_file) {
	die("Configuration file \"$config_file\" not found..\n".
		"Better provide full path to config file.\n")

}


readConfiguration($config_file);



if (defined $opt_c) {
	$maxNumberOfChunks = $opt_c;
}



################################
# create project sub-directories:

my $last_char = substr($project_dir,-1,1);
if ($last_char ne "/") {
	$project_dir .= "/";
}

my $search_input_dir = $project_dir."search_input/";

my $search_hmm_output_dir = $project_dir."search_hmm_output/";
my $search_blastx_output_dir = $project_dir."search_blastx_output/";
my $search_blastn_output_dir = $project_dir."search_blastn_output/";
my $search_blastp_output_dir = $project_dir."search_blastp_output/";


my $classify_hmm_input_dir = $project_dir."classify_hmm_input/";
my $classify_blastx_input_dir = $search_blastx_output_dir;
my $classify_blastn_input_dir = $search_blastn_output_dir;
my $classify_blastp_input_dir = $search_blastp_output_dir;


my $classify_hmm_output_dir = $project_dir."classify_hmm_output/";
my $classify_blastx_output_dir = $project_dir."classify_blastx_output/";
my $classify_blastn_output_dir = $project_dir."classify_blastn_output/";
my $classify_blastp_output_dir = $project_dir."classify_blastp_output/";


unless(-d $search_input_dir) {
	mkdir $search_input_dir or die "Couldn't create dir: [$search_input_dir] ($!)";
}


if ((defined $opt_m) || (defined $opt_M)) {
	makeDirectory($search_hmm_output_dir);
	makeDirectory($classify_hmm_input_dir);
	makeDirectory($classify_hmm_output_dir);
}


if (defined $opt_x) {
	makeDirectory($search_blastx_output_dir);
	makeDirectory($classify_blastx_output_dir);
}

if (defined $opt_P) {
	makeDirectory($search_blastp_output_dir);
	makeDirectory($classify_blastp_output_dir);
}


if (defined $opt_n) {
	makeDirectory($search_blastn_output_dir);
	makeDirectory($classify_blastn_output_dir);
}

################################
#read fasta-file, distribute sequences into chunks:


my $realNumberOfFASTAChunks = SplitFasta::splitFasta($maxNumberOfChunks, $fasta_file, $search_input_dir, $prefix);

if ($realNumberOfFASTAChunks == -1) {
	print STDERR $message."No FASTA-sequence found!\n";
	die;
}

################################
#submit search jobs:



if (defined $opt_x) {

	# ---------------------------------------------------------------------------------------------------------------
	# search-blastx:
	
	my $all_exist = all_files_exist($search_blastx_output_dir, $prefix, ".blastx", $realNumberOfFASTAChunks);

	if ($all_exist == 1) {
		print STDERR "BLASTX results already there...\n";
	} else {
		print STDERR "Submit jobs for BLASTX search...\n";
		my $blastx_search_command = $carma_binary.
						$config_overlay_arg.
						" --error ".$search_blastx_output_dir.$prefix."\$SGE_TASK_ID.stderr".
						" --config ".$config_file.
						" --local --delay \$SGE_TASK_ID".
						" --blast --type n --database p --input ".$search_input_dir.$prefix."\$SGE_TASK_ID.fas".
						" --output ".$search_blastx_output_dir.$prefix."\$SGE_TASK_ID.blastx";

		print STDERR "blastx_search_command: $blastx_search_command\n";
		my $ret = submit_SGE_jobs_meta($blastx_search_command, $search_blastx_output_dir, $realNumberOfFASTAChunks);
		
		if ($ret == 1) {
			print STDERR "Since blastx crashed, I delete all results...\n";
			my $rm_command = "rm ".$search_blastx_output_dir.$prefix."\*\.blastx";
			system($rm_command);
			die;
		}
	}
	
	# ---------------------------------------------------------------------------------------------------------------
	# classify-blastx (blastx):
	
	my $all_exist = all_files_exist($classify_blastx_output_dir, $prefix, ".tax", $realNumberOfFASTAChunks);

	if ($all_exist == 1) {
		print STDERR "blastx tax results already there...\n";
	} else {
		my $blastx_classify_command = $carma_binary.
						$config_overlay_arg.
						" --error ".$classify_blastx_output_dir.$prefix."\$SGE_TASK_ID.stderr".
						" --local --delay \$SGE_TASK_ID".
						" --config ".$config_file.
						" --classify-blast --type n --database p --input ".$classify_blastx_input_dir.$prefix."\$SGE_TASK_ID.blastx".
						" --blast-egts ".$classify_blastx_output_dir.$prefix."\$SGE_TASK_ID.blastxegt".
						" --output ".$classify_blastx_output_dir.$prefix."\$SGE_TASK_ID.tax";

		print STDERR "blastx_classify_command: $blastx_classify_command\n";
		my $ret = submit_SGE_jobs_meta($blastx_classify_command, $classify_blastx_output_dir,$realNumberOfFASTAChunks, 'vf' => '20G');
		if ($ret == 1) {
			die;
		}
	}

	
}

if (defined $opt_P) {

	# ---------------------------------------------------------------------------------------------------------------
	# search-blastp:
	
	my $all_exist = all_files_exist($search_blastp_output_dir, $prefix, ".blastp.gz", $realNumberOfFASTAChunks);

	if ($all_exist == 1) {
		print STDERR "BLASTP results already there...\n";
	} else {
		print STDERR "Submit jobs for BLASTP search...\n";
		my $blastp_search_command = $carma_binary.
						$config_overlay_arg.
						" --error ".$search_blastp_output_dir.$prefix."\$SGE_TASK_ID.stderr".
						" --config ".$config_file.
						" --local --delay \$SGE_TASK_ID".
						" --blast --type p --database p --input ".$search_input_dir.$prefix."\$SGE_TASK_ID.fas".
						" --gzip".
						" --output ".$search_blastp_output_dir.$prefix."\$SGE_TASK_ID.blastp.gz";
						

		print STDERR "blastp_search_command: $blastp_search_command\n";
		my $ret = submit_SGE_jobs_meta($blastp_search_command, $search_blastp_output_dir, $realNumberOfFASTAChunks);
		
		if ($ret == 1) {
			print STDERR "Since blastp crashed, I delete all results...\n";
			my $rm_command = "rm ".$search_blastp_output_dir.$prefix."\*\.blastp";
			system($rm_command);
			die;
		}
	}
	
	# ---------------------------------------------------------------------------------------------------------------
	# classify-blastp (blastp):
	
	my $all_exist = all_files_exist($classify_blastp_output_dir, $prefix, ".tax", $realNumberOfFASTAChunks);

	if ($all_exist == 1) {
		print STDERR "blastp tax results already there...\n";
	} else {
		my $blastp_classify_command = $carma_binary.
						$config_overlay_arg.
						" --error ".$classify_blastp_output_dir.$prefix."\$SGE_TASK_ID.stderr".
						" --local --delay \$SGE_TASK_ID".
						" --config ".$config_file.
						" --classify-blast --type p --database p --input ".$classify_blastp_input_dir.$prefix."\$SGE_TASK_ID.blastp.gz".
						" --fasta-input ".$search_input_dir.$prefix."\$SGE_TASK_ID.fas".
						" --output ".$classify_blastp_output_dir.$prefix."\$SGE_TASK_ID.tax";
						

		print STDERR "blastp_classify_command: $blastp_classify_command\n";
		my $ret = submit_SGE_jobs_meta($blastp_classify_command, $classify_blastp_output_dir,$realNumberOfFASTAChunks);
		if ($ret == 1) {
			die;
		}
	}

	
}



if (defined $opt_n) {

	# ---------------------------------------------------------------------------------------------------------------
	# search-blastn:
	
	my $all_exist = all_files_exist($search_blastn_output_dir, $prefix, ".blastn", $realNumberOfFASTAChunks);


	if ($all_exist == 1) {
		print STDERR "BLASTN results already there...\n";
	} else {
		print STDERR "Submit jobs for BLASTN search...\n";
		my $blastn_search_command = $carma_binary.
						$config_overlay_arg.
						" --error ".$search_blastn_output_dir.$prefix."\$SGE_TASK_ID.stderr".
						" --config ".$config_file.
						" --local --delay \$SGE_TASK_ID".
						" --blast --type n --database n --input ".$search_input_dir.$prefix."\$SGE_TASK_ID.fas".
						" --output ".$search_blastn_output_dir.$prefix."\$SGE_TASK_ID.blastn";
		
		print STDERR "blastn_search_command: $blastn_search_command\n";
		my $ret = submit_SGE_jobs_meta($blastn_search_command, $search_blastn_output_dir, $realNumberOfFASTAChunks);
		if ($ret == 1) {
			print STDERR "Since blastn crashed, I delete all results...\n";
			my $rm_command = "rm ".$search_blastn_output_dir.$prefix."\*\.blastn";
			system($rm_command);
			die;
		}

	}
	
	
	
	
	# ---------------------------------------------------------------------------------------------------------------
	# classify-blastn (blastn):
	
	my $all_exist = all_files_exist($classify_blastn_output_dir, $prefix, ".tax", $realNumberOfFASTAChunks);

	if ($all_exist == 1) {
		print STDERR "blastn tax results already there...\n";
	} else {

		my $blastn_classify_command = $carma_binary.
						$config_overlay_arg.
						" --error ".$classify_blastn_output_dir.$prefix."\$SGE_TASK_ID.stderr".
						" --local --delay \$SGE_TASK_ID".
						" --config ".$config_file.
						" --classify-blast --type n --database n --input ".$classify_blastn_input_dir.$prefix."\$SGE_TASK_ID.blastn".
						" --output ".$classify_blastn_output_dir.$prefix."\$SGE_TASK_ID.tax";

		print STDERR "blastn_classify_command: $blastn_classify_command\n";
		my $ret = submit_SGE_jobs_meta($blastn_classify_command, $classify_blastn_output_dir, $realNumberOfFASTAChunks);
		if ($ret == 1) {
			die;
		}
	}

}




if ((defined $opt_m) || (defined $opt_M)) {

	my $all_exist = all_files_exist($search_hmm_output_dir, $prefix, ".egt", $realNumberOfFASTAChunks);
	
	my $blast_egts_param = "";
	
	if ((defined $opt_x) && (defined $opt_m)) {
		$blast_egts_param = " --blast-egts ".$classify_blastx_output_dir.$prefix."\$SGE_TASK_ID.blastxegt";
	}

	my $type = "p";
	if (defined $opt_m) {
		$type = "n";
	}

	if ($all_exist == 1) {
		print STDERR "EGTs already there...\n";
	} else {
		print STDERR "Submit jobs for EGT search...\n";
		my $egt_search_command = $carma_binary.
				$config_overlay_arg.
				" --error ".$search_hmm_output_dir.$prefix."\$SGE_TASK_ID.stderr".
				" --local --delay \$SGE_TASK_ID".
				" --config ".$config_file.
				" --hmmer --type ".$type." --input ".$search_input_dir.$prefix."\$SGE_TASK_ID.fas".
				$blast_egts_param .
				" --output ".$search_hmm_output_dir.$prefix."\$SGE_TASK_ID.egt";

		print STDERR "egt_search_command: $egt_search_command\n";
		my $ret = submit_SGE_jobs_meta($egt_search_command, $search_hmm_output_dir , $realNumberOfFASTAChunks);
	
		if ($ret == 1) {
			die;
		}
	}





	################################
	#collect EGTs:  (don't need to collect blast results, they can directly be used for classification)

	my $egt_result_file = $project_dir."result.egt";

	if ((defined $opt_m) || (defined $opt_M)) {
		my $egt_result_file_filesize = &collectEGTs($egt_result_file, $search_hmm_output_dir);
	
		if ($egt_result_file_filesize <= 5) {
			print STDERR $message."No EGTs were found.\n";
		
			unless (defined $opt_x || defined $opt_n) {
				print STDERR $message."No homology could be detected.\n";
				exit(0);
			}
		}
	}
	
	################################
	#EGTs: sort by family:

	my $egt_by_family_hash;

	if ((defined $opt_m) || (defined $opt_M)) {	
		$egt_by_family_hash = get_egt_by_family_hash($egt_result_file);	
	}

	################################
	#EGTs: make chunks:

	my $realNumberOfEGTChunks;
	if ((defined $opt_m) || (defined $opt_M)) { 
		$realNumberOfEGTChunks = createNewEGTchunks($egt_by_family_hash, $classify_hmm_input_dir);
	}	


	################################
	#submit classification jobs:





	my $all_exist = all_files_exist($classify_hmm_output_dir, $prefix, ".tax", $realNumberOfEGTChunks);

	if ($all_exist == 1) {
		print STDERR "Hmm tax results already there...\n";
	} else {
		my $egt_classify_command = $carma_binary.
						$config_overlay_arg.
						" --error ".$classify_hmm_output_dir.$prefix."\$SGE_TASK_ID.stderr".
						" --local --delay \$SGE_TASK_ID".
						" --config ".$config_file.
						" --classify-egt --type p --input ".$classify_hmm_input_dir.$prefix."\$SGE_TASK_ID.egt".
						" --output ".$classify_hmm_output_dir.$prefix."\$SGE_TASK_ID.tax";

		print STDERR "egt_classify_command: $egt_classify_command\n";
		my $ret = submit_SGE_jobs_meta($egt_classify_command, $classify_hmm_output_dir, $realNumberOfEGTChunks);
		if ($ret == 1) {
			die;
		}
	}
	
}







################################
#collect results:


my $tax_hmm_result_file = $project_dir."hmm_result.tax";
my $tax_blastx_result_file = $project_dir."blastx_result.tax";
my $tax_blastp_result_file = $project_dir."blastp_result.tax";
my $tax_blastn_result_file = $project_dir."blastn_result.tax";
my $tax_blast_result_file;

if ((defined $opt_m) || (defined $opt_M)) {
	my $ret = mergeTaxFiles($classify_hmm_output_dir, $tax_hmm_result_file);
	if ($ret != 0) {
		print STDERR $message."No classification results for the HMM-based approach found.\n";
	}
}

if (defined $opt_x) {
	$tax_blast_result_file = $tax_blastx_result_file;
	my $ret = mergeTaxFiles($classify_blastx_output_dir, $tax_blastx_result_file);
	if ($ret != 0) {
		print STDERR $message."No classification results for the BLASTX-based approach found.\n";
	}
} elsif (defined $opt_P) {
	$tax_blast_result_file = $tax_blastp_result_file;
	my $ret = mergeTaxFiles($classify_blastp_output_dir, $tax_blastp_result_file);
	if ($ret != 0) {
		print STDERR $message."No classification results for the BLASTX-based approach found.\n";
	}
}

if (defined $opt_n) {
	my $ret = mergeTaxFiles($classify_blastn_output_dir, $tax_blastn_result_file);
	if ($ret != 0) {
		print STDERR $message."No classification results for the BLASTN-based approach found.\n";
	}
}


################################
#result.tax:


#if ((defined $opt_x) && (defined $opt_n) && (defined $opt_m)) {
#
#	my $tax_result_file = $project_dir."result.tax";
#	combineTaxResults($tax_result_file, $tax_blastx_result_file, $tax_blastn_result_file, $tax_hmm_result_file);
#}



######################################################################################
# Everything is done, now produce some nice pictures:

#### functional profile:
my $result_ps=$project_dir."/functional_profile.ps";
my $result_tsv=$project_dir."/functional_profile.tsv";
my $cmd;

my $files_wildcard;
my @files;
my $result_file;


if ((defined $opt_m) || (defined $opt_M)) {
	unless (-e $result_tsv) {
		$cmd=$FindBin::Bin."/tools/getFunctionalProfile.pl";
		$cmd.= " -o $result_tsv -g $result_ps -s -l 40 $project_dir/result.egt"; #for gnuplot output

		unless (system($cmd) == 0) {
			print STDERR "ERROR: getFunctionalProfile.pl stopped for some unknown reason.\n";
			print STDERR "command was: \"$cmd\"\n";
			exit(1);
		}
	} else {
		print STDERR "File \"$result_ps\" already exists, skipping... \n";
	}
}

#### taxonomic profile:
$result_tsv=$project_dir."/taxonomic_profile.tsv";


unless (-e $result_tsv) {
	if (-e $tax_blast_result_file) {
		$cmd=$FindBin::Bin."/tools/getTaxonomicProfile.pl";
		$cmd.= " -l 40 -o $result_tsv -g $project_dir $tax_blast_result_file";

		unless (system($cmd) == 0) {
			print STDERR "ERROR: getTaxonomicProfile.pl stopped for some unknown reason.\n";
			print STDERR "command was: \"$cmd\"\n";
			exit(1);
		}
	}
} else {
	print STDERR "(At least) \"$result_tsv\" already exists, skipping... \n";
}

##### convert ps2pdf


unless (-e "$project_dir/superkingdom.pdf") {

	my $ps2pdf = `which ps2pdf`;
	chomp($ps2pdf);

	if ($ps2pdf) {

		$files_wildcard = "$project_dir/*.ps";

		@files = glob($files_wildcard);
	
		foreach my $file (@files) {
			my $outputfile = $file;
			$outputfile =~ s/ps$/pdf/;
			$cmd = "$ps2pdf $file $outputfile";
			if (system($cmd) != 0) {
				print STDERR "ERROR calling \"$cmd\", returned: $?\n";
				exit(1);
			}

		}
		
		#delete postscripts, keep only pdf:
		system ("rm $project_dir/*.ps");
		
	} else {
		print STDERR "ERROR: ps2pdf not found\n";
	}
	
	

} else {
	print STDERR "(At least) \"$project_dir/superkingdom.pdf\" already exists, skipping... \n";
}





my $search_input_dir = $project_dir."search_input/";

my $search_hmm_output_dir = $project_dir."search_hmm_output/";
my $search_blastx_output_dir = $project_dir."search_blastx_output/";
my $search_blastp_output_dir = $project_dir."search_blastp_output/";
my $search_blastn_output_dir = $project_dir."search_blastn_output/";
my $classify_hmm_input_dir = $project_dir."classify_hmm_input/";
my $classify_blastx_input_dir = $search_blastx_output_dir;
my $classify_blastn_input_dir = $search_blastn_output_dir;
my $classify_blastp_input_dir = $search_blastp_output_dir;
my $classify_hmm_output_dir = $project_dir."classify_hmm_output/";
my $classify_blastx_output_dir = $project_dir."classify_blastx_output/";
my $classify_blastp_output_dir = $project_dir."classify_blastp_output/";
my $classify_blastn_output_dir = $project_dir."classify_blastn_output/";




if (defined $opt_d) {
	if (-d $search_input_dir) {
		rmtree ($search_input_dir);
	}


	if (-d $search_hmm_output_dir) {
		rmtree ($search_hmm_output_dir);
	}
	if (-d $search_blastx_output_dir) {
		rmtree ($search_blastx_output_dir);
	}
	if (-d $search_blastp_output_dir) {
		rmtree ($search_blastp_output_dir);
	}
	if (-d $search_blastn_output_dir) {
		rmtree ($search_blastn_output_dir);
	}
	
	
	if (-d $classify_hmm_input_dir) {
		rmtree ($classify_hmm_input_dir);
	}
	if (-d $classify_blastx_input_dir) {
		rmtree ($classify_blastx_input_dir);
	}
	if (-d $classify_blastp_input_dir) {
		rmtree ($classify_blastp_input_dir);
	}
	if (-d $classify_blastn_input_dir) {
		rmtree ($classify_blastn_input_dir);
	}
	
	
	if (-d $classify_hmm_output_dir) {
		rmtree ($classify_hmm_output_dir);
	}
	if (-d $classify_blastx_output_dir) {
		rmtree ($classify_blastx_output_dir);
	}
	if (-d $classify_blastp_output_dir) {
		rmtree ($classify_blastp_output_dir);
	}
	
	if (-d $classify_blastn_output_dir) {
		rmtree ($classify_blastn_output_dir);
	}
	
}


print STDERR "CARMA pipeline finished successfully!\n";
mailme("CARMA pipeline finished successfully!");











