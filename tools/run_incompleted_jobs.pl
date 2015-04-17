#!/usr/bin/env perl
# Written by Wolfgang Gerlach

use strict;
use warnings;

unless (@ARGV) {
	print "usage: runincl oldprojectdir tempprojectdir [email]\n";
	exit(0);
}

my $old_project_dir = $ARGV[0];

my $temp_project_dir = $ARGV[1];


unless (-d $old_project_dir) {
	print SDTERR $old_project_dir." not found\n";
	exit(1);
	
}

unless (-d $temp_project_dir) {
	print SDTERR $temp_project_dir." not found\n";
	exit(1);
	
}


#print $old_project_dir."\n";

my $old_search_input_dir = $old_project_dir."/search_input/";

unless (-d $old_search_input_dir) {
	print SDTERR $old_search_input_dir." not found\n";
	exit(1);
	
}

my @old_search_input_dir_files = <$old_search_input_dir/chunk*.fas>; 

my @old_search_input_dir_files_number;

foreach my $file (@old_search_input_dir_files) {
	chomp($file);
	#print $file . "\n";
	
	my ($num) = $file =~ /chunk_(\d+).fas$/;
	#print $num."\n";
	push (@old_search_input_dir_files_number, $num);
	
}

print "found ". @old_search_input_dir_files_number ." fasta files\n";


my $old_search_blastx_output_dir = $old_project_dir."/search_blastx_output/";

unless (-d $old_search_blastx_output_dir) {
	print SDTERR $old_search_blastx_output_dir." not found\n";
	exit(1);
	
}


my @old_search_blastx_output_dir_files = <$old_search_blastx_output_dir/chunk_*.blastx>; 

my %old_search_blastx_output_dir_files_numberhash;

foreach my $file (@old_search_blastx_output_dir_files) {
	chomp($file);
	#print $file . "\n";
	
	my ($num) = $file =~ /chunk_(\d+).blastx$/;
	#print $num."\n";
	
	$old_search_blastx_output_dir_files_numberhash{$num}=1;
	#push (@old_search_input_dir_files_number, $num);
	
}

my @list_of_missing_jobs;

foreach my $num (@old_search_input_dir_files_number) {
	
	unless (defined $old_search_blastx_output_dir_files_numberhash{$num} ) {
		push(@list_of_missing_jobs, $num);
		
	}
	
}

print  @list_of_missing_jobs ." jobs are missing\n";

# -----------------------------------------------
my $total = @list_of_missing_jobs;
my $count = 0;

my $fail_count = 0;

#foreach my $num (@list_of_missing_jobs) {
for (my $i=0 ; $i < $total; $i++ ) {
	my $num = $list_of_missing_jobs[$i];
	$count++;
	
	print "process jobid $num. $count of $total\n";
	my $jobdir = $temp_project_dir."/chunk_".$num."/";
	print $jobdir."\n";
	if (-e $jobdir."blastx_result.tax") {
		next;
	}
	if (-e $jobdir) {
		# delete content
		system ("rm -rf ".$jobdir."*");
	} else {
		mkdir $jobdir;
	}
	my $inputfile = $old_search_input_dir . "chunk_".$num.".fas";
	my $command = "/vol/carma/CARMA3/carma3-SGE.pl -p $jobdir -a /vol/carma/CARMA3/carma.cfg -x -c 300";
	if (($i == $total - 1) && (defined $ARGV[2])) {
		$command .= " -e ".$ARGV[2];
	}
	$command .= " $inputfile";
		
	print $command."\n";
	system($command);
		
	
	unless (-e $jobdir."blastx_result.tax") {
		print STDERR "warning: ".$jobdir."blastx_result.tax not found\n";
		$fail_count++;
	}
	
}


if ($fail_count == 0) {
	print "finished. no fails. total jobs: $total\n";
} else {
	print "finished with ".$fail_count." fails. total jobs: $total\n";
	exit(1);
}



# check if tax files have been created...
for (my $i=0 ; $i < $total; $i++ ) {
	my $num = $list_of_missing_jobs[$i];
	my $jobdir = $temp_project_dir."/chunk_".$num."/";
	print $jobdir."\n";
	unless (-e $jobdir."blastx_result.tax") {
		
		print "file ".$jobdir."blastx_result.tax"." is missing\n";
		exit(1);
	}
}



for (my $i=0 ; $i < $total; $i++ ) {
	my $num = $list_of_missing_jobs[$i];
	#print $num."\n";
	my $jobdir = $temp_project_dir."/chunk_".$num."/";
	print $jobdir."\n";
	#unless (-e $jobdir."blastx_result.tax") {
	my $newblast = $old_project_dir."/search_blastx_output/chunk_".$num.".blastx";
	unless (-e $newblast) {
		my $command = "cat ".$temp_project_dir."/chunk_".$num."/search_blastx_output/chunk_*.blastx > $newblast";
		print $command."\n";
		system($command);
	}
	
	#my $newtaxfile = $old_project_dir."/classify_blastx_output/chunk_".$num.".tax";
	#print $newtaxfile."\n";
	#unless (-e $newtaxfile) {
#		my $command = "cp ".$temp_project_dir."/chunk_".$num."/blastx_result.tax $newtaxfile";
#		print $command."\n";
		#system($command);
		
#	}
	
	
}


