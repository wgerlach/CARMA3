package SplitFasta;


1;



# returns number of files (chunks) produced

sub splitFasta{
	($numberOfChunks, $fasta_file, $output_dir, $prefix) = @_;
	
	my $CHUNK_FILE;

	$last_char = substr($output_dir,-1,1);
	if ($last_char ne "/") {
		$output_dir .= "/";
	}

	

	# count sequences:
	my $numberOfSequences = 0;

	open(FILE, '-|', "zcat -f '$fasta_file'")  or die $!;
	while (<FILE>) { 
		#print $_; 
		if ($_ =~ /^>/) {
			$numberOfSequences++;
		}
	}

	close (FILE);

	if ($numberOfSequences == 0) {
		print STDERR "splitFasta: no sequences found !?\n";
		print STDERR "input file was: ".$fasta_file."\n";
		return -1;
	}

	my $first_chunk_file = $output_dir.$prefix."1.fas";
	
	if (-e $first_chunk_file) {
		print STDERR "First chunk file $first_chunk_file already exists, so skip all...\n";
		print STDERR "To force creation of new chunks delete old chunk files first!\n";
		if ($numberOfSequences < $numberOfChunks) {
			return $numberOfSequences;
		} else {
			return $numberOfChunks; 
		}
	}

	my $sequencesPerChunk = int($numberOfSequences / $numberOfChunks);
	my $rest = $numberOfSequences - ($numberOfChunks * $sequencesPerChunk);

	print STDERR "splitFasta information:\n";
	print STDERR "numberOfSequences: $numberOfSequences\n";
	print STDERR "numberOfChunks: $numberOfChunks\n";
	print STDERR "sequencesPerChunk: $sequencesPerChunk\n";
	print STDERR "rest: $rest\n";

	# create chunks:
	my $chunk_id=1;
	my $sequences_in_chunk=0;
	if ($chunk_id <= $rest) {
		$sequences_in_chunk--;
		}

	my $chunk_file_name = $output_dir.$prefix.$chunk_id.".fas";
	
	
	
	open $CHUNK_FILE, ">>$chunk_file_name" or die $!;


	open(FILE, '-|', "zcat -f '$fasta_file'")  or die $!;
	while ($line =<FILE>)
	{

		if ($line =~ /^>/) {
	
			# check if a new chunk has to be opened:
			if ($sequences_in_chunk >= $sequencesPerChunk) {
				if (defined $CHUNK_FILE) {
					close($CHUNK_FILE);
				}
				$chunk_id++;
				$sequences_in_chunk=0;
				if ($chunk_id <= $rest) {
					$sequences_in_chunk--;
					}
				
				my $chunk_file_name = $output_dir.$prefix.$chunk_id.".fas";
				open $CHUNK_FILE, ">>$chunk_file_name" or die $!;
			}
		
			$sequences_in_chunk++;
		}
	
		unless (defined $CHUNK_FILE) {
			die;
		}
		print $CHUNK_FILE $line;

	}

	if (defined $CHUNK_FILE) {
		close($CHUNK_FILE);
	}
	close (FILE);

	if ($numberOfSequences < $numberOfChunks) {
		return $numberOfSequences;
	} else {
		return $numberOfChunks; 
	}

	
}
