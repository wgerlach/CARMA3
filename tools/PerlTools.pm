package PerlTools;


1;


sub init_CONSIDERED_RANKS_HASH() {
	my $CONSIDERED_RANKS_ref = shift(@_);
	my @CONSIDERED_RANKS = @_;
	
	
	@{$CONSIDERED_RANKS_ref} = @CONSIDERED_RANKS;

	
#	split(/,/,$optr);
	
	my %hash=();
		
	foreach my $rank (@CONSIDERED_RANKS) {
	#print "RANK $rank\n";
			$hash{lc($rank)}=1;
		}
	
	unless (%hash) {
		die;
		}

	return \%hash;
}


sub init_NCBI_NODES(){
	my ($NODES_FILE, $MERGED_FILE) = @_;
	my $NCBI_NODES;
	
	open (STREAM, $NODES_FILE) or die "cant open file ". $NODES_FILE ."\n";
	
	while (my $line = <STREAM>) {
		#print $line;
		my ($ncbi_id, $parent_id, $rank);
		unless ( ($ncbi_id, $parent_id, $rank)=split(/\s*\|\s*/, $line) ) {
			die;
		}
		#print "(".$ncbi_id.",".$parent_id.",".$rank.")\n";
		$NCBI_NODES->{$ncbi_id}=[$parent_id, $rank];
	}

	close(STREAM);


	open (STREAM2, $MERGED_FILE) or die "cant open file ". $MERGED_FILE ."\n";

	while (my $line = <STREAM2>) {
		#print $line;
		my ($old_id, $new_id);
		unless ( ($old_id, $new_id) =split(/\s*\|\s*/, $line) ) {
			die;
		}
		
		if (defined $NCBI_NODES->{$new_id} ) {
		
			unless (defined $NCBI_NODES->{$old_id} ) {
				
				#print "(".$old_id.",".$new_id.")\n";
				
				$NCBI_NODES->{$old_id}=$NCBI_NODES->{$new_id};
			}
		} 
		
		
		
	}
	
	
	close(STREAM2);	


	unless ($NCBI_NODES) {
		die;
	}
	
	return $NCBI_NODES;
}


sub extractTaxID{
	my ($d_ref, $NAMES_TO_TAXID) = @_;
	
	
	my $tax_name;
	(undef, $tax_name, undef, undef) = $$d_ref =~ /^>(\S+) (.+)/;
	
	if (not defined $tax_name) {
		print STDERR "error: tax_name undefined\n";
		exit(1);
	}

	my $semi = index($tax_name, ";");
	$tax_name = substr($tax_name, 0, $semi);
	

	$tax_name = lc($tax_name);
	
	#print "d: " . $$d_ref."\n";
	#print "tax_name: " .$tax_name."\n";
	
	my $taxid;
	
	while (not defined $NAMES_TO_TAXID->{$tax_name}) {
		
		my $pos = rindex($tax_name," ");
		
		if ($pos == -1) {
			print STDERR "error: pos == -1\n";
			die;
		}
		$tax_name = substr($tax_name, 0, $pos);
		
		#print "\"".$tax_name."\"\n";
		
		if (defined $NAMES_TO_TAXID->{$tax_name}) {
			last;
			#$taxid = $NAMES_TO_TAXID->{$tax_name};
		} else {
			#print STDERR "error: NOT defined NAMES_TO_TAXID->{tax_name}\n";
			#die;
		}
		
		#die;
	}
	
	$taxid = $NAMES_TO_TAXID->{$tax_name};
	return $taxid;
	
}

sub getRanks() {

			my ($ncbi_code, $NCBI_NODES, $NCBI_NAMES, $CONSIDERED_RANKS_HASH)=@_;
		
			my %ranks=();
		
			if ($ncbi_code==0) {
				return \%ranks;
			}
			
			while ($ncbi_code!=1) {
		
				
				unless (defined @{$NCBI_NODES->{$ncbi_code}}) {
					print STDERR "ncbi_code \"$ncbi_code\" not found in hash ncbi_nodes\n";
					print STDERR "Maybe you need a new taxdump.tar.gz ?\n";
					return 0;
					exit(1);
				}
				
				unless ( ($parent, $rank) = @{$NCBI_NODES->{$ncbi_code}} ) {
					die;
				}

				if ($CONSIDERED_RANKS_HASH->{$rank}==1) {
					$ranks{$rank}=$NCBI_NAMES->{$ncbi_code};
				}
				
				$ncbi_code= $parent;
			}
			return \%ranks;
}


sub getRanksTaxId() {
	
	my ($ncbi_code, $NCBI_NODES, $NCBI_NAMES, $CONSIDERED_RANKS_HASH)=@_;
	
	my %ranks=();
	
	if ($ncbi_code==0) {
		return \%ranks;
	}
	
	while ($ncbi_code!=1) {
		
		
		unless (defined @{$NCBI_NODES->{$ncbi_code}}) {
			print STDERR "Got nothing for ncbi_code in ncbi_nodes: \"$ncbi_code\"\n";
			print STDERR "Maybe you need a new taxdump.tar.gz ?\n";
			return 0;
			exit(1);
		}
		
		unless ( ($parent, $rank) = @{$NCBI_NODES->{$ncbi_code}} ) {
			die;
		}
		
		if (defined($CONSIDERED_RANKS_HASH->{$rank}) && $CONSIDERED_RANKS_HASH->{$rank}==1) {
			$ranks{$rank}=$ncbi_code;
		}
		
		$ncbi_code= $parent;
	}
	return \%ranks;
}

# maps name to taxid (warning: not unique!)
sub init_NCBI_NAMES_TO_TAXID(){
	my ($NAMES_FILE) = @_;
	my $NCBI_NAMES_TO_TAXID;
	
	open (STREAM, $NAMES_FILE);
	
	while (my $line = <STREAM>) {
		#print $line;
		my ($ncbi_id, $name, $nametype);
		unless (($ncbi_id, $name, undef, undef) =split(/\s*\|\s*/, $line) ) {
			die;
		}
		#print "GOT: [$ncbi_id, $name, $nametype]\n";
		#if ($nametype eq "scientific name") {
		
#		if (defined $NCBI_NAMES_TO_TAXID->{lc($name)}) {
#		
#			if ($NCBI_NAMES_TO_TAXID->{lc($name)} != $ncbi_id) {
#				print STDERR "already defined: $name and different tax_id!\n";
#				print STDERR "$line\n";
#				print STDERR "continue...\n";
#				#exit(1);
#			}
#		}
		
		if ($name =~ /^environmental samples$/) {
			next;
		}

		# here we check for non unique names.
		if (defined $NCBI_NAMES_TO_TAXID->{lc($name)}) {
			my $i = 2;

			while (defined $NCBI_NAMES_TO_TAXID->{lc($name)."|".$i}) {
				$i++;
			}
			
			#print STDERR "insert: \"".lc($name)."|".$i."\"\n";
			$NCBI_NAMES_TO_TAXID->{lc($name)."|".$i}=$ncbi_id;
		} else {

			$NCBI_NAMES_TO_TAXID->{lc($name)}=$ncbi_id;
		}
		
		# e.g.:  "Rickettsia prowazekii" -> "R. prowazekii"
		
		#print "W ".$name."\n";
		my @parts = split( / / , $name );
		my ($init) = $parts[0] =~ /^(\S)/;
		
		
		if ((defined $init) && scalar(@parts) > 1 ) {
			$parts[0] = $init.".";
			
			my $variant = join(" ", @parts);
			#print $variant1."\n";
			unless (defined $NCBI_NAMES_TO_TAXID->{lc($variant)}) {
				$NCBI_NAMES_TO_TAXID->{lc($variant)} = $ncbi_id;
			} 
			substr($variant, 2, 1) = "";
			#print "v2 ".$variant1."\n";
			#exit(1);
			unless (defined $NCBI_NAMES_TO_TAXID->{lc($variant)}) {
				$NCBI_NAMES_TO_TAXID->{lc($variant)} = $ncbi_id;
			}
		
		}
		
			#$NCBI_NAMES->{$ncbi_id}=lc($name);	
			#print $name."\n";
		#}
	}

	close(STREAM);
	return $NCBI_NAMES_TO_TAXID;
}


# maps taxid to name
sub init_NCBI_NAMES(){

	my ($NAMES_FILE, $MERGED_FILE) = @_;
	my $NCBI_NAMES;

	open (STREAM, $NAMES_FILE) or die "cant open file ". $NAMES_FILE ."\n";
	
	while (my $line = <STREAM>) {
		#print $line;
		my ($ncbi_id, $name, $nametype);
		unless (($ncbi_id, $name, undef, $nametype) =split(/\s*\|\s*/, $line) ) {
			die;
		}
		#print "GOT: [$ncbi_id, $name, $nametype]\n";
		if ($nametype eq "scientific name") {
			$NCBI_NAMES->{$ncbi_id}=$name;	
			#print $name."\n";
		}
	}

	close(STREAM);

	unless ($NCBI_NAMES) {
		die;
	}
	
	# include MERGED_FILE
	
	
	open (STREAM2, $MERGED_FILE) or die "cant open file ". $MERGED_FILE ."\n";

	while (my $line = <STREAM2>) {
		#print "LINE!!!! ".$line;
		my ($old_id, $new_id);
		unless (($old_id, $new_id) =split(/\s*\|\s*/, $line) ) {
			die;
		}
		#print "LINE: $line";
		#print "mapping: $old_id, $new_id\n";
		if (defined $NCBI_NAMES->{$new_id} ) {
				#print STDERR "new $new_id exists\n";		
			unless (defined $NCBI_NAMES->{$old_id} ) {
				#print "(".$old_id.",".$new_id.")\n"; 
				$NCBI_NAMES->{$old_id}=$NCBI_NAMES->{$new_id};
			}
		} 
		
		
		
	}
	
	
	close(STREAM2);	
	
	return $NCBI_NAMES;
}

