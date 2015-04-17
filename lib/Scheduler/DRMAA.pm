package Scheduler::DRMAA;

# ----------------------------- COPYRIGHT -------------------------------
#
#	Copyright (C) 2008 CeBiTec, Bielefeld University
#					 <bioinfo@cebitec.uni-bielefeld.de>
#
#	 This file is part of CeBiTec Common Perl Libraries.
#
#	CeBiTec Common Perl Libraries is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	CeBiTec Common Perl Libraries is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with CeBiTec Common Perl Libraries; see the file LICENSE.  
#	If not, see <http://www.gnu.org/licenses/>.
#
# -----------------------------------------------------------------------

=head1 NAME

Scheduler::DRMAA

=head1 DESCRIPTION

scheduling of GenDB jobs using a DRMAA-compliant interface

=cut

use strict;
use warnings;
use Scheduler qw(:JOB_FLAGS);
use base qw(Scheduler);

use Schedule::DRMAAc qw(:all);
use IO::File;
use Carp qw(confess);

# the amount of seconds to sleep between creating the input
# file and submitting the job to allow the NFS server to sync
use constant NFS_SYNC_TIME => 30;

our ($available, $global_instance, $initialized);
$available = undef;
$global_instance = undef;
$initialized = undef;

1;

sub support_asynchronous {
	return 1;
}

sub _initialize {
	my ($self, $keep_session) = @_;
	my ($error, $diag) = drmaa_init(undef);
	if ($error) {
	$self->log(_drmaa_error("unable to initialize DRMAA library, disabling DRMAA scheduler", $error, $diag)) if ($error);
	return 0;
	}
	else {
	$initialized = 1;
	#$self->log("DRMAA library initialized");
	($error, my $impl, $diag) = drmaa_get_DRMAA_implementation();
	if ($error) {
		#$self->log(_drmaa_error("no information about implementation available: ",$error, $diag));
	}
	else {
		#$self->log("using $impl");
	}
	_drmaa_exit() unless ($keep_session);
	return 1;
	}
}

sub _drmaa_exit {
	if ($initialized) {
	my ($error, $diag) = drmaa_exit();
	print STDERR _drmaa_error("error while shutting down DRMAA session: ", $error, $diag) if ($error);
	$initialized = 0;
	}
	$available = undef;
}

# helper methods to log DRMAA error messages

sub _drmaa_error {
	my ($message, $error, $diag) = @_;
	if ($error =~ /^\d+$/) {
	$error = drmaa_strerror($error);
	}
	return $message.": $error".($diag ? ', diagnostic: '.$diag : "");
}

sub log_and_die {
	my ($self, $message, $error, $diag) = @_;
	if ($error) {
	   $message = _drmaa_error($message, $error, $diag);
	   $self->log($message."\n");
	   confess $message;
	}
}

sub available {
	my ($self) = @_;
	unless (defined ($available)) {
	$available = _initialize($self,0);
	}
	return $available;
}
	
sub new {
	my ($class, %parameters) = @_;
	if (defined ($global_instance)) {
	die "creating second instance of DRMAA scheduler is not permitted";
	}
	my $self = $class->SUPER::new(%parameters);
	bless $self, ref $class || $class;
	
	die unless ($self->_initialize(1));
	# create a job template
	my ($error, $jt, $diag) = drmaa_allocate_job_template();
	$self->log_and_die('allocating job template failed', $error, $diag);
	$self->{job_template} = $jt;
	$self->{job_ids} = [];
	$self->{job_info} = {};
	$self->{job_to_seq} = {};
	$global_instance = $self;
	$available = 1;
	return $self;
}

sub DESTROY {
	my ($self) = @_;
	
	# deallocate template
	if ($self->{job_template}) {
	my ($error, $diag) = drmaa_delete_job_template($self->{job_template});
	$self->log(_drmaa_error('deleting job template failed', $error, $diag)) if ($error);
	$self->{job_template} = undef;
	}
	_drmaa_exit() if ($available);
	$global_instance = undef;
	$available = undef;
}

sub END {
	_drmaa_exit() if ($available);
}

sub set_native_option {
	my ($self, $value) = @_;
	my ($error, $diag) = drmaa_set_attribute($self->{job_template}, 
						 $DRMAA_NATIVE_SPECIFICATION, 
						 $value);
	$self->log_and_die('error setting native option attribute', $error, $diag);
}

sub _fill_job_template {
	my ($self) = @_;

	# default values
	my ($error, $diag) = drmaa_set_attribute($self->{job_template}, $DRMAA_WD, $self->{temp_dir});
	$self->log_and_die('error setting working directory attribute', $error, $diag);
	
	foreach (keys %{$self->{options}}) {
		if (ref($self->{options}->{$_}) eq 'ARRAY') {
			($error, $diag) = drmaa_set_vector_attribute($self->{job_template},
								 $_, $self->{options}->{$_});
			$self->log_and_die('error setting vector attribute $_', $error, $diag);
		}
		else {
			($error, $diag) = drmaa_set_attribute($self->{job_template},
							  $_, $self->{options}->{$_});
			$self->log_and_die('error setting attribute $_', $error, $diag);
		}
	}
}


sub _substitute_string {
	my ($match, $subst) = @_;
	# escape meta characters in match and subst
	$match =~ s/\$/\\\$/g;
	$subst =~ s/\$/\\\$/g;
	return "s!$match!$subst!g";
}

sub submit2 {
	my ($self, $counter, $use_output_files) = @_;
	
	#print "counter: ".$counter."\n";
	
	$self->_fill_job_template();
	#my @names = (keys %{$self->{sequences}}, keys %{$self->{data_sets}});

	my @names = (1 .. $counter);

	if ($use_output_files == 1) {
		my $i = 1;
		foreach (@names) {
		
			$self->{job_info}->{$_}->{output} = $self->{temp_dir}."/$i.stdout";	
			$i++;
		}

		drmaa_set_attribute($self->{job_template}, $DRMAA_OUTPUT_PATH,
				':'.$self->{temp_dir}.'/'.$DRMAA_PLACEHOLDER_INCR.'.stdout');

		$self->{error_file}  = $self->{temp_dir} .$Scheduler::DRMAA::DRMAA_PLACEHOLDER_INCR.".stderr";

		if ($self->{error_file}) {
			# redirect STDERR to the given file
			#print "DRMAA_ERROR_PATH: ".$self->{error_file}."\n"; 
			drmaa_set_attribute($self->{job_template}, $DRMAA_ERROR_PATH,
					':'.$self->{error_file});
		}
		else {
			#print "DRMAA_ERROR_PATH: "."/dev/null"."\n";
			drmaa_set_attribute($self->{job_template}, $DRMAA_ERROR_PATH,
					':/dev/null');
		}
	} else {
	
		my $i = 1;
		foreach (@names) {
		
			$self->{job_info}->{$_}->{output} = "/dev/null";	
			$i++;
		}

		drmaa_set_attribute($self->{job_template}, $DRMAA_OUTPUT_PATH,
				':/dev/null');

		$self->{error_file}  = "/dev/null";

		if ($self->{error_file}) {
			# redirect STDERR to the given file
			#print "DRMAA_ERROR_PATH: ".$self->{error_file}."\n"; 
			drmaa_set_attribute($self->{job_template}, $DRMAA_ERROR_PATH,
					':'.$self->{error_file});
		}
		else {
			#print "DRMAA_ERROR_PATH: "."/dev/null"."\n";
			drmaa_set_attribute($self->{job_template}, $DRMAA_ERROR_PATH,
					':/dev/null');
		}
	}
	# set the command line
	my ($command, @arguments) = @{$self->{commandline}};
	
	drmaa_set_attribute($self->{job_template}, $DRMAA_REMOTE_COMMAND,
			$command);
 
 
	@arguments = () unless (@arguments);
	if (scalar(@arguments)) {
	# substitute the input/output placeholders in the arguments
	if (defined ($self->{input_placeholder})) {
		@arguments = map { ($_ && $_ eq $self->{input_placeholder}) ? 
				   '/dev/fd/0' : $_} @arguments;
	}
	if (defined ($self->{output_placeholder})) {
		# set the output place holder to /dev/fd/1 aka stdout
		@arguments = map { ($_ && $_ eq $self->{output_placeholder}) ? 
				   '/dev/fd/1' : $_} @arguments;
	}
	drmaa_set_vector_attribute($self->{job_template}, $DRMAA_V_ARGV, \@arguments);
	}

	# sleep to allow the NFS server to sync the files to the
	# cluster hosts
	sleep(NFS_SYNC_TIME);
 
	# submit the jobs
	# since the DRMAA master is overloaded from time to time,
	# we need to retry several times the submit fails
	my $retries = 3;
	while ($retries > 0) {
	my ($error, $job_iterator, $diag) = drmaa_run_bulk_jobs($self->{job_template},
								1, $counter, 1);
	$retries--;
	if ($error) {
		if ($error == $DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE) {
		# we had an error contacting the qmaster...
		if ($retries > 0) {
			# retry it
			$self->log('submitting bulk job failed with communication error, retrying after 10 seconds: '._drmaa_error($error, $diag)."\n");
			sleep(10);
			next;
		}
		}
		# give up
		$self->log_and_die('submitting bulk job failed: ',$error, $diag);
	} else {
		print "seems to have submitted..\n";
	}

	foreach(@names) {
		my ($error, $job_id) = drmaa_get_next_job_id($job_iterator);
		if ($error) {
			last;
		}
		$self->{job_info}->{$_}->{job_id} = $job_id;
		$self->{job_to_seq}->{$job_id} = $_;
		push @{$self->{job_ids}}, $job_id;
	}
	drmaa_release_job_ids($job_iterator);
	
	return @{$self->{job_ids}} if (wantarray);
	return $self->{job_ids};
	}
}



sub _drmaa_state_to_job_state {
	my ($state) = @_;
	# the flags defined in the scheduler class to not know
	# about DRMAA job state code
	# we have to transform the codes
	if ($state == $DRMAA_PS_UNDETERMINED) {
	return JOB_UNKNOWN;
	}
	elsif ($state == $DRMAA_PS_QUEUED_ACTIVE) {
	return JOB_QUEUED;
	}
	elsif ($state == $DRMAA_PS_SYSTEM_ON_HOLD) {
	return JOB_SYSTEM_HOLD;
	}
	elsif ($state == $DRMAA_PS_USER_ON_HOLD) {
	return JOB_USER_HOLD;
	}
	elsif ($state == $DRMAA_PS_USER_SYSTEM_ON_HOLD) {
	return JOB_SYSTEM_USER_HOLD;
	}
	elsif ($state == $DRMAA_PS_RUNNING) {
	return JOB_RUNNING;
	}
	elsif ($state == $DRMAA_PS_SYSTEM_SUSPENDED) {
	return JOB_SYSTEM_SUSPENDED;
	}
	elsif ($state == $DRMAA_PS_USER_SUSPENDED) {
	return JOB_USER_SUSPENDED;
	}
	elsif ($state == $DRMAA_PS_USER_SYSTEM_SUSPENDED) {
	return JOB_SYSTEM_USER_SUSPENDED;
	}
	elsif ($state == $DRMAA_PS_DONE) {
	return JOB_FINISHED;
	}
	elsif ($state == $DRMAA_PS_FAILED) {
	return JOB_FAILED;
	}
	else {
	return undef;
	}
}

sub get_status {
	my ($self, $job_id) = @_;
	
	my @states;
	my @jobs = (defined ($job_id) ? ($job_id) : @{$self->{job_ids}});
	foreach (@jobs) {
	if (defined ($self->{job_states}->{$_})) {
		push @states, $self->{job_states}->{$_};
		next;
	}
	my ($error, $ps, $diag) = drmaa_job_ps($_);
	if ($error) {
		$self->log(_drmaa_error("error retrieving job status for job $_",
					$error, $diag)) if ($error);
		push @states, JOB_UNKNOWN;
		next;
	}
	my $state = _drmaa_state_to_job_state($ps);
	if (defined ($state)) {
		push @states, $state;
	}
	else {
		$self->log("unknown log state $ps");
		push @states, JOB_UNKNOWN;
	}
	}
	if (defined ($job_id)) {
	return $states[0];
	}
	return @states if (wantarray);
	return \@states;
}

sub wait_for {
	my ($self, $job_id) = @_;
	my @jobs = ($job_id ) ? $job_id : @{$self->{job_ids}};
	
	while (scalar(@jobs)) {
	my ($error, $diag) = drmaa_synchronize(\@jobs,
						   $DRMAA_TIMEOUT_WAIT_FOREVER, 0);
	$self->log_and_die("error waiting for jobs: ",$error, $diag) if ($error);

	# synchronize should do the job, but we need to check the
	# single jobs again to reap their exit status
	for(my $i =0; $i < scalar @jobs;) {
		my ($error, $ps, $diag) = drmaa_job_ps($jobs[$i]);
		my $job_state = _drmaa_state_to_job_state($ps);
		if (defined ($job_state)) {
		if ((($job_state == JOB_FINISHED) ||
			 ($job_state == JOB_FAILED))) {
			$self->{job_states}->{$jobs[$i]} = $job_state;
			splice @jobs, $i, 1;
			next;
		}
		$i++;
		}
		else {
		$self->log("unknown log state $ps");
		$i++;
		}
		
	}
	}
}

sub cancel {
	my ($self, $job_id) = @_;
	my @jobs = (defined ($job_id) ? ($job_id) : @{$self->{job_ids}});
	foreach (@jobs) {
	my ($error, $diag) = drmaa_control($_, $DRMAA_CONTROL_TERMINATE);
	$self->log_and_die("error canceling job $_", $error, $diag);
	}
}

sub get_output {
	my ($self, $name) = @_;
	return undef unless ($self->{job_info}->{$name});
	unless (defined ($self->{job_info}->{$name}->{output})) {
	$self->log('undefined output file for sequence $name, not submitted yet ?');
	return undef;
	}
	my $handle = IO::File->new($self->{job_info}->{$name}->{output}, '<');
	return $handle;
}

sub get_job_id {
	my ($self, $name) = @_;
	return undef unless (defined ($self->{job_info}->{$name}));
	return $self->{job_info}->{$name}->{job_id};
}

sub dispose {
	my ($self, $seq_name) = @_;
	die "calling dispose() without submitted jobs" unless (defined($self->{job_info}));
	return unless (defined ($self->{job_info}->{$seq_name}));

	# remove input and output files
	foreach ($self->{job_info}->{$seq_name}->{input}, 
		 $self->{job_info}->{$seq_name}->{output}) {
	unlink $_ if (-f $_);
	}

	# cleanup internal management infos
	delete $self->{job_to_seq}->{$self->{job_info}->{$seq_name}->{job_id}};
	delete $self->{job_info}->{$seq_name};
}


sub dispose2 {
	my ($self, $seq_name) = @_;
	die "calling dispose() without submitted jobs" unless (defined($self->{job_info}));
	return unless (defined ($self->{job_info}->{$seq_name}));

	# remove input and output files
	#foreach ($self->{job_info}->{$seq_name}->{input}, 
#		 $self->{job_info}->{$seq_name}->{output}) {
#	unlink $_ if (-f $_);
 #   }

	# cleanup internal management infos
	delete $self->{job_to_seq}->{$self->{job_info}->{$seq_name}->{job_id}};
	delete $self->{job_info}->{$seq_name};
}


sub iterate2 {
	my ($self, $success_cb, $error_cb) = @_;
	die "no success callback given" unless (ref($success_cb) eq 'CODE');
	die "no error callback given" unless (ref($error_cb) eq 'CODE');

	while (scalar (keys %{$self->{job_info}})) {
	# wait for next job to be finished
	my ($error, $job_id_out, $stat, $rusage, $diag) =
		drmaa_wait($DRMAA_JOB_IDS_SESSION_ANY, 
			   $DRMAA_TIMEOUT_WAIT_FOREVER);
	if ($error) {
		if ($error == $DRMAA_ERRNO_INVALID_JOB) {
			# our jobs have been removed (e.g. by qdel)
			foreach (values (%{$self->{job_to_seq}})) {
				# set jobs to failed
				&$error_cb($_, 'job was removed during/before running on the cluster');
			}
			$self->log_and_die("error while waiting for jobs: jobs have been removed, aborting iteration.", $error, $diag);
		}
		elsif ($error == $DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE) {
			$self->log(_drmaa_error("error while waiting for jobs: ", $error, $diag));
			$self->log('waiting for 5 seconds and retrying');
			sleep(5);
			next;
		}
		elsif ($error == $DRMAA_ERRNO_EXIT_TIMEOUT) {
			$self->log('got timeout during wainting for jobs, retrying');
			sleep(2);
			next;
		}
		else {
			#$self->log_and_die("error while waiting for jobs: ", $error, $diag);
			$self->log('critical error? continue.. ;)');
		}
	}
	# we don't need the usage information
	drmaa_release_attr_values( $rusage );
	my $seq_name = $self->{job_to_seq}->{$job_id_out};
	#my $seq_name = "seq_name";
	my ($my_error, $aborted, undef) = drmaa_wifaborted($stat);
	
	
	my $exit_status;
	if ($aborted) {
		my $diagnosis;
		my $core_dumped;
		my $aborted;
		my $exited;
		my $signaled;
		my $termsig;
		print "stat: $stat\n";
		( $error, $core_dumped, $diagnosis ) = drmaa_wcoredump( $stat );
		print "core_dumped: $core_dumped\n";
		( $error, $aborted, $diagnosis ) = drmaa_wifaborted( $stat );
		print "aborted: $aborted\n";
		( $error, $exit_status, $diagnosis ) = drmaa_wexitstatus( $stat );
		print "exit_status: $exit_status\n";
		( $error, $exited, $diagnosis ) = drmaa_wifexited( $stat );
		print "exited: $exited\n";
		( $error, $signaled, $diagnosis ) = drmaa_wifsignaled( $stat );
		print "signaled: $signaled\n";
		( $error, $termsig, $diagnosis ) = drmaa_wtermsig( $stat );
		print "termsig: $termsig\n";
		
		&$error_cb($seq_name,'job failed before running on the cluster');
		print "error: $my_error\n";
		(undef,$exit_status,undef) = drmaa_wexitstatus($stat);
		print "job terminated with exit status $exit_status\n";
	} 
	else {
		my (undef, $exited, undef) = drmaa_wifexited($stat);
		if ($exited) {
		(undef,$exit_status,undef) = drmaa_wexitstatus($stat);
		if ($exit_status == 0) {
			# job run was ok
			#&$success_cb($seq_name, $self->get_output($seq_name));
			&$success_cb($seq_name);
		}
		else {
			&$error_cb($seq_name, "job terminated with exit status $exit_status");
		}			
		} 
		else {
		#$self->log(sprintf ("job exit status: %x\n",$stat));
		my (undef, $signaled, undef) = drmaa_wifsignaled($stat);
		if ($signaled) {
			my (undef, $termsig, undef) = drmaa_wtermsig($stat);
			&$error_cb($seq_name, "job finished due to signal $termsig");
		}
		else {
			&$error_cb($seq_name, 'job finished with unclear conditions');
		}
		}
	}
	$self->dispose2($seq_name);
	}
}

sub can_cancel {
	return 1;
}
