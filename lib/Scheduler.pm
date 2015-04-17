package Scheduler;

# ----------------------------- COPYRIGHT -------------------------------
#
#    Copyright (C) 2008 CeBiTec, Bielefeld University
#                     <bioinfo@cebitec.uni-bielefeld.de>
#
#     This file is part of CeBiTec Common Perl Libraries.
#
#    CeBiTec Common Perl Libraries is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    CeBiTec Common Perl Libraries is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with CeBiTec Common Perl Libraries; see the file LICENSE.  
#    If not, see <http://www.gnu.org/licenses/>.
#
# -----------------------------------------------------------------------

=head1 NAME

Scheduler

=head1 DESCRIPTION

abstract superclass and factory for scheduler implementations used with BRIDGE

=head1 Factory methods

=over 4

=item * get_scheduler_class(tool)

given a GenDB tool, this factory method returns the preferred scheduler for the tool or undef if no scheduler is available.

=item * job_based()

returns true if the scheduler works based on jobs instead of sequences. These schedulers manages jobs and sequences internally and their interface differ from sequence based schedulers. Default to sequence based scheduling. Scheduler classes need to overwrite this method to change the behaviour.

=back

=head1 Methods

=over 4

=item * available()

returns true if a scheduler class is available on a system

=item * new(options)

create a new scheduler instance. An instance encapsulates all relevant information about a job that is do be scheduled, e.g. command line, environment variables, scheduler options etc. See section about setting up jobs below. If creating a new scheduler instance fails, the new() method die with an error message.

=back

=head1 Sequence scheduler

A sequence scheduler handles handles raw input sequences like DNA or AA sequences. These schedulers may support two different modes of operation, synchronous and asynchronous scheduling. All implementation of scheduler classes have to support synchronous scheduling. The capability to support asynchronous scheduling can be queried with the support_asynchronous() method.

=head2 Common methods

=over 4 

=item * support_asynchronous()

returns true if the scheduler class supports asynchronous scheduling.

=item * add_input_sequence(name, sequence)

adds a input sequence to a job. The name has to be unique within a bulk job. Returns false if a name is already used, true otherwise.

=item * add_input_data(name, data)

similar to add_input_sequence, but the data may be an arbitary formatted string. Use this method pass multiple fasta file or binary data. Returns false if a name is already used, true otherwise.

=item * set_option(option name, option value)

sets a scheduler options. Option names and possible values vary from scheduler class to scheduler. See documentation about the supported classes.

=item * submit()

submits the current jobs and returns the an array containing the job ids. If an error occurs during submission, the method dies with an error message.

=item * get_status(job id)

returns the status for the current job. If no job_id is given, this method returns the status of all jobs submitted in the current session.

=item * get_output(name)

returns a filehandle to the output of a job. If no job id is given, this method returns an array of filehandle, one for each job.

=item * cancel(job_id)

cancels the current job and remove it from the queuing system. If an error occurs during the call, the method dies with an error message.

=item * dispose(sequence name)

clean up temporary files that were created by the scheduler and removes the job information from the scheduler objects. This method should be called after a job terminated and the job output was processed, e.g. by calling get_output(). This method is called internally by the asynchronous scheduling methods.

=back

=head2 Synchronous scheduling

"Synchronous" refers to the way job processing and job analysis are made. With synchronous calls job arrays are handled in a monolithic way.

=over 4

=item * wait_for(job_id)

waits until the current job terminates. If no job_id is given, this method waits until all jobs submitted in the current session are terminated. If an error occurs during waiting for jobs, the method dies with an error message.

=item * get_output(name)

returns a filehandle to the output of a job. If no job id is given, this method returns an array of filehandle, one for each job.

=back

=head2 Asynchronous scheduling

Instead of process a job arrays as a whole this scheduling interface focusses on single jobs within an array. Every time a single job finishes its results are processed at once, in contrast to the bulk processing done with synchronous scheduling.

=over 4

=item * iterate(success_callback, error_callback)

waits until one of the submitted job finishes. If the job is finished successfully, the success_callback is invoked, using the sequence name and a filehandle tied to the job output as parameters. If the job has failed, the error callback is invoked, with the sequence name and an error message as arguments.

=back

=head1 Job based scheduler

=over 4

=item * add_jobs (@jobs)

adds the given jobs to the internal list of jobs to be processed

=item * execute($annotate_flag, $project_name)

executes the jobs. Depending on the kind of scheduler the jobs may be executed locally, one-by-one, or submitted to some queuing system.

=back

=head1 Job states

The scheduler class defines a number for states a job may be in. These states are modelled after the DRMAA specification and a scheduler class may not support all states.

=over 4

=item * JOB_UNKNOWN

information about the job is not available (e.g. invalid job id)

=item * JOB_QUEUED

the job has been submitted to the queuing system and is waiting for execution

=item * JOB_SYSTEM_HOLD

a system hold has been placed on the job

=item * JOB_USER_HOLD

a user hold has been placed on the job

=item * JOB_SYSTEM_USER_HOLD

a system and user hold has been placed on the job

=item * JOB_RUNNING

the job is currently being executed

=item * JOB_SYSTEM_SUSPENDED

job has been suspended by the system

=item * JOB_USER_SUSPENDED

job has been suspended by a user

=item * JOB_SYSTEM_USER_SUSPENDED

job has been suspended by the system and a user

=item * JOB_FINISHED

the job has been finished successfully

=item * JOB_FAILED

the job has been failed

=back

These flags are available by importing the JOB_FLAGS export tags.

=head1 Checking scheduler capabilities

Since some scheduler classes may not support certain functions, the following methods may be used to check for them:

=over 4

=item * can_cancel()

returns true if a scheduler class is able to cancel a job. The job may either be running or submitted.

=back

=head1 Internal methods

These methods are mostly convinient methods for use in the scheduler classes.

=over 4

=item * set_logger(code ref)

set a new logging method. Everytime a log message is written it is passed as single parameter to the given code ref. The default logging method print the message to STDERR.

=item * log(message)

log a message

=back

=head1 Setting up a job

The parameter of the new() method defines what kind of job is to be created and how this jobs is going to be executed. The parameters are passed as a hash with certain keys and value:

=over 4

=item * input_placeholder =E<gt> string

a variable in the command line that is substituted with the name of the input file.

=item * output_placefolder =E<gt> string

a placeholder for an output file in the command line.

=item * commandline =E<gt> string

the command to be executed.

=item * temp_dir =E<gt> string

directory for storing temporary files. Some scheduler classes require this directory to be readable across the network.

=item * error_file =E<gt> string

if given the output of STDERR will be redirected to this file, if the scheduler class supports error redirection.

=back

=cut

    use strict;
use warnings;
use base qw (Exporter);

# job status flags
use constant JOB_UNKNOWN => 0;
use constant JOB_QUEUED => 1;
use constant JOB_SYSTEM_HOLD => 2;
use constant JOB_USER_HOLD => 3;
use constant JOB_SYSTEM_USER_HOLD => 4;
use constant JOB_RUNNING => 5;
use constant JOB_SYSTEM_SUSPENDED => 6;
use constant JOB_USER_SUSPENDED => 7;
use constant JOB_SYSTEM_USER_SUSPENDED => 8;
use constant JOB_FINISHED => 9;
use constant JOB_FAILED => 10;

our(@EXPORT_OK);
@EXPORT_OK = qw(JOB_UNKNOWN JOB_QUEUED JOB_SYSTEM_HOLD JOB_USER_HOLD
		 JOB_SYSTEM_USER_HOLD JOB_RUNNING JOB_SYSTEM_SUSPENDED
		 JOB_USER_SUSPENDED JOB_SYSTEM_USER_SUSPENDED JOB_FINISHED
		 JOB_FAILED);
our (%EXPORT_TAGS);
%EXPORT_TAGS = (JOB_FLAGS => [qw(JOB_UNKNOWN JOB_QUEUED JOB_SYSTEM_HOLD 
				 JOB_USER_HOLD JOB_SYSTEM_USER_HOLD
				 JOB_RUNNING JOB_SYSTEM_SUSPENDED
				 JOB_USER_SUSPENDED JOB_SYSTEM_USER_SUSPENDED 
				 JOB_FINISHED JOB_FAILED)]);

1;

sub BEGIN {
    # create abstract methods
    foreach (qw(submit get_status wait_for cancel get_output 
		can_cancel)) {
	eval "sub $_ { die \"calling abstract method $_ in \".__PACKAGE__;}";
    }
}

# this class is abstract, so no instances are available...
sub available {
    return 0;
}

sub job_based {
    return 0;
}

sub support_asynchronous {
    return 0;
}

sub _logger {
    print STDERR map {"[$$]: $_\n"} split("\n",$_[0]);
}

sub log {
    my ($self, $message) = @_;
    if (ref ($self)) {
	&{$self->{logger}}($message);
    }
    else {
	_logger($message);
    }
}

sub new {
    my ($class, %parameters) = @_;
    
    # IC IC no error checking yet, add this later
    my $self = {logger => \&_logger,
		sequences => {},
		data_sets => {},
		jobs => [],
		options => {},
		commandline => $parameters{commandline},
		input_placeholder => $parameters{input_placeholder} || undef,
		output_placeholder => $parameters{output_placeholder} || undef,
		temp_dir => $parameters{temp_dir},
		error_file => $parameters{error_file}};
    
    bless $self, ref $class || $class;
    return $self;
}

sub add_input_sequence {
    my ($self, $name, $sequence) = @_;
    return 0 if ($self->{sequences}->{$name});
    return 0 if ($self->{data_sets}->{$name});
    $self->{sequences}->{$name} = $sequence;
}

sub add_input_data {
    my ($self, $name, $data) = @_;
    return 0 if ($self->{sequences}->{$name});
    return 0 if ($self->{data_sets}->{$name});
    $self->{data_sets}->{$name} = $data;
}

sub add_jobs {
    my ($self, @jobs) = @_;
    push @{$self->{jobs}}, @jobs;
}

sub set_option {
    my ($self, $name, $value) = @_;
    $self->{options}->{$name} = $value;
}

sub set_native_option {
    die "calling abstract method set_native_option in package ".__PACKAGE__;
}

sub set_logger {
    my ($self, $coderef) = @_;
    $self->{logger} = $coderef;
}
