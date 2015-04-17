package FileUtils::TempDir;

=head1 NAME

FileUtils::TempDir

=head1 DESCRIPTION

wrapper class for temporary directories.

This module mimics File::Temp, providing control over a temporary directory. As soon as a FileUtils::TempDir object leaves the scope, the destructor of this class remove the temporary directory and all contained files.

This module is not intended to be used as a stand-alone module. Use the get_tmp_dir_object() method defined in L<FileUtils> instead !

=head1 Methods

=over 4

=item * FileUtils::TempDir->new(parameters)

creates a new FileUtils::TempDir objects. Parameters are passed to the tempdir() function. See L<File::Temp> for a description. The constructor also ensures that the UNLINK flag is set to zero.

=item * $tempdir->dirname()

return the name of the temporary directory

=item * $tempdir->tempfile(parameters)

uses the tempfile() function from L<File::Temp> to create a temporary file in the temporay directory.

=back

=cut

    use strict;
use warnings;
use File::Temp;
use File::Path qw(rmtree);

1;

sub new {
    my ($class,$template, %parameters) = @_;

    # ensure that the directory is not unlinked by other means than
    # our destructor...
    $parameters{CLEANUP}=0;
    my $dirname = File::Temp::tempdir ($template, %parameters);
    my $self = {dirname => $dirname};
    bless $self, ref $class || $class;
    return $self;
}

sub dirname {
    return $_[0]->{dirname};
}

sub tempfile {
    my $self = shift;
    
    my $template;
    my %parameters;
    unless (ref $_[0]) {
	$template = shift;
    }
    %parameters = @_;
    $parameters{DIR} = $self->dirname();
    $parameters{UNLINK} = 0;
    
    if ($template) {
	return File::Temp::tempfile($template, %parameters);
    }
    else {
	return File::Temp::tempfile(%parameters);
    }
}

sub DESTROY {
    my ($self) = @_;
    rmtree ($self->dirname);
}
