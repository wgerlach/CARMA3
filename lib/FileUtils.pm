package FileUtils;

# several methods to work with files and directories

# $Id: FileUtils.pm,v 1.1 2004/11/05 13:27:48 blinke Exp $

# $Log: FileUtils.pm,v $
# Revision 1.1  2004/11/05 13:27:48  blinke
# Initial revision
#
# Revision 1.11  2004/09/14 15:55:51  blinke
# fixed problem with choosing the wrong directory for temporary files in get_tmp_file*
# added get_tmp_dir_object
#
# Revision 1.10  2004/09/14 09:04:48  blinke
# fixed a problem with get_tmp_dir creating directories in the wrong location
#
# Revision 1.9  2004/08/24 13:38:02  blinke
# added get_tmp_dir
#
# Revision 1.8  2004/08/24 12:32:05  blinke
# rewrote get_tmp_file and get_tmp_file_object from scratch
#
# Revision 1.7  2004/08/23 16:29:46  blinke
# code cleanup
#
# Revision 1.6  2004/07/26 11:37:37  blinke
# undeleted after accidental removal
#
# Revision 1.4  2003/12/10 17:44:51  alice
# *** empty log message ***
#
# Revision 1.3  2003/12/05 13:14:37  alice
# added method get_tmp_file
#
# Revision 1.2  2003/05/12 11:00:10  blinke
# removed error in remove_dir
#
# Revision 1.1  2002/03/26 16:16:50  blinke
# moved to common modules
#
# Revision 1.2  2002/03/15 15:18:57  blinke
# added documentation
#
# Revision 1.1.1.1  2002/03/06 09:02:27  agoesman
# initial CVS import for  GENDB 2.0
#
# Revision 1.2  2002/02/26 16:20:51  blinke
# added create_dir
#
# Revision 1.1  2002/02/22 13:52:52  blinke
# Initial revision
#


=head1 NAME

FileUtils

=head1 DESCRIPTION

Collection of simple functions to handle files and directories

=head2 Synopsis

# create the directories a, a/b and a/b/c
FileUtils::create_dir("a/b/c");

# remove a complete directory hierarchie
FileUtils::remove_dir("a");

=head2 Defined functions

=over 4

=item * create_dir (dir name, permission mask)   [DEPRECATED]

Creates a sub directory. Each element of the path is created if it does not exists. The permission mask defaults to 0755 and is applied to each newly created directory.

=item * remove_dir (dir name)                    [DEPRECATED]

Removes a directory, including all files and sub directories.

=item * get_tmp_file (directory [, suffix])

creates a temporary file. In scalar context, this function returns a file handle to the file. In list context, it returns the file handle and the file name. An optional suffix is appended to the file name upon creation. Be aware that this may introduce race conditions.

Temporary files created with this function are not unlinked at program termination. The programmer is responsible for handling this !

The function is based on the L<File::Temp> module. For hints how to use the temporary file properly read the L<File::Temp> documentation !

=item * get_tmp_file_object(directory [,suffix])

similar to get_tmp_file, this method creates a temporary file encapsulated in a L<File::Temp> object. See the L<File::Temp> manpage for information about these objects. In contrast to the temporary files create by get_tmp_file, the temporary file is removed at object destruction.

=item * get_tmp_dir(directory)

creates a temporary sub directory in the given directory and returns the complete path. The directory is not removed automatically.

=back

=cut

use strict;
use warnings;
use File::Temp qw(tempfile tempdir);
use File::Spec;
use File::Basename;
use FileUtils::TempDir;
use Carp qw(croak);
use base qw(Exporter);

our @EXPORT = qw(create_dir remove_dir
		 get_tmp_file get_tmp_file_object
		 get_tmp_dir get_tmp_dir_object );

use constant TEMPFILE_TEMPLATE => 'XXXXXXXX';

1;

# recursibly create a directory hierarchie

sub create_dir {
    my ($dirname, $mode) = @_;

    warn "create_dir in module ".__PACKAGE__." is deprecated\n";
    if (!defined ($mode)) {
	$mode = 0775;
    }
    my $dir = "";
    foreach (split /\//, $dirname) {
	$dir = $dir."/$_";
	if (! -d $dir) {
	    mkdir $dir, $mode || croak "Cannot create directory $dir: $@";
	}
    }
}

# recursivly remove a directory
sub remove_dir {
    my ($dirname) = @_;

    warn "remove_dir in module ".__PACKAGE__." is deprecated\n";
    if (!opendir (DIR, $dirname)) {
	print STDERR "unable to open directory $dirname in remove_dir()\n";
	return;
    };
    for my $file (readdir (DIR)) {
        next if ($file =~ /^(\.|\.\.)$/); # skip pseudo directories . and ..
        if (-d $dirname."/".$file) {
            remove_dir ($dirname."/".$file);
        }
        else {
            unlink ($dirname."/".$file);
        }
    }
    closedir (DIR);
    rmdir ($dirname);
}

sub get_tmp_file {
    my ($dirname, $suffix) = @_;
    defined($suffix) or $suffix = '';

    # due to problems with the automounter 
    # have have to stat() the directory we want to create
    # temporary files in
    stat($dirname.'/foobar') if (defined ($dirname));

    # get the handle and the filename
    my ($tmphandle, $tmpname) = (defined $dirname) ?
	tempfile (TEMPFILE_TEMPLATE,
		  DIR => $dirname,
		  UNLINK => 0,
		  SUFFIX => $suffix) :
		      tempfile (TEMPFILE_TEMPLATE,
			    UNLINK => 0,
				SUFFIX => $suffix);
    if (wantarray) {
	return ($tmphandle, $tmpname);
    }
    else {
	return ($tmpname);
    }
}

sub get_tmp_file_object{
    my ($dirname, $suffix) = @_;
    defined($suffix) or $suffix = '';
    
    # due to problems with the automounter 
    # have have to stat() the directory we want to create
    # temporary files in
    stat($dirname.'/foobar') if (defined ($dirname));

    if ($dirname) {
	return File::Temp->new(TEMPLATE => TEMPFILE_TEMPLATE,
			       DIR => $dirname,
			       SUFFIX => $suffix);
    }
    else {
	return File::Temp->new(TEMPLATE => TEMPFILE_TEMPLATE,
			       SUFFIX => $suffix,
			       DIR => File::Spec->tmpdir);
    }
}

sub get_tmp_dir {
    my ($dirname) = @_;
    
    if (defined ($dirname)) {
	stat($dirname.'/foobar');
	return tempdir(TEMPFILE_TEMPLATE,
		       DIR => $dirname, CLEANUP => 0);
    }
    else {
	return tempdir(TEMPFILE_TEMPLATE,
		       TMPDIR => 1, CLEANUP => 0);
    }
}

sub get_tmp_dir_object {
    my ($dirname) = @_;
    if (defined ($dirname)) {
	stat($dirname.'/foobar');
	return FileUtils::TempDir->new(TEMPFILE_TEMPLATE,
				       DIR => $dirname);
    }
    else {
	return FileUtils::TempDir->new(TEMPFILE_TEMPLATE,
				       DIR => File::Spec->tmpdir);
    }
}
