#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use File::Path;
use Digest::MD5;
use JSON::PP;
use Process_cmd;


my $usage = <<__EOUSAGE__;

##################################################################################################
#
#  Reverts a CTAT genome lib supplemented by supplement_ctat_genome_lib.pl to its original state:
#  restores each <file>.bkup over the modified file and removes <genome_lib_dir>/supplemental/
#
#  Required:
#
#   --genome_lib_dir <string>     path to ctat_genome_lib_build_dir
#
#  Optional:
#
#   --keep_supplement_files       move supplemental/ to supplemental.reverted.<time> instead of removing it
#
#   --force                       restore even if a .bkup file's checksum differs from the one recorded
#                                 when it was created
#
##################################################################################################


__EOUSAGE__

    ;


my $help_flag;
my $genome_lib_dir;
my $keep_supplement_files_flag = 0;
my $force_flag = 0;

&GetOptions ( 'h' => \$help_flag,
              'genome_lib_dir=s' => \$genome_lib_dir,
              'keep_supplement_files' => \$keep_supplement_files_flag,
              'force' => \$force_flag,
    );

if ($help_flag) {
    die $usage;
}

unless ($genome_lib_dir) {
    die $usage;
}

if (@ARGV) {
    die "Error, unrecognized parameters: @ARGV\n";
}

$genome_lib_dir = &ensure_full_path($genome_lib_dir);

my $SUPP_DIR = "$genome_lib_dir/supplemental";
my $MANIFEST = "$SUPP_DIR/manifest.json";
my $LOCK_DIR = "$SUPP_DIR/.lock";

my $have_lock = 0;

main: {

    unless (-e $MANIFEST) {
        my @stray_bkups = glob("$genome_lib_dir/*.bkup");
        if (@stray_bkups) {
            die "Error, no $MANIFEST found, but the lib contains .bkup files:\n  "
                . join("\n  ", @stray_bkups) . "\n"
                . "These were not recorded by supplement_ctat_genome_lib.pl; restore them manually if appropriate.\n";
        }
        die "Error, $genome_lib_dir has not been supplemented (no $MANIFEST)\n";
    }

    unless (mkdir($LOCK_DIR)) {
        die "Error, $LOCK_DIR exists: a supplement or revert may be running.\n"
            . "  If not, remove that directory and try again.\n";
    }
    $have_lock = 1;

    open(my $fh, $MANIFEST) or die "Error, cannot read $MANIFEST";
    my $manifest = JSON::PP->new->decode(do { local $/; <$fh> });
    close $fh;

    my %backups = %{$manifest->{backups}};

    ## verify every backup before restoring any of them
    my @problems;
    foreach my $file (sort keys %backups) {
        my $bkup = "$genome_lib_dir/$file.bkup";
        unless (-e $bkup) {
            push (@problems, "missing backup: $bkup");
            next;
        }
        print STDERR "-verifying $file.bkup\n";
        my $md5 = &md5_file($bkup);
        if ($md5 ne $backups{$file}) {
            push (@problems, "checksum of $bkup ($md5) differs from the one recorded at backup ($backups{$file})");
        }
    }
    if (@problems) {
        if ($force_flag && ! grep { /^missing/ } @problems) {
            print STDERR "-WARNING, proceeding due to --force:\n  " . join("\n  ", @problems) . "\n";
        }
        else {
            die "Error, cannot safely revert; no files changed:\n  " . join("\n  ", @problems) . "\n";
        }
    }

    foreach my $file (sort keys %backups) {
        print STDERR "-restoring original $file\n";
        rename("$genome_lib_dir/$file.bkup", "$genome_lib_dir/$file")
            or die "Error, cannot rename $genome_lib_dir/$file.bkup to $genome_lib_dir/$file";
    }

    rmdir($LOCK_DIR);
    $have_lock = 0;

    if ($keep_supplement_files_flag) {
        my $dest = "$SUPP_DIR.reverted." . time();
        rename($SUPP_DIR, $dest) or die "Error, cannot rename $SUPP_DIR to $dest";
        print STDERR "-supplement files moved to $dest\n";
    }
    else {
        rmtree($SUPP_DIR);
        print STDERR "-removed $SUPP_DIR\n";
    }

    print STDERR "\n-done. $genome_lib_dir restored to its original (pre-supplement) state.\n\n";

    exit(0);
}


END {
    if ($have_lock) {
        rmdir($LOCK_DIR);
    }
    ## consistent failure status (die otherwise exits with whatever errno was last set)
    if ($? != 0) {
        $? = 1;
    }
}


####
sub md5_file {
    my ($file) = @_;

    open(my $fh, $file) or die "Error, cannot open $file";
    binmode($fh);
    my $md5 = Digest::MD5->new->addfile($fh)->hexdigest;
    close $fh;

    return($md5);
}
