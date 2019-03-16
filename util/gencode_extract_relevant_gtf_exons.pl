#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 gencode.annotation.gtf\n\n";

my $gtf_file = $ARGV[0] or die $usage;

main: {

    open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
    while (<$fh>) {
        my $line = $_;
        if (/^\#/) { next; }
        my @x = split(/\t/);

        # just exon features
        unless ($x[2] eq "exon") { next; }

        my $info = $x[8];
        if ($info =~ /gene_type \"([^\"]+)\"/) {
            my $gene_type = $1;

            if ($gene_type =~ /^(protein_coding|lincRNA)$/) { 
                print $line;
            }
        }
    }
    close $fh;

    exit(0);
}

