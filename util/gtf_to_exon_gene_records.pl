#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 annot.gtf\n\n";

my $gtf_file = $ARGV[0] or die $usage;

main: {

    open(my $fh, $gtf_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        my @x = split(/\t/);

        unless ($x[2] eq "exon") { next; }

        my $gene_id = "";
        if ($x[8] =~ /gene_id \"?([^\"]+)\"?;?/) {
            $gene_id = $1;
        }
        else {
            next;
        }

        if ($x[8] =~ /gene_name \"?([^\"]+)\"?;?/) {
            my $gene_name = $1;
            $gene_id = "$gene_name^$gene_id";
        }

        $x[8] = "gene_id \"$gene_id\"";
        
        print join("\t", @x) . "\n";
    }

    close $fh;

    exit(0);
}

