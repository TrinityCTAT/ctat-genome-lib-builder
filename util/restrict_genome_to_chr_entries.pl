#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\tusage: $0 genome.fa\n\n";

my $genome_fasta_file = $ARGV[0] or die $usage;

open(my $fh, $genome_fasta_file) or die "Error, cannot open file: $genome_fasta_file";

my $print_flag = 0;
while(my $line = <$fh>) {
    
    if ($line =~ /^>/) {
        if ($line =~ />chr/i) {
            $print_flag = 1;
        }
        else {
            $print_flag = 0;
        }
    }

    if ($print_flag) {
        print $line;
    }
}

exit(0);


    
    
