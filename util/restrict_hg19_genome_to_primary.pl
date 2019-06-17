#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\tusage: $0 genome.fa\n\n";

my $genome_fasta_file = $ARGV[0] or die $usage;

open(my $fh, $genome_fasta_file) or die "Error, cannot open file: $genome_fasta_file";

my $print_flag = 0;
while(my $line = <$fh>) {
    
    if ($line =~ /^>/) {

        $print_flag = 0;
        
        if ($line =~ /^>(\S+) (\S+)/) {
            my $valA = $1;
            my $valB = $2;

            if ($valA =~ /^chr/i) {
                $print_flag = 1;
            }
            else {
            
                $valA =~ s/chr//i;
            
                if ($valA eq $valB) {
                    $print_flag = 1;
                    # not a patch
                }
            }
        }
    }
    
    if ($print_flag) {
        print $line;
    }
}

exit(0);

