#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use TiedHash;
use Carp;

my $usage = "\n\n\tusage: output_index_filename blast_pairs.gene_sym.gz [blast_pairs_supp.gene_sym.gz ...]\n\n";

my $output_index_filename = $ARGV[0] or die $usage;
shift @ARGV;
my @blast_pairs_files = @ARGV;

unless (@blast_pairs_files) {
    die $usage;
}


main: {

    my $idx = new TiedHash( { create => $output_index_filename } );

    my %best_bit_score;
    
    foreach my $blast_pairs (@blast_pairs_files) {

        print STDERR "-processing blast pairs from: $blast_pairs\n";
        
        open (my $fh, "gunzip -c $blast_pairs |") or die "Error, cannot read $blast_pairs  ";
        
        my $total_linecount = `gunzip -c $blast_pairs | wc -l `;
        chomp $total_linecount;
        $total_linecount = int($total_linecount);
        
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my ($geneA, $geneB, $per_id, $Evalue, $score) = ($x[0], $x[1], $x[2], $x[10], $x[11]);
                
            if ($geneA eq $geneB) { next; }
            
            my $gene_pair_tok = join("^", sort ($geneA, $geneB));
            
            if ( (! exists $best_bit_score{$gene_pair_tok}) || $best_bit_score{$gene_pair_tok} < $score) {

                $best_bit_score{$gene_pair_tok} = $score;
                
            
                my $token = "{$geneA|$geneB|pID:$per_id|E:$Evalue}";
                
                $idx->store_key_value("$geneA--$geneB", $token);
                $idx->store_key_value("$geneB--$geneA", $token);
                
            }
        }
    }


    $idx->store_key_value("geneABC--geneXYZ", "__placeholder_testval__"); # for sanity checking
    
    print STDERR "-done creating index file: $output_index_filename\n\n";

    exit(0);
    
}

