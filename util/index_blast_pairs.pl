#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use TiedHash;
use Carp;
use DelimParser;
use List::Util qw(max);

my $usage = "\n\n\tusage: output_index_filename blast_pairs.gene_sym.gz [min_pct_overlap=25]\n\n";

my $output_index_filename = $ARGV[0] or die $usage;
my $blast_pairs_file = $ARGV[1] or die $usage;
my $min_pct_overlap = $ARGV[2];
if (! defined $min_pct_overlap) {
    $min_pct_overlap = 25;
}

main: {

    my $idx = new TiedHash( { create => $output_index_filename } );

    print STDERR "-processing blast pairs from: $blast_pairs_file\n";
    
    open (my $fh, "gunzip -c $blast_pairs_file |") or die "Error, cannot read $blast_pairs_file  ";
    my $tab_reader = new DelimParser::Reader($fh, "\t");
    
    my $num_stored = 0;
    my $num_skipped = 0;

    while (my $row = $tab_reader->get_row() ) {
        
        my $geneA = $tab_reader->get_row_val($row, "#query_gene_acc");
        my $geneB = $tab_reader->get_row_val($row, "target_gene_acc");

        my $per_id = $tab_reader->get_row_val($row, "avg_per_id");
        my $Evalue = $tab_reader->get_row_val($row, "min_Evalue");

        my $pct_query_len = $tab_reader->get_row_val($row, "pct_query_len");
        my $pct_target_len = $tab_reader->get_row_val($row, "pct_target_len");
        
        my $max_per_len = max($pct_query_len, $pct_target_len);
        if ($max_per_len < $min_pct_overlap) {
            $num_skipped += 1;
            next; 
        }

        my $token = "{$geneA|$geneB|pID:$per_id|E:$Evalue}";
        
        $idx->store_key_value("$geneA--$geneB", $token);
        $idx->store_key_value("$geneB--$geneA", $token);
        
        $num_stored += 1;
    }
    
    
    $idx->store_key_value("geneABC--geneXYZ", "__placeholder_testval__"); # for sanity checking

    print STDERR "-stored $num_stored blast pairs;  skipped $num_skipped as having < $min_pct_overlap align pct overlap.\n";
    
    print STDERR "-done creating index file: $output_index_filename\n\n";
    
    
    exit(0);
    
}

