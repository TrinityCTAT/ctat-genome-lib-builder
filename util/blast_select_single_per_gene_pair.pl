#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib";
use DelimParser;

use List::Util qw(max);

my $usage = "\n\tusage: $0 allvsall.outfmt6.wGeneSym.gz\n\n";


my $blast_info = $ARGV[0] or die $usage;

unless ($blast_info =~ /\.gz$/) {
    die "Error, input file must be gzipped";
}


main: {

    my %gene_pair_to_match_info;

    open(my $fh, "gunzip -c $blast_info | ") or die "Error, cannot open file: $blast_info";
    my $tab_reader = new DelimParser::Reader($fh, "\t");
    while(my $row = $tab_reader->get_row()) {
        
        my $geneA = $tab_reader->get_row_val($row, "#query_gene_acc");
        my $geneB = $tab_reader->get_row_val($row, "target_gene_acc");
        
        my $per_len_A = $tab_reader->get_row_val($row, "pct_query_len");
        my $per_len_B = $tab_reader->get_row_val($row, "pct_target_len");

        my $gene_pair_key = join("$;", sort($geneA, $geneB));
        
        my $max_per_len = max($per_len_A, $per_len_B);

        my $min_e_val = $tab_reader->get_row_val($row, "min_Evalue");

        push (@{$gene_pair_to_match_info{$gene_pair_key}}, { row => $row,
                                                             max_per_len => $max_per_len,
                                                             min_e_val  => $min_e_val } );
        
    }

    ## output single best per gene pair.
    
    my @column_headers = $tab_reader->get_column_headers();
    
    my $tab_writer = new DelimParser::Writer(*STDOUT, "\t", \@column_headers);

    foreach my $gene_pair (sort keys %gene_pair_to_match_info) {
        
        my @hits = @{$gene_pair_to_match_info{$gene_pair}};
        
        @hits = sort {$b->{max_per_len} <=> $a->{max_per_len}
                      ||
                          $a->{min_e_val} <=> $b->{min_e_val} } @hits;


        my $best_hit = shift @hits;

        $tab_writer->write_row($best_hit->{row});

    }


    exit(0);
}


