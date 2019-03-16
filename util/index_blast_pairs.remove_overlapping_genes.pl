#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Carp;
use Process_cmd;
use DB_File;

my $usage = "\n\n\tusage: $0 ctat_genome_lib_dir\n\n";

my $ctat_genome_lib_dir = $ARGV[0] or die $usage;


main: {

    my $blast_idx = "$ctat_genome_lib_dir/blast_pairs.idx";
    unless (-s $blast_idx) {
        die "Error, cannot locate file: $blast_idx";
    }

    # make a dated copy of the original:
    my $backup = $blast_idx . ".prev." . time();
    print STDERR "-making a backup of the original index at: $backup\n";
    &process_cmd("cp $blast_idx $backup");
        
    my $gene_spans_file = "$ctat_genome_lib_dir/ref_annot.gtf.gene_spans";
    unless (-s $gene_spans_file) {
        die "Error, cannot locate file: $gene_spans_file";
    }
        
    my %idx;
    
    tie (%idx, 'DB_File', $blast_idx, O_RDWR, 0666, $DB_BTREE);
    

    my %gene_to_span_info = &parse_gene_spans($gene_spans_file);

    my $counter = 0;
    foreach my $gene_pair (keys %idx) {
        $counter++;
        print STDERR "\r[$counter]   " if $counter % 100 == 0;
        my ($geneA, $geneB) = split(/--/, $gene_pair);
        if (&overlap($geneA, $geneB, \%gene_to_span_info)) {
            print STDERR "-pruning: $gene_pair\n";
            delete($idx{$gene_pair});
        }
    }
    
    untie(%idx);
        
    print STDERR "-done updating blast index\n\n";
    
    exit(0);
    
}


####
sub parse_gene_spans {
    my ($gene_spans_file) = @_;


    my %gene_spans;
    
    open(my $fh, $gene_spans_file) or die "Error, cannot open file: $gene_spans_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($chr, $lend, $rend, $gene_id) = ($x[1], $x[2], $x[3], $x[5]);

        ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
        
        $gene_spans{$gene_id} = [$chr, $lend, $rend];
    }
    close $fh;

    return(%gene_spans);

}


####
sub overlap {
    my ($geneA, $geneB, $gene_to_span_info_href) = @_;

    
    
    my $geneA_info_aref = $gene_to_span_info_href->{$geneA};

    unless ($geneA_info_aref) {
        print STDERR "-warning, no record of gene span for: $geneA, skipping pair($geneA,$geneB)\n";
        return(0);
    }
    
    my $geneB_info_aref = $gene_to_span_info_href->{$geneB};
    unless ($geneB_info_aref) {
        print STDERR "-warning, no record of gene span for: $geneB, skipping pair($geneA,$geneB)\n";
        return(0);
    }
        
    my ($chr_A, $lend_A, $rend_A) = @$geneA_info_aref;
        
    my ($chr_B, $lend_B, $rend_B) = @$geneB_info_aref;

    if ($chr_A eq $chr_B
        &&
        $lend_A < $rend_B && $rend_A > $lend_B) {
        return(1);
    }
    else {
        return(0);
    }
}

