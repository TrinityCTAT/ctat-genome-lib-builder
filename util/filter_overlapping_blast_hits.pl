#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Carp;


my $usage = "\n\n\tusage: $0 gene_blast_pairs.gz ref_annot.gtf.gene_spans > gene_blast_pairs_overlaps_filt.gz\n\n";

my $blast_pairs_gz = $ARGV[0] or die $usage;
my $gene_spans_file = $ARGV[1] or die $usage;

main: {

    my %gene_to_span_info = &parse_gene_spans($gene_spans_file);

    my %gene_pair_overlap_info;

    my $count_removed = 0;

    open(my $fh, "gunzip -c $blast_pairs_gz | ") or die "Error reading $blast_pairs_gz";
    while(my $line = <$fh>) {
        my ($geneA, $transA, $geneB, @rest) = split(/\t/, $line);
        
        my $gene_pair = join("--", sort ($geneA, $geneB));
        my $overlaps = $gene_pair_overlap_info{$gene_pair};
        if (! defined $overlaps) {
            $overlaps = (&overlap($geneA, $geneB, \%gene_to_span_info)) ? "YES" : "NO";
            $gene_pair_overlap_info{$overlaps} = $overlaps;
        }
        if ($overlaps eq "YES") {
            $count_removed += 1;
        }
        else {
            print $line;
        }
    }
    close $fh;

    print STDERR "-removed $count_removed blast entries involving overlapping genes\n";

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

