#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../lib");
use Pipeliner;

my $usage = "\n\n\tusage: $0 /path/to/ctat_genome_lib_build_dir\n\n";


## expects you're starting with a plug-n-play library 


my $ctat_genome_lib = $ARGV[0] or die $usage;


my $UTILDIR = "$FindBin::Bin";


main: {

    chdir($ctat_genome_lib) or die "Error, cannot cd to $ctat_genome_lib";
    
    
    my @required_inputs = ("blast_pairs.dat.gz",
                           "trans.blast.dat.gz",
                           "PFAM.domtblout.dat.gz",
                           "ref_genome.fa",
                           "ref_annot.gtf",
                           "ref_annot.gtf.gene_spans");
    
    # verify we have what we need.
    foreach my $required_file (@required_inputs) {
        unless (-s $required_file) {
            die "Error, not finding required file at $ctat_genome_lib named $required_file";
        }
    }
    

    ## rebuild fusion annot db index
    my $cmd = "$UTILDIR/build_fusion_annot_db_index.pl --gene_spans ref_annot.gtf.gene_spans --out_db_file fusion_annot_lib.idx";
    if (-s "fusion_annot_lib.gz") {
        $cmd .= " --key_pairs fusion_annot_lib.gz ";
    }

    &Pipeliner::process_cmd($cmd);
    
    
    ## index blast homologous pairs
    $cmd = "$UTILDIR/index_blast_pairs.pl blast_pairs.idx blast_pairs.dat.gz";
    &Pipeliner::process_cmd($cmd);
    
    ## index cdna pairwise alignments
    $cmd = "$UTILDIR/build_chr_gene_alignment_index.pl --blast_genes trans.blast.dat.gz  --out_prefix trans.blast.align_coords";
    &Pipeliner::process_cmd($cmd);

    $cmd = "$UTILDIR/index_pfam_domain_info.pl --pfam_domains PFAM.domtblout.dat.gz --genome_lib_dir . ";
    &Pipeliner::process_cmd($cmd);

    print STDERR "\n\n\t-done, all indices rebuilt\n\n";
    
    exit(0);
}


### $idx->store_key_value("geneABC--geneXYZ", "__placeholder_testval__");
