#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Gene_obj;
use Gene_obj_indexer;
use GTF_utils;
use Fasta_retriever;
use Carp;

my $usage = "\n\nusage: $0 gtf_file genome_db seqType=[cDNA|CDS|prot|CDSplus]\n\n";

my $gtf_file = $ARGV[0] or die $usage;
my $fasta_db = $ARGV[1] or die $usage;
my $seqType = $ARGV[2] or die $usage;

unless ($seqType =~ /^(cDNA|CDS|prot|CDSplus)$/) {
    die "Error, don't recognize seqType: $seqType, must be: cDNA|CDS|prot  ";
}


my $gene_obj_indexer = {};
    
## associate gene identifiers with contig id's.
&GTF_utils::index_GTF_gene_objs($gtf_file, $gene_obj_indexer);


my $fasta_retriever = new Fasta_retriever($fasta_db);


## associate all gene_ids with contigs
my @all_gene_ids = keys %$gene_obj_indexer;

my %contig_to_gene_list;
foreach my $gene_id (@all_gene_ids) {
    my $gene_obj = $gene_obj_indexer->{$gene_id};
    
    my $contig = $gene_obj->{asmbl_id} 
    or croak "Error, can't find contig id attached to gene_obj($gene_id) as asmbl_id val\n" 
        . $gene_obj->toString();
    
    
    my $gene_list_aref = $contig_to_gene_list{$contig};
    unless (ref $gene_list_aref) {
        $gene_list_aref = $contig_to_gene_list{$contig} = [];
    }

    push (@$gene_list_aref, $gene_id);

}



foreach my $asmbl_id (sort keys %contig_to_gene_list) {
    
    my $genome_seq = $fasta_retriever->get_seq($asmbl_id);
    
    my @gene_ids = @{$contig_to_gene_list{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer->{$gene_id};

        $gene_obj_ref->create_all_sequence_types(\$genome_seq);
        
        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
            
            my $isoform_id = $isoform->{Model_feat_name};
            my $gene_id = $isoform->{TU_feat_name};
            
            my $seq = "";
            if ($seqType eq 'cDNA') {
                $seq = $isoform->get_cDNA_sequence();
            }
            elsif ($seqType =~ /CDS/) {
                if ($isoform->is_coding_gene()) {
                    $seq = $isoform->get_CDS_sequence();
                }
                elsif ($seqType eq 'CDSplus') {
                    $seq = $isoform->get_cDNA_sequence();
                }
            }
            elsif ($seqType eq 'prot') {
                if ($isoform->is_coding_gene()) {
                    # protein
                    $seq = $isoform->get_protein_sequence();
                }
            }
            else {
                die "Error, shouldn't reach here...";
            }
            
            unless($seq) { next; }
            
            $seq =~ s/(\S{60})/$1\n/g; # make fasta format
            chomp $seq;
            
            my $com_name = $isoform->{com_name} || "";
            
            my $gene_name = $isoform->{gene_name};
            my $header = ">$isoform_id $gene_id";
            if ($gene_name) {
                $header .= " $gene_name ";
            }
            if ($com_name) {
                $gene_name .= " $com_name";
            }
            
            print "$header\n$seq\n";
        }
    }
}


exit(0);

