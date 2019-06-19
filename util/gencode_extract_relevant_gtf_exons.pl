#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\n\tusage: $0 gencode.annotation.gtf\n\n";

my $gtf_file = $ARGV[0] or die $usage;



    #cat gencode.v19.annotation.gtf.noreadthrus | perl -lane 'if (/gene_type \"([^\"]+)\".*transcript_type \"([^\"]+)/) { print "$1\t$2";}' | token_counter.pl 1 | tee feature_types.txt
my %WANT = (

    'protein_coding' => 1,
    
    #protein_coding	protein_coding	1848709.00
    #protein_coding	nonsense_mediated_decay	280806.00
    #protein_coding	processed_transcript	137779.00
    #protein_coding	retained_intron	132311.00
    

    'lincRNA' => 1,
    
    #lincRNA	lincRNA	51518.00
    
    #antisense	antisense	41237.00

    'pseudogene' => 1,
    
    #pseudogene	processed_pseudogene	22972.00
    #pseudogene	pseudogene	15239.00
    #pseudogene	unprocessed_pseudogene	11202.00
    #pseudogene	processed_transcript	10912.00
    #miRNA	miRNA	9165.00
    
    'miRNA' => 1,
    
    #pseudogene	transcribed_unprocessed_pseudogene	7086.00
    

    'misc_RNA' => 1,
    
    #misc_RNA	misc_RNA	6102.00
    #snRNA	snRNA	5748.00
    #processed_transcript	processed_transcript	5493.00

    #snoRNA	snoRNA	4371.00
    #sense_intronic	sense_intronic	3145.00
    #processed_transcript	lincRNA	3047.00
    #processed_transcript	antisense	2960.00
    #processed_transcript	retained_intron	1854.00
    #polymorphic_pseudogene	polymorphic_pseudogene	1718.00
    #rRNA	rRNA	1581.00
    #pseudogene	unitary_pseudogene	1430.00
    #sense_overlapping	sense_overlapping	1350.00
    
    'IGV_V_gene' => 1,
    #IG_V_gene	IG_V_gene	1123.00
    
    #pseudogene	transcribed_processed_pseudogene	1091.00
    #protein_coding	non_stop_decay	1070.00
    #pseudogene	retained_intron	1050.00
    #polymorphic_pseudogene	protein_coding	1017.00
    
    'TR_V_gene' => 1,
    
    #TR_V_gene	TR_V_gene	763.00
    #IG_V_pseudogene	IG_V_pseudogene	681.00
    #polymorphic_pseudogene	processed_transcript	426.00
    
    'TR_J_gene' => 1,
    
    #TR_J_gene	TR_J_gene	300.00
    #lincRNA	retained_intron	190.00
    
    'IG_C_gene' => 1,
    
    #IG_C_gene	IG_C_gene	185.00
    #polymorphic_pseudogene	retained_intron	184.00
    #antisense	retained_intron	173.00
    #polymorphic_pseudogene	nonsense_mediated_decay	170.00
    
    'IG_D_gene' => 1,
    
    #IG_D_gene	IG_D_gene	152.00
    #3prime_overlapping_ncrna	3prime_overlapping_ncrna	100.00
    #TR_V_pseudogene	TR_V_pseudogene	99.00
    
    'IG_J_gene' => 1,
    
    #IG_J_gene	IG_J_gene	82.00
    #processed_transcript	snoRNA	74.00
    #processed_transcript	sense_overlapping	67.00
    #Mt_tRNA	Mt_tRNA	66.00
    #lincRNA	miRNA	66.00
    
    'TR_C_gene' => 1,
    
    #TR_C_gene	TR_C_gene	58.00
    #lincRNA	snoRNA	42.00
    #IG_C_pseudogene	IG_C_pseudogene	36.00
    #processed_transcript	miRNA	26.00
    #antisense	miRNA	24.00
    #lincRNA	misc_RNA	20.00
    #lincRNA	processed_transcript	18.00
    #antisense	snoRNA	14.00
    #sense_intronic	snoRNA	14.00
    #antisense	misc_RNA	12.00
    #TR_J_pseudogene	TR_J_pseudogene	12.00
    #lincRNA	snRNA	12.00
    
    'TR_D_gene' => 1,
    
    #TR_D_gene	TR_D_gene	12.00
    #sense_intronic	retained_intron	10.00
    #IG_J_pseudogene	IG_J_pseudogene	9.00
    #lincRNA	rRNA	8.00
    #Mt_rRNA	Mt_rRNA	6.00
    #sense_intronic	miRNA	4.00
    #pseudogene	translated_processed_pseudogene	3.00
    #processed_transcript	sense_intronic	3.00
    #sense_overlapping	miRNA	2.00
    #sense_intronic	snRNA	2.00
    
    );



main: {

    open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
    while (<$fh>) {
        my $line = $_;
        if (/^\#/) { next; }
        my @x = split(/\t/);

        unless ($x[2] =~ /^(exon|CDS)$/) { next; }
        
        my $info = $x[8];
        if ($info =~ /gene_type \"([^\"]+)\"/) {
            my $gene_type = $1;

            if ($WANT{$gene_type}) {
                print $line;
            }
        }
    }
    close $fh;

    exit(0);
}

