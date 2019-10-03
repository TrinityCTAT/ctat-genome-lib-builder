#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use TiedHash;
use Carp;
use Overlap_piler;
use Data::Dumper;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use JSON::XS;
require "$FindBin::Bin/../lib/isoform_to_genome_coord_conversion.ph";


my $usage = <<__EOUSAGE;

###########################################################
#
# --blast_outfmt6 <string>      blast pairs 
#
# --gtf <string>                genome annotation gtf file
#
###########################################################

__EOUSAGE

    ;

my $blast_file;
my $gtf_file;

my $help_flag;


&GetOptions ( 'h' => \$help_flag,
              'blast_outfmt6=s' => \$blast_file,
              'gtf=s' => \$gtf_file,
    );


unless ($blast_file && $gtf_file) {
    die $usage;
}




main: {

    my %isoform_info = &parse_isoform_gtf($gtf_file);

    

    print STDERR "-starting feature mapping.\n";
    
    my %gene_pair_to_coordsets;
        
    open(my $fh, $blast_file) or die "Error, cannot open file: $blast_file";
    
    my $num_blast_hits = `cat $blast_file | wc -l`;
    chomp $num_blast_hits;

    print STDERR "-num alignments: $num_blast_hits\n";
    
    my $counter = 0;
    while (<$fh>) {

        $counter++;
        if ($counter % 1000 == 0) {
            my $pct_done = sprintf("%.2f", $counter / $num_blast_hits * 100);
            print STDERR "\r[$pct_done % done]   ";
        }
        
        chomp;
        
        my $line = $_;
        
        my @x = split(/\t/);
        my $isoform_A = $x[0];
        my $isoform_B = $x[1];

        if ($isoform_A eq $isoform_B) { next; } # no selfies

        my $per_id = $x[2];
        my $E_value = $x[10];

        
        my ($lend_A, $rend_A) = sort {$a<=>$b} ($x[6], $x[7]);
        my ($lend_B, $rend_B) = sort {$a<=>$b} ($x[8], $x[9]);

        my $struct_A = $isoform_info{$isoform_A} or die "Error, no isoform struct for [$isoform_A]";
        my $struct_B = $isoform_info{$isoform_B} or die "Error, no isoform struct for [$isoform_B]";

        my $gene_A = $struct_A->{gene_id};
        my $gene_B = $struct_B->{gene_id};
        if ($gene_A eq $gene_B) { next; } # no selfies

        my $genome_A_lend = &translate_to_genomic_coord($struct_A, $lend_A);
        my $genome_A_rend = &translate_to_genomic_coord($struct_A, $rend_A);

        
        my $genome_B_lend = &translate_to_genomic_coord($struct_B, $lend_B);
        my $genome_B_rend = &translate_to_genomic_coord($struct_B, $rend_B);


        my $chr_A = $struct_A->{chr};
        my $chr_B = $struct_B->{chr};
        
        if ($gene_A gt $gene_B) {

            # swap everything for simplicity
            
            ($isoform_A, $gene_A, $chr_A, $genome_A_lend, $genome_A_rend,
             $isoform_B, $gene_B, $chr_B, $genome_B_lend, $genome_B_rend)  = 
                 
                ($isoform_B, $gene_B, $chr_B, $genome_B_lend, $genome_B_rend,
                 $isoform_A, $gene_A, $chr_A, $genome_A_lend, $genome_A_rend);
                
        }

        print join("\t",
                   $isoform_A, $gene_A, $chr_A, $genome_A_lend, $genome_A_rend,
                   $isoform_B, $gene_B, $chr_B, $genome_B_lend, $genome_B_rend,
                   $per_id, $E_value) . "\n";
        
                
                
    }
        
    print STDERR "-done\n";
    
        
    exit(0);
}

