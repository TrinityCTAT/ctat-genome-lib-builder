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

my $usage = <<__EOUSAGE;

###########################################################
#
# --dfam_coords <string>         dfam genomic coords file (geneated by: isoform_dfam_gene_chr_conversion.pl)
#
# --out_prefix                  output file prefix
#
###########################################################

__EOUSAGE

    ;

my $dfam_coords;
my $out_prefix;

my $help_flag;


&GetOptions ( 'h' => \$help_flag,
              'dfam_coords=s' => \$dfam_coords,
              'out_prefix=s' => \$out_prefix,
    );


unless ($dfam_coords && $out_prefix) {
    die $usage;
}


main: {
        
    my %gene_to_coordsets;
    
    open(my $fh, $dfam_coords) or die "Error, cannot open file: $dfam_coords";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($gene, $genome_lend, $genome_rend) = ($x[3], $x[5], $x[6]);
        
        push (@{$gene_to_coordsets{$gene}}, [$genome_lend, $genome_rend]);

    }
    close $fh;
    

    
    &create_dfam_region_db(\%gene_to_coordsets, $out_prefix);
    
    print STDERR "\n\nDone.\n";
    

    exit(0);
}


####
sub create_dfam_region_db {
    my ($gene_to_coordsets_href, $out_prefix) = @_;
    
    my $genome_collapsed_coords_file = "$out_prefix.genomic_regions.dat";
    open(my $ofh, ">$genome_collapsed_coords_file") or die "Error, cannot write to file: $genome_collapsed_coords_file";
    my $idx = new TiedHash( { create => "$out_prefix.genomic_regions.dbm" } );
        
    foreach my $gene (keys %{$gene_to_coordsets_href}) {
        
        my @coords = @{$gene_to_coordsets_href->{$gene}};
        
        my @ranges = &Overlap_piler::simple_coordsets_collapser(@coords);

        my $data_store_struct = { 
            'coords' => \@ranges,
        };
        
        my $json = &encode_json($data_store_struct);

        # print "$json\n";

        $idx->store_key_value($gene, $json);
        
        print $ofh join("\t", $gene, $json) . "\n";
        
    }

    close $ofh;
      
    return;
}


####
sub set_rel_coords {
    my ($struct) = @_;

    my @exons = sort {$a->{lend} <=> $b->{lend}} @{$struct->{exons}};
                      
    if ($struct->{orient} eq '-') {
        @exons = reverse @exons;
    }

    my $last_rel_pos = 0;
    
    foreach my $exon (@exons) {

        my ($lend, $rend) = ($exon->{lend}, $exon->{rend});
        
        my $exon_len = $rend - $lend + 1;

        $exon->{rel_lend} = $last_rel_pos + 1;
        $exon->{rel_rend} = $last_rel_pos + $exon_len;

        $last_rel_pos += $exon_len;
    }

    return;
}
    
