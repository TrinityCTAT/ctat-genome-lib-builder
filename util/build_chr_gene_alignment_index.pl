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
# --blast_genes <string>        blast alignments in chr/gene format
#                                (must be sorted by gene pair! 
#                                  ex. sort -k2,2 -k7,7 blast.toGenes > blast.toGenes.sorted )
#
# --out_prefix                  output file prefix
#
###########################################################

__EOUSAGE

    ;

my $blast_file;
my $out_prefix;


my $help_flag;


&GetOptions ( 'h' => \$help_flag,
              'blast_genes=s' => \$blast_file,
              'out_prefix=s' => \$out_prefix,
    );


unless ($blast_file && $out_prefix) {
    die $usage;
}



my $genome_collapsed_coords_file = "$out_prefix.align_coords.dat";
open(my $ofh, ">$genome_collapsed_coords_file") or die "Error, cannot write to file: $genome_collapsed_coords_file";
my $idx = new TiedHash( { create => "$out_prefix.align_coords.dbm" } );


main: {

    
    my %gene_pair_to_coordsets;


    my $fh;
    my $num_blast_hits;
    if ($blast_file =~ /\.gz$/) {
        $num_blast_hits = `gunzip -c $blast_file | wc -l`;
        open($fh, "gunzip -c $blast_file | ") or die "Error, cannot open pipe to gunzip -c $blast_file";
    }
    else {
        $num_blast_hits = `cat $blast_file | wc -l`;
        open($fh, $blast_file) or die "Error, cannot open file: $blast_file";
    }
    chomp $num_blast_hits;

    unless ($num_blast_hits =~ /\d/ && $num_blast_hits > 0) {
        die "Error, no blast hits from file: $blast_file";
    }
    
    my %seen_pairs;
    
    my $prev_gene_pair = "";
    
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

        my ($gene_A, $genome_A_lend, $genome_A_rend,
            $gene_B, $genome_B_lend, $genome_B_rend) = ($x[1], $x[3], $x[4],
                                                        $x[6], $x[8], $x[9]);

        if ($gene_A eq $gene_B) { next; } # no selfies

        ($genome_A_lend, $genome_A_rend) = sort {$a<=>$b} ($genome_A_lend, $genome_A_rend);
        ($genome_B_lend, $genome_B_rend) = sort {$a<=>$b} ($genome_B_lend, $genome_B_rend);

        my $gene_pair = "$gene_A--$gene_B";
        
        if ($gene_pair ne $prev_gene_pair) {

            if ($seen_pairs{$gene_pair}) {
                confess "Error, already processed gene pair: [$gene_pair] ... are the data not sorted? ";
            }
            
            if (%gene_pair_to_coordsets) {
                &summarize_and_store_gene_pair_info(%gene_pair_to_coordsets);
            }
            %gene_pair_to_coordsets = (); # reinit
            $seen_pairs{$gene_pair} = 1;
        }
        
        push (@{$gene_pair_to_coordsets{$gene_pair}->{A_coords}}, [$genome_A_lend, $genome_A_rend]);
        push (@{$gene_pair_to_coordsets{$gene_pair}->{B_coords}}, [$genome_B_lend, $genome_B_rend]);
        

        $prev_gene_pair = $gene_pair;
    }

    close $fh;

    if (%gene_pair_to_coordsets) {
        &summarize_and_store_gene_pair_info(%gene_pair_to_coordsets);
    }



    close $ofh;
    
    print STDERR "\n\nDone.\n";
    

    exit(0);
}


####
sub summarize_and_store_gene_pair_info {
    my (%gene_pair_to_coordsets) = @_;
    
    ## now summarize and store results.
    
    foreach my $gene_pair (keys %gene_pair_to_coordsets) {

        my @A_coords = @{$gene_pair_to_coordsets{$gene_pair}->{A_coords}};
        my @B_coords = @{$gene_pair_to_coordsets{$gene_pair}->{B_coords}};

        my @A_ranges = &Overlap_piler::simple_coordsets_collapser(@A_coords);
        my @B_ranges = &Overlap_piler::simple_coordsets_collapser(@B_coords);

        #print "$gene_pair: " . Dumper(\@A_ranges) . Dumper(\@B_ranges);

        my ($gene_A, $gene_B) = split(/--/, $gene_pair); 

        my $data_store_struct = {  'gene_A' => $gene_A,
                                   'coords_A' => \@A_ranges,

                                   'gene_B' => $gene_B,
                                   'coords_B' => \@B_ranges,
        };

        my $json = &encode_json($data_store_struct);

        # print "$json\n";

        $idx->store_key_value($gene_pair, $json);

        print $ofh join("\t", $gene_A, &encode_json(\@A_ranges), $gene_B, &encode_json(\@B_ranges)) . "\n";
        
    }
      
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
    
