#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

####
sub parse_isoform_gtf {
    my ($gtf_file) = @_;

    print STDERR "-parsing GTF file: $gtf_file\n";
    my %isoform_info;
    
    open(my $fh, $gtf_file) or die "Error, cannot open file: $gtf_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        my $line = $_;
        my @x = split(/\t/);
        my $chr = $x[0];
        my $feat_type = $x[2];

        unless($feat_type eq 'exon') { next; }
                
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        
        my $info = $x[8];

        my $gene_name;
        if ($info =~ /gene_name \"([^\"]+)/) {
            $gene_name = $1;
        }
        elsif ($info =~ /gene_id \"([^\"]+)/) {
            $gene_name = $1;
        }
        else {
            print STDERR "-not finding gene_id or gene_name for entry: $line\nskipping...\n";
            next;
        }

        my $transcript_id;
        if ($info =~ /transcript_id \"([^\"]+)/) {
            $transcript_id = $1;
        }
        else {
            print STDERR "-not finding transcript_id for line $line\nskipping...\n";
            next;
        }
        

        
        my $struct = $isoform_info{$transcript_id};
        unless ($struct) {
            $struct = $isoform_info{$transcript_id} = { gene_id => $gene_name,
                                                        transcript_id => $transcript_id,
                                                        chr => $chr,
                                                        orient => $orient,
                                                        exons => [],
            };
        }

        push (@{$struct->{exons}}, { lend => $lend,
                                     rend => $rend,
              } );

    }
    close $fh;
    

    foreach my $struct (values (%isoform_info)) {
        &set_rel_coords($struct);
    }

    return(%isoform_info);
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
    

####
sub translate_to_genomic_coord {
    my ($struct, $coord) = @_;

    my $orient = $struct->{orient};

    my @exons = @{$struct->{exons}};
    foreach my $exon (@exons) {

        if ($exon->{rel_lend} <= $coord && $exon->{rel_rend} >= $coord) {

            my $delta = $coord - $exon->{rel_lend};

            if ($orient eq '+') {
                return($exon->{lend} + $delta);
            }
            else {
                # reverse strand
                return($exon->{rend} - $delta);
            }
        }
    }


    confess "Error, didn't locate coordinate $coord within " . Dumper($struct);
}


1;
