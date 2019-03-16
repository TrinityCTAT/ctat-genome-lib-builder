#!/usr/bin/env perl

use strict;
use warnings;
use Set::IntervalTree;
use Data::Dumper;

my $usage = "\nusage: $0 annots.gtf min_long_intron_length\n\n";

my $annots_gtf = $ARGV[0] or die $usage;
my $min_long_intron_length = $ARGV[1] or die $usage;

main: {

    my %data;

    open (my $fh, $annots_gtf) or die "Error, cannot open file $annots_gtf";
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        unless (/\w/) { next; }
        my @x = split(/\t/);
        my $chr = $x[0];
        my $feat_type = $x[2];

        unless ($feat_type eq 'exon') { next; }
        
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        
        my $info = $x[8];
        
        $info =~ /gene_id \"([^\"]+)\"/ or die "Error, cannot extract gene_id from $info";
        my $gene_id = $1;

        $info =~ /transcript_id \"([^\"]+)\"/ or die "Error, cannot extract transcript_id from $info";
        my $transcript_id = $1;
        
        push (@{$data{$transcript_id}->{coords}}, [$lend, $rend]);
        $data{$transcript_id}->{chr} = $chr;
        $data{$transcript_id}->{orient} = $orient;
        $data{$transcript_id}->{gene_id} = $gene_id;
        $data{$transcript_id}->{info} .= "$info\n";
        
    }
    close $fh;

    ## find the long introns.
    my %transcript_to_long_introns;

    my %chr_to_exon_interval_tree;
    
    foreach my $transcript_id (keys %data) {

        my $transcript_struct = $data{$transcript_id};

        my $chr = $transcript_struct->{chr};
        my $coords_aref = $transcript_struct->{coords};
        my @coords = sort {$a->[0]<=>$b->[0]} @$coords_aref;

        my @long_introns;
        for (my $i = 1; $i <= $#coords; $i++) {
            my $prev_exon_rend = $coords[$i-1]->[1];
            my $curr_exon_lend = $coords[$i]->[0];
            my $intron_lend = $prev_exon_rend + 1;
            my $intron_rend = $curr_exon_lend - 1;

            my $intron_length = $intron_rend - $intron_lend + 1;
            if ($intron_length > $min_long_intron_length) {
                push (@long_introns, [$intron_lend, $intron_rend]);
            }
        }

        if (@long_introns) {
            $transcript_to_long_introns{$transcript_id} = \@long_introns;
        }
        else {
            ## add exons to interval tree
            my $itree = $chr_to_exon_interval_tree{$chr};
            unless ($itree) {
                $itree = $chr_to_exon_interval_tree{$chr} =  Set::IntervalTree->new;
            }
            foreach my $coordset (@coords) {
                my ($lend, $rend) = @$coordset;
                if (abs($rend-$lend) > 0) {
                    $itree->insert($transcript_id, $lend, $rend);
                }
            }
        }
    }


    my %transcripts_to_remove;
    
    foreach my $transcript_w_long_intron (keys %transcript_to_long_introns) {
        my $transcript_struct = $data{$transcript_w_long_intron};

        my $chr = $transcript_struct->{chr};
        my $orient = $transcript_struct->{orient};
        
        my $itree = $chr_to_exon_interval_tree{$chr};

        
        
        my $long_introns_aref = $transcript_to_long_introns{$transcript_w_long_intron};
        foreach my $long_intron (@$long_introns_aref) {
            my ($lend, $rend) = @$long_intron;

            ## search left side
            my $overlaps_left_aref = $itree->fetch($lend-2, $lend-1);

            ## search right side
            my $overlaps_right_aref = $itree->fetch($rend+1, $rend+2);
            
            #print "Long intron: $transcript_w_long_intron [$lend-$rend] overlaps left: " . Dumper($overlaps_left_aref) . " and right: " . Dumper($overlaps_right_aref) . "\n";

          gene_pair_search:
            foreach my $left_overlap (@$overlaps_left_aref) {
                my $left_struct = $data{$left_overlap};
                
                foreach my $right_overlap (@$overlaps_right_aref) {
                    my $right_struct = $data{$right_overlap};

                    if ($left_struct->{orient} eq $orient && $left_struct->{orient} eq $right_struct->{orient}
                        &&
                        $left_struct->{gene_id} ne $right_struct->{gene_id}) {

                        ## possible readthru
                        $transcripts_to_remove{$transcript_w_long_intron} = 1;

                        last gene_pair_search;
                    }
                }
            }
            
   
        }
    }
    
    print STDERR "Removing possible readthru transcripts: \n";
    foreach my $transcript_to_remove (keys %transcripts_to_remove) {
        my $struct = $data{$transcript_to_remove};
        print STDERR Dumper($struct);
    }


    {
        # output the records that passed.

        open (my $fh, $annots_gtf) or die "Error, cannot open file $annots_gtf";
        while (<$fh>) {
            if (/transcript_id \"([^\"]+)\"/) {
                my $transcript_id = $1;
                if ($transcripts_to_remove{$transcript_id}) {
                    next;
                }
            }
            print;
        }
    }

    print STDERR "Done. Removed: " . scalar(keys %transcripts_to_remove) . " likely readthru transcripts.\n";
    
    exit(0);
}

