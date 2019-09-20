#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Fasta_reader;


my $usage = <<__EOUSAGE__;

All pseudogene regions get masked from the genome.

######################################################
#
# --gencode_gtf <str>     gencode GTF annotation file
#
# --genome_fa <str>       reference genome sequence
#
# --out_masked <str>      output masked genome file
#
#
#######################################################


__EOUSAGE__

    ;


my $help_flag;

my $gencode_gtf;
my $genome_fa;
my $out_masked;


&GetOptions ( 'h' => \$help_flag,

              'gencode_gtf=s' => \$gencode_gtf,
              'genome_fa=s' => \$genome_fa,
              'out_masked=s' => \$out_masked,
    );

if ($help_flag) {
    die $usage;
}

unless ($gencode_gtf && $genome_fa && $out_masked) {
    die $usage;
}


main: {

## get list of regions to mask out:
    my %chr_to_mask_regions;
    {
        open(my $fh, $gencode_gtf) or die "Error, cannot open file: $gencode_gtf";
        while(<$fh>) {
            my @x = split(/\t/);
            if (scalar(@x) > 8 && $x[2] eq "exon") {
                if ($x[8] =~ /gene_type \"([^\"]+)\"/) {
                    my $gene_type = $1;
                    
                    ## masking out pseudogenes.
                    if ($gene_type =~ /pseudogene/) {
                        my ($chr, $lend, $rend) = ($x[0], $x[3], $x[4]);
                        push (@{$chr_to_mask_regions{$chr}}, [$lend, $rend]);
                    }
                }
            }
        }
    }

    ## do masking:
    my $fasta_reader = new Fasta_reader($genome_fa);
    
    open(my $ofh, ">$out_masked") or die "Error, cannot write to $out_masked";
    
    while (my $seq_obj = $fasta_reader->next()) {
        
        my $chr = $seq_obj->get_accession();
        my $header = $seq_obj->get_header();
        my $sequence = $seq_obj->get_sequence();

        my $regions_aref = $chr_to_mask_regions{$chr};
        if (ref $regions_aref) {
            my @regions = @$regions_aref;
            print STDERR "-masking " . scalar(@regions) . " pseudogene exonic regions for chr: $chr\n";
            my @seqchars = split(//, $sequence);
            my $counter = 0;
            foreach my $region (@regions) {
                my ($lend, $rend) = sort {$a<=>$b} @$region;
                $counter++;
                print STDERR "\t-[$counter] masking $chr $lend-$rend\n";
                if ($rend - $lend > 1e6) {
                    die "Error, this pseudogene looks too long!";
                }
                
                for (my $i = $lend; $i <= $rend; $i++) {
                    $seqchars[$i-1] = "N";
                }
            }
            $sequence = join("", @seqchars);
        }
        
        print STDERR "-writing masked sequence for chr: $chr.\n";
        $sequence =~ s/(\S{60})/$1\n/g;
        chomp $sequence;
        
        print $ofh ">$header\n$sequence\n";
    }
    
    close $ofh;

    print STDERR "-done refining genome sequence => $out_masked\n";
    
    
    exit(0);
}


