#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Fasta_reader;
use Pipeliner;

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

my $UTILDIR = "$FindBin::Bin";

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
    #my %chr_to_mask_regions = &get_pseudogene_coordinates($gencode_gtf);  # maybe not such a great idea after all. Must be more targeted - specific ones.
    
    my %chr_to_mask_regions;
    
    ## extract non-pseudogene cdna sequences
    my $tmpdir = "__paralogs";
    unless (-d $tmpdir) {
        mkdir($tmpdir);
    }

    my $pipeliner = new Pipeliner(-verbose => 2,
                                  -checkpoint_dir => $tmpdir);

    my $paramask_file = "para.mask.list";
    my $cmd = "$UTILDIR/define_paralog_mask_targets.pl --gtf $gencode_gtf --genome_fa $genome_fa --tmpdir $tmpdir --output $paramask_file";
    $pipeliner->add_commands(new Command($cmd, "define_para_mask.ok"));
    $pipeliner->run();
    
    &append_mask_regions($paramask_file, \%chr_to_mask_regions);

    ## incorporate the pseudoautosomal regions on chrY
    &append_PAR_regions($gencode_gtf, \%chr_to_mask_regions);
    
    
    ## do making:
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




####
sub get_pseudogene_coordinates {
    my ($gencode_gtf) = @_;
    
    my %chr_to_mask_regions;

    open(my $fh, $gencode_gtf) or die "Error, cannot open file: $gencode_gtf";
    while(<$fh>) {
        my @x = split(/\t/);
        if (scalar(@x) > 8 && $x[2] eq "exon") {
            if ($x[8] =~ /gene_type \"([^\"]+)\"/) {
                my $gene_type = $1;
                
                ## customization for DUX4 pseudogenes labeled as coding genes in hg19 gencode v19
                my $gene_name;
                if ($x[8] =~ /gene_name \"([^\"]+)\"/) {
                    $gene_name = $1;
                }
                
                ## masking out pseudogenes.
                if ($gene_type =~ /pseudogene/) {
                    my ($chr, $lend, $rend) = ($x[0], $x[3], $x[4]);
                    push (@{$chr_to_mask_regions{$chr}}, [$lend, $rend]);
                }
            }
        }
    }
    
    return(%chr_to_mask_regions);
}

####
sub append_mask_regions {
    my ($paramask_file, $chr_to_mask_regions_href) = @_;
    
    open(my $fh, $paramask_file) or die "Error, cannot open file: $paramask_file";
    while(<$fh>) {
        chomp;
        my ($chr, $lend, $rend) = split(/\t/);
        push (@{$chr_to_mask_regions_href->{$chr}}, [$lend, $rend]);
    }
    close $fh;

    return;
}


####
sub append_PAR_regions {
    my ($gencode_gtf, $chr_to_mask_regions_href) = @_;
 

    my $found_PAR = 0;
   
    open(my $fh, $gencode_gtf) or die "Error, cannot open file: $gencode_gtf";
    while(<$fh>) {
        my $line = $_;
        chomp;
        
        my @x = split(/\t/);
        unless (scalar(@x) > 8) { next; }
        
        unless ($x[2] eq "gene") { next; }
        
        if ($x[8] =~ /tag \"PAR\"/) {
            unless ($x[0] eq "chrY") {
                die "Error! found PAR tag feature and not on chrY:  [$line]  ";
            }
            push (@{$chr_to_mask_regions_href->{$x[0]}}, [$x[3], $x[4]]);
            $found_PAR = 1;
        }
    }
    close $fh;

    unless($found_PAR) {
        die "Error, didn't locate PAR features in $gencode_gtf";
    }
    
    return;
}
