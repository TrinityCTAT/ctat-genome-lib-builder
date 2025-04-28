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
# --is_grch38             Is GRCh38
#
#######################################################


__EOUSAGE__

    ;


my $help_flag;

my $gencode_gtf;
my $genome_fa;
my $out_masked;
my $is_grch38 = 0;

my $UTILDIR = "$FindBin::Bin";

&GetOptions ( 'h' => \$help_flag,

              'gencode_gtf=s' => \$gencode_gtf,
              'genome_fa=s' => \$genome_fa,
              'out_masked=s' => \$out_masked,
              'is_grch38=s' => \$is_grch38,

    );

if ($help_flag) {
    die $usage;
}

unless ($gencode_gtf && $genome_fa && $out_masked) {
    die $usage;
}


my $IS_GRCH38 = ($genome_fa =~ /GRCh38/ || $is_grch38  ) ? 1 : 0;

main: {

    ## get list of regions to mask out:
    
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

        my $seqlen = length($sequence);

        my $regions_aref = $chr_to_mask_regions{$chr};
        if (ref $regions_aref) {
            my @regions = @$regions_aref;
            print STDERR "-masking " . scalar(@regions) . " pseudogene exonic regions for chr: $chr\n";
            my @seqchars = split(//, $sequence);
            my $counter = 0;
            foreach my $region (@regions) {
                my ($lend, $rend) = sort {$a<=>$b} @$region;

                # ensure in range:
                if ($lend < 1) { $lend = 1; }
                if ($rend > $seqlen) { $rend = $seqlen; }
                
                $counter++;
                print STDERR "\t-[$counter] masking $chr $lend-$rend\n";
                if ($rend - $lend > 1e6 && $chr ne "chrY") {
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

        my $contig = $x[0];
        
        if ($x[8] =~ /tag \"PAR\"/) {
            unless ($x[0] eq "chrY") {
                die "Error! found PAR tag feature and not on chrY:  [$line]  ";
            }


            my $lend = $x[3];
            my $rend = $x[4];
            
            # include 50 kb upstream/downstream regions
            $lend -= 50000;
            $rend += 50000;
            
            if ($lend < 1) { $lend = 1; }
            
            push (@{$chr_to_mask_regions_href->{$contig}}, [$lend, $rend]);
            $found_PAR = 1;
        }
    }
    close $fh;

    unless($found_PAR) {
        if ($IS_GRCH38) {
            # just add manually. 
            # see: https://useast.ensembl.org/info/genome/genebuild/human_PARS.html

            # chromosome:GRCh38:Y:10001 - 2781479 is shared with X: 10001 - 2781479 (PAR1)
            push (@{$chr_to_mask_regions_href->{"chrY"}}, [10001, 2781479]);
            
            # chromosome:GRCh38:Y:56887903 - 57217415 is shared with X: 155701383 - 156030895 (PAR2)
            push (@{$chr_to_mask_regions_href->{"chrY"}}, [56887903, 57217415]);
            
        }
        else {
            die "Error, didn't locate PAR features in $gencode_gtf";
        }
    }
    
    return;
}
