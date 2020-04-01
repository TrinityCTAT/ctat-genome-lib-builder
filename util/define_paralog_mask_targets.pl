#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Pipeliner;
use File::Basename;
use Data::Dumper;

my @MANUALLY_DEFINED_REFERENCE_GENES = qw( DUX4 SEPT14 );

my $CPU = 4;

my $usage = <<__EOUSAGE__;

###########################################################
#
# --gtf <str>          input gtf file
#
# --genome_fa <str>    genome_fa 
#
# --tmpdir <str>       destination for intermediate outputs
#
# --output <str>       output list of targets as redundant.
#
#  optional
# 
# --CPU <int>          default: $CPU
#
###########################################################


__EOUSAGE__

    ;


my $UTILDIR = $FindBin::Bin;

my $help_flag;

my $in_gtf;
my $genome_fa;
my $output_file;

my $tmpdir = $ENV{TMPDIR} || "/tmp";


&GetOptions ( 'h' => \$help_flag,

              'gtf=s' => \$in_gtf,
              'genome_fa=s' => \$genome_fa,
              
              'tmpdir=s' => \$tmpdir,
              'output=s' => \$output_file,
              
              'CPU=i' => \$CPU,
    );

if ($help_flag) {
    die $usage;
}

unless ($in_gtf && $genome_fa && $tmpdir && $output_file) {
    die $usage;
}


my $CONFOUNDERS_FASTA = "$UTILDIR/data/confounders.fa";

main: {
    
    if (! -d $tmpdir) {
        mkdir($tmpdir) or die "Error, cannot mkdir $tmpdir";
    }
    my $chckpts_dir = "$tmpdir/__chckpts";
    
    my $pipeliner = new Pipeliner('-verbose' => 2,
                                  '-checkpoint_dir' => $chckpts_dir);
    
    
    my $blast_outfile = "$tmpdir/confounders.blastn.outfmt6";
    my $cmd = "blastn -db $genome_fa -query $CONFOUNDERS_FASTA -outfmt 6 -evalue 1e-10 > $blast_outfile";
    $pipeliner->add_commands(new Command($cmd, "blast_confounders.ok"));
    $pipeliner->run();
    
    my %chr_to_exclusion_zones = &get_gene_exclusion_zones($in_gtf, \@MANUALLY_DEFINED_REFERENCE_GENES);
    

    open(my $ofh, ">$output_file") or confess "Error, cannot write to file: $output_file";
    open(my $fh, $blast_outfile) or confess "Error, cannot read file: $blast_outfile";
    while(my $line = <$fh>) {
        chomp $line;
        my @x = split(/\t/, $line);
        
        my $chr = $x[1];
        my $query_acc = $x[0];
        my $lend = $x[8];
        my $rend = $x[9];
        
        ($lend, $rend) = sort {$a<=>$b} ($lend, $rend);
        
        if (&in_exclusion_zone($chr, $lend, $rend, \%chr_to_exclusion_zones)) {
            print STDERR "-blast match to $query_acc in exclusion zone: $chr:$lend-$rend\n";
        } 
        else {
            print STDERR "-blast match to $query_acc TARGETED FOR MASKING: $chr:$lend-$rend\n";
            print $ofh join("\t", $chr, $lend, $rend) . "\n";
        }
    }
    close $fh;
    close $ofh;
    
    
    exit(0);
}


####
sub get_gene_exclusion_zones {
    my ($gtf, $exclusion_genes_aref) = @_;
        
    my %exclusion_genes = map { + $_ => 1 } @$exclusion_genes_aref;
        
    my %chr_to_exclusion_zones;
        
    open(my $fh, $gtf) or die "Error, cannot open file: $gtf";
    while(<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        my @x = split(/\t/);
        if (scalar(@x) < 9) { next; }
        
        unless ($x[2] eq "transcript") { next; }

        my $chr = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        
        my $info = $x[8];
        
        my ($gene_name, $gene_id, $transcript_id);

        if ($info =~ /gene_name\s+\"([^\"]+)\"/) {
            $gene_name = $1;
        }
        
        if ($info =~ /gene_id\s+\"([^\"]+)\"/) {
            $gene_id = $1;
        }
        
        if ($info =~ /transcript_id\s+\"([^\"]+)\"/) {
            $transcript_id = $1;
        }
        
                

        unless ($transcript_id && ($gene_id || $gene_name)) {
            print STDERR "-warning, missing transcript_id or (gene_name|gene_id) info in $info\n";
            next;
        }
        
        if ( ($gene_id && $exclusion_genes{$gene_id}) 
             ||
             ($gene_name && $exclusion_genes{$gene_name}) ) {
            
            
            push(@{$chr_to_exclusion_zones{$chr}}, [$lend, $rend]);
            
            my $gene_val_report = ($gene_id && $exclusion_genes{$gene_id}) ? $gene_id : $gene_name;
            
            print STDERR "-adding mask exclusion zone for $gene_val_report $chr:$lend-$rend\n";
        }
    }
    close $fh;
    
    return(%chr_to_exclusion_zones);
    
}


####
sub add_isoform_eclusion_regions {
    my ($gene_name, $chr_to_selected_regions_href, $gene_name_to_trans_href, $trans_to_chr_coordspans_href) = @_;

    my @transcripts = keys %{$gene_name_to_trans_href->{$gene_name}};
    foreach my $transcript (@transcripts) {
        my ($chr, $lend, $rend) = @{$trans_to_chr_coordspans_href->{$transcript}};
        push (@{$chr_to_selected_regions_href->{$chr}}, [$lend, $rend]);
    }
    
    return;
}


####
sub get_gtf_info {
    my ($gtf_file, $gene_name_to_trans_href, $trans_to_chr_coordspans_href) = @_;

    ## TODO: - use reusable code here...
    
    my %trans_id_to_chr;
    my %trans_id_to_coords;
    
    open(my $fh, $gtf_file) or die "Error, cannot open file: $gtf_file";
    while(<$fh>) {
        chomp;
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        my @x = split(/\t/);
        unless (scalar(@x) > 8) { next; }
        unless ($x[2] eq "exon") { next; }
        
        my $chr = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        
        my $info = $x[8];
        my ($transcript_id, $gene_id, $gene_name);
        if ($info =~ /transcript_id \"([^\"]+)\"/) {
            $transcript_id = $1;
        }
        if ($info =~ /gene_id \"([^\"]+)\"/) {
            $gene_id = $1;
        }
        if ($info =~ /gene_name \"([^\"]+)\"/) {
            $gene_name = $1;
        }
        
        if ($transcript_id) {
            if ($gene_id) {
                $gene_name_to_trans_href->{$gene_id}->{$transcript_id} = 1;
            }
            if ($gene_name && $gene_id ne $gene_name) {
                $gene_name_to_trans_href->{$gene_name}->{$transcript_id} = 1;
            }
            
            $trans_id_to_chr{$transcript_id} = $chr;
            push(@{$trans_id_to_coords{$transcript_id}}, $lend, $rend);
        }
    }
    close $fh;

    ## set transcript coordinate spans.
    foreach my $transcript_id (keys %trans_id_to_coords) {
        my @coords = sort {$a<=> $b} @{$trans_id_to_coords{$transcript_id}};
        my $lend = shift @coords;
        my $rend = pop @coords;
        
        my $chr = $trans_id_to_chr{$transcript_id};
        $trans_to_chr_coordspans_href->{$transcript_id} = [$chr, $lend, $rend];
    }

    return;
}
            
####
sub in_exclusion_zone {
    my ($chr, $lend, $rend, $chr_to_exclusion_zones_href) = @_;

    if (exists $chr_to_exclusion_zones_href->{$chr}) {
        
        foreach my $coordset (@{$chr_to_exclusion_zones_href->{$chr}}) {
            my ($zone_lend, $zone_rend) = @$coordset;
            if ($lend < $zone_rend && $rend > $zone_lend) {
                return(1); # yes, in exclusion zone.
            }
        }
    }

    return(0); # no, not in exclusion zone
}
