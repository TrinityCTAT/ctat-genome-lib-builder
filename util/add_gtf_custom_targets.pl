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

my $usage = <<__EOUSAGE__;

Add custom targets to better capture fusions relevant to cancer but
that do not have breakpoints that involve reference isoforms.

Modificiations include:
-addition of an unprocessed CRLF2 transcript that extends upstream by 50kb


######################################################
#
# --in_gtf <str>     input gtf file
#
# --out_gtf <str>    output gtf file
#
# optional:
#
# --tmpdir <str>      destination for intermediate outputs
#
#######################################################


__EOUSAGE__

    ;


my $UTILDIR = $FindBin::Bin;

my $help_flag;

my $in_gtf;
my $out_gtf;

my $tmpdir = $ENV{TMPDIR} || "/tmp";


&GetOptions ( 'h' => \$help_flag,

              'in_gtf=s' => \$in_gtf,
              'out_gtf=s' => \$out_gtf,
              
              'tmpdir=s' => \$tmpdir,

    );

if ($help_flag) {
    die $usage;
}

unless ($in_gtf && $out_gtf) {
    die $usage;
}


main: {

    my %chr_gene_to_span_n_orient = &parse_gene_spans($in_gtf);

    
    my $custom_targets = "$tmpdir/custom_targets.$$.gtf";
    open(my $ofh, ">$custom_targets") or die "Error, cannot write to file: $custom_targets";
    
    { 
        
        ## CRLF2
 
        ## chrX:CRLF2 - add 50kb upstream
        #  crlf2       csf2ra
        # 1212750     1268800   hg19 coords.           
        #       (56050)
        
        print STDERR "-adding custom CRLF2 entry\n";
        my $crlf2_struct = $chr_gene_to_span_n_orient{"chrX:CRLF2"} or confess "Error, cannot find gene struct for chrX:CRLF2 based on $in_gtf";
        my $custom_gtf = &get_custom_RNA_gtf($crlf2_struct, "chrX", "CRLF2", "-", 50000, 0);
        print $ofh "$custom_gtf";
    }
    
    {
        ## MALT1 chr18  - add 40kb upstream

       
        #    ALPK2        MALT1
        #  58629091      58671465
        #        (42374)

        print STDERR "-adding custom MALT1 entry\n";
        my $malt1_struct = $chr_gene_to_span_n_orient{"chr18:MALT1"} or confess "Error, cannot find gene struct for chr18:MALT1 based on $in_gtf";
        my $custom_gtf = &get_custom_RNA_gtf($malt1_struct, "chr18", "MALT1", "+", 40000, 0);
        print $ofh "$custom_gtf";
        
    }

    
    {
        ## DUX4 chr4 - add 10kb upstream
        
        print STDERR "-adding custom DUX4 entry\n";
        my $dux4_struct = $chr_gene_to_span_n_orient{"chr4:DUX4"} or confess "Error, cannot find gene struct for chr4:DUX4 based on $in_gtf";
        my $custom_gtf = &get_custom_RNA_gtf($dux4_struct, "chr4", "DUX4", "+", 10000, 0);
        print $ofh $custom_gtf;
    }
    
    close $ofh;

    my $cmd = "cat $in_gtf $custom_targets > $out_gtf";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }
    print STDERR "-customized gtf file written to: $out_gtf\n";
    
    
    exit(0);
}


####
sub parse_gene_spans {
    my ($gtf) = @_;

    my %chr_gene_name_to_coords;
    my %chr_gene_name_to_orient;
    my %chr_gene_name_to_gene_id;

    print STDERR "-add_gtf_custom_targets: parsing gtf file: $gtf\n";

    open(my $fh, $gtf) or die "Error, cannot open file: $gtf";
    while(<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        my @x = split(/\t/);
        if (scalar(@x) < 9) { next; }

        my $chr = $x[0];
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        my $info = $x[8];
        

        if ($info =~ /gene_name\s+\"([^\"]+)\"/) {
            my $gene_name = $1;
            
            my $chr_gene_name = "$chr:$gene_name";
            $chr_gene_name_to_orient{$chr_gene_name} = $orient;
            push (@{$chr_gene_name_to_coords{$chr_gene_name}}, $lend, $rend);
            
            if ($info =~ /gene_id\s+\"([^\"]+)\"/) {
                my $gene_id = $1;
                $chr_gene_name_to_gene_id{$chr_gene_name} = $gene_id;
            }
        }
    }
    close $fh;
    
    print STDERR "-add_gtf_custom_targets: getting gene span info.\n";

    my %chr_gene_name_to_span_n_orient;

    foreach my $chr_gene_name (keys %chr_gene_name_to_coords) {
        my @coords = sort {$a<=>$b} @{$chr_gene_name_to_coords{$chr_gene_name}};

        my $lend = shift @coords;
        my $rend = pop @coords;

        $chr_gene_name_to_span_n_orient{$chr_gene_name} = { 
            chr_gene_name => $chr_gene_name,
            gene_id => $chr_gene_name_to_gene_id{$chr_gene_name},
            lend => $lend,
            rend => $rend,
            orient => $chr_gene_name_to_orient{$chr_gene_name},
        };
        
    }

    return(%chr_gene_name_to_span_n_orient);

}


####
sub get_custom_RNA_gtf {
    my ($gene_bounds_struct, $chr, $gene_name, $expected_orient, $upstream_extension, $downstream_extension) = @_;
    
    
    my ($gene_id, $lend, $rend, $orient) = ($gene_bounds_struct->{gene_id},
                                            $gene_bounds_struct->{lend},
                                            $gene_bounds_struct->{rend},
                                            $gene_bounds_struct->{orient});
    
    if ($orient ne $expected_orient) {
        confess "Error, $gene_name should be on ($expected_orient) strand! " . Dumper($gene_bounds_struct);
    }

    if ($upstream_extension) {
        if ($orient eq '+') {
            $lend -= $upstream_extension;
        }
        else {
            $rend += $upstream_extension;
        }
    }
    if ($downstream_extension) {
        if ($orient eq "+") {
            $rend += $downstream_extension;
        }
        else {
            $lend -= $downstream_extension;
        }
    }
        
    my $custom_gtf = join("\t", 
                          "$chr",
                          "CUSTOM",
                          "transcript",
                          $lend,
                          $rend,
                          ".",
                          $orient,
                          ".",
                          "gene_id \"$gene_id\"; transcript_id \"$gene_id.trans.custom\"; gene_type \"protein_coding\"; gene_name \"$gene_name\"; transcript_type \"misc_RNA\";\n");

    
    $custom_gtf .= join("\t", 
                        "$chr",
                        "CUSTOM",
                        "exon",
                        $lend,
                        $rend,
                        ".",
                        $orient,
                        ".",
                        "gene_id \"$gene_id\"; transcript_id \"$gene_id.trans.custom\"; gene_type \"protein_coding\"; gene_name \"$gene_name\"; transcript_type \"misc_RNA\";\n");
    
    
    return($custom_gtf);

}
