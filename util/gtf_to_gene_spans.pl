#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\nusage: $0 annots.gtf\n\n";

my $annots_gtf = $ARGV[0] or die $usage;

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

        my $gene_name = "$gene_id";
        if ($info =~ /gene_name \"([^\"]+)\"/) {
            $gene_name = $1;
        }
        
        my $gene_type = "";
        if ($info =~ /gene_type \"([^\"]+)\"/) {
            $gene_type = $1;
        }
        elsif ($info =~ /gene_biotype \"([^\"]+)\"/) {
            $gene_type = $1;
        }
        
        push (@{$data{$gene_id}->{coords}}, $lend, $rend);
        $data{$gene_id}->{chr} = $chr;
        $data{$gene_id}->{orient} = $orient;
        $data{$gene_id}->{gene_name} = $gene_name if $gene_name;
        $data{$gene_id}->{gene_type} = $gene_type if $gene_type;
    }
    close $fh;

    foreach my $gene (keys %data) {
        my $chr = $data{$gene}->{chr};
        my $orient = $data{$gene}->{orient};

        my @coords = @{$data{$gene}->{coords}};

        @coords = sort {$a<=>$b} @coords;

        my $lend = shift @coords;
        my $rend = pop @coords;
        my $gene_name = $data{$gene}->{gene_name} || ".";
        my $gene_type = $data{$gene}->{gene_type} || ".";
        
        print join("\t", $gene, $chr, $lend, $rend, $orient, $gene_name, $gene_type) . "\n";
    }
    
    exit(0);
}

