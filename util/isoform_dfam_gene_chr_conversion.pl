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
require "$FindBin::Bin/../lib/isoform_to_genome_coord_conversion.ph";


my $usage = <<__EOUSAGE;

##################################################################
#
# --dfam_results <string>       dfam result file 
#
# --gtf <string>                genome annotation gtf file
#
#################################################################

__EOUSAGE

    ;

my $dfam_results;
my $gtf_file;

my $help_flag;

&GetOptions ( 'h' => \$help_flag,
              
              'dfam_results=s' => \$dfam_results,
              'gtf=s' => \$gtf_file,
    );


unless ($dfam_results && $gtf_file) {
    die $usage;
}




main: {

    my %isoform_info = &parse_isoform_gtf($gtf_file);

    open(my $fh, $dfam_results) or die "Error, cannot open file: $dfam_results";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        s/^\s+//; # trim leading ws.
        
        my @x = split(/\s+/);
        my ($lend, $rend) = sort {$a<=>$b} ($x[11], $x[12]); # using the envelope region
        my $isoform_acc = $x[2];
        
        my $isoform_struct = $isoform_info{$isoform_acc} or die "Error, no isoform struct for [$isoform_acc]";
        my $gene_id = $isoform_struct->{gene_id};
        my $chr = $isoform_struct->{chr};
        
        my $genome_lend = &translate_to_genomic_coord($isoform_struct, $lend);
        my $genome_rend = &translate_to_genomic_coord($isoform_struct, $rend);
        
        ($genome_lend, $genome_rend) = sort {$a<=>$b} ($genome_lend, $genome_rend);
        
        print join("\t", $isoform_acc,$lend, $rend, $gene_id, $chr, $genome_lend, $genome_rend) . "\n";
                        
    }
        
    print STDERR "-done\n";
    
        
    exit(0);
}

