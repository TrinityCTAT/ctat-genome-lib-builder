#!/usr/bin/env perl

use strict;
use warnings;

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use TiedHash;


my $usage = <<__EOUSAGE__;

#####################################################################
#
# --gene_spans <string>     gene spans file
#
# --out_db_file <string>    output filename for database
#
# optional:
#
# --key_pairs <string>      annotations in key(tab)value format
#
#                           format:
#                           key(tab)simple_annot(tab)complex_annot
#
#                                      where complex_annot is optional, and if included
#                                      is specified as structured json notation.
#
#####################################################################

__EOUSAGE__

    ;


my $help_flag;

my $gene_spans_file;
my $out_db_file;
my $key_pairs_file;

&GetOptions ( 'h' => \$help_flag,

              'gene_spans=s' => \$gene_spans_file,
              'out_db_file=s' => \$out_db_file,
              'key_pairs=s' => \$key_pairs_file,
    );

if ($help_flag) { die $usage; }

unless ($gene_spans_file && $out_db_file) { die $usage; }
              

main: {


    my $idx = new TiedHash( { create => $out_db_file } );

    ## store placeholder to ensure it's working in applications that leverage it.
    $idx->store_key_value("geneABC--geneXYZ", "__placeholder_testval__");
    
    open(my $fh, $gene_spans_file) or die "Error, cannot open file: $gene_spans_file";
    while (<$fh>) {
        chomp;
        my ($gene_id, $chr, $lend, $rend, $orient, $gene_name, $gene_type) = split(/\t/);
        if ($gene_name && $gene_name ne ".") {
            $gene_id = $gene_name;
        }

        $idx->store_key_value("$gene_id$;COORDS", "$chr:$lend-$rend:$orient"); # hacky way of specifying coordinate info for direct coordinate info lookups.
        if ($gene_type) {
            unless($idx->get_value($gene_id)) {
                # store at least the gene type info for the annotation string.
                $idx->store_key_value("$gene_id$;GENE_TYPE", $gene_type);
            }
        }
    }
    close $fh;
    

    ## load in generic annotations:

    if ($key_pairs_file) {

        my $num_records;
        
        if ($key_pairs_file =~ /\.gz$/) {
            $num_records = `gunzip -c $key_pairs_file | wc -l `;
            open($fh, "gunzip -c $key_pairs_file | ") or die "Error, cannot open( gunzip -c $key_pairs_file  )";
            
        }
        else {
            $num_records = `cat $key_pairs_file | wc -l`;
            open($fh, $key_pairs_file) or die "Error, cannot read file $key_pairs_file";
        }

        chomp $num_records;
        my $counter = 0;
        while (<$fh>) {
            $counter++;
            if ($counter % 100000 == 0) {
                my $pct_done = sprintf("%.2f", $counter/$num_records * 100);
                print STDERR "\r[$pct_done\%  done]    ";
            }
            chomp;
            my ($gene_pair, $simple_annot, $complex_annot) = split(/\t/);
            
            if ($simple_annot ne '.') {
                $idx->store_key_value("$gene_pair$;SIMPLE", $simple_annot);
            }
            if ($complex_annot ne '.') {
                $idx->store_key_value("$gene_pair$;COMPLEX", $complex_annot);
            }
        }
        close $fh;
    }
    
    print STDERR "\n\nDone building annot db: $out_db_file\n";
    
    exit(0);
}
           
