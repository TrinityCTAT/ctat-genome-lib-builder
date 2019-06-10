#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use DB_File;

my $usage = "\n\n\tusage: $0 fusion_annot_lib.idx\n\n";

my $fusion_lib = $ARGV[0] or die $usage;

main: {

    my %tied_idx;

    tie (%tied_idx, 'DB_File', $fusion_lib, O_RDWR, 0, $DB_BTREE);
    
    my @tokens = ("BCR--ABL1$;SIMPLE", "geneABC--geneXYZ");

    foreach my $token (@tokens) {
            
        if (my $annot = $tied_idx{$token}) { 
            print "Found $token annot already exists as: $annot\n";
        }
        else {
            print "Storing placeholder: $token.\n";
            $tied_idx{$token} = "__placeholder__";
        }
    }
    untie(%tied_idx);
    
    exit(0);
}


