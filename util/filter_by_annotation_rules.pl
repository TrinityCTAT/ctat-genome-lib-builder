#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use DelimParser;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my $help_flag;


my $usage = <<__EOUSAGE__;

####################################################################################
#
#  --fusions <string>    fusion predictions (must be annotated via FusionAnnotator)
#
#  --genome_lib_dir <string>   CTAT genome lib dir
#
# Optional:
#
#  --min_FFPM <float>     min FFPM 
#
#
####################################################################################

__EOUSAGE__

    ;


my $predictions_file;
my $genome_lib_dir;
my $min_FFPM;

&GetOptions ( 'h' => \$help_flag,

              'fusions=s' => \$predictions_file,
              'genome_lib_dir=s' => \$genome_lib_dir,
   
              'min_FFPM=f' => \$min_FFPM,
    );


unless ($genome_lib_dir) {
    if (exists $ENV{CTAT_GENOME_LIB}) {
        $genome_lib_dir = $ENV{CTAT_GENOME_LIB};
    }
}

unless ($predictions_file && $genome_lib_dir) {
    die $usage;
}


main: {
    
    my $annot_filt_module = "$genome_lib_dir/AnnotFilterRule.pm";
    unless (-s $annot_filt_module) {
        die "Error, cannot locate required $annot_filt_module  ... be sure to use a more modern version of the companion CTAT_GENOME_LIB ";
    }

    require $annot_filt_module;

    open(my $fh, $predictions_file) or die "Error, cannot open file: $predictions_file ";
    
    my $delim_parser = new DelimParser::Reader($fh, "\t");

    my $pass_file = "$predictions_file.pass";
    open(my $ofh_pass, ">$pass_file") or die "Error, cannot write to file: $pass_file";

    
    my $fail_file = "$predictions_file.annot_filt";
    open(my $ofh_fail, ">$fail_file") or die "Error, cannot write to file: $fail_file";

    my @column_headers = $delim_parser->get_column_headers();

    my $pass_writer = new DelimParser::Writer($ofh_pass, "\t", \@column_headers);
    
    my $fail_reason_header = 'Annot_Fail_Reason';
    my $fail_writer = new DelimParser::Writer($ofh_fail, "\t", [@column_headers, $fail_reason_header]);


    while(my $pred = $delim_parser->get_row()) {

        if (my $filter_reason = &AnnotFilterRule::examine_fusion_prediction($pred)) {
            $pred->{$fail_reason_header} = $filter_reason;
            $fail_writer->write_row($pred);
            next;
        }

        if (defined($min_FFPM)) {
            my $ffpm = $pred->{FFPM};
            unless (defined $ffpm) {
                confess "Error, no FFPM value defined in file: $predictions_file, row: " . Dumper($pred);
            }
            if ($ffpm < $min_FFPM) {
                $pred->{$fail_reason_header} = " FFPM $ffpm < min_FFPM required of $min_FFPM ";
                $fail_writer->write_row($pred);
                next;
            }
        }
        
        
        # all good
        $pass_writer->write_row($pred);
        
    }


    close $ofh_pass;
    close $ofh_fail;

    print STDERR "-done, see $pass_file\n";
    
    exit(0);
}


