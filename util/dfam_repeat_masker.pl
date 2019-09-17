#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Pipeliner;
use File::Basename;
use Fasta_reader;

use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $CPU = 4;

my $usage = <<__EOUSAGE__;

########################################################
#
# --dfam_hmm <string>          path to Dfam.hmm library (get it at: http://dfam.org/releases/Dfam_3.1/families/Dfam.hmm.gz)
#
# --target_fa <string>         path to target fasta file
#
# --out_masked <string>        path to masked fasta file
#
# --CPU <int>                  num threads (default: $CPU)
#
#######################################################

__EOUSAGE__

    ;


my $help_flag;
my $dfam_hmm;
my $target_fa;
my $out_masked_fa;

&GetOptions ( 'h' => \$help_flag,
      
              'dfam_hmm=s' => \$dfam_hmm,
              'target_fa=s' => \$target_fa,
              'out_masked=s' => \$out_masked_fa,
              'CPU=i' => \$CPU,
    );

if ($help_flag) { 
    die $usage;
}

unless ($dfam_hmm && $target_fa && $out_masked_fa) {
    die $usage;
}

main: {
    
    ## run dfamscan.pl
    my $chckpts_dir = "__dfam_" . basename($target_fa);
    my $repeat_regions_file = "$chckpts_dir/dfam.out";
    my $cmd = "dfamscan.pl -fastafile $target_fa -hmmfile $dfam_hmm -dfam_outfile $repeat_regions_file --masking_thresh --cpu $CPU";
    
    
    my $pipeliner = new Pipeliner('-verbose' => 2,
                                  '-checkpoint_dir' => $chckpts_dir);

    $pipeliner->add_commands(new Command($cmd, "dfamscan.ok"));
    $pipeliner->run();
    
    my %repeat_regions = &parse_repeat_regions($repeat_regions_file);
    
    
    open(my $ofh, ">$out_masked_fa") or die "Error, cannot write to file: $out_masked_fa";
    my $fasta_reader = new Fasta_reader($target_fa);
    while (my $seq_obj = $fasta_reader->next()) {
        my $acc = $seq_obj->get_accession();
        my $header = $seq_obj->get_header();
        my $sequence = $seq_obj->get_sequence();
        
        if (my $repeat_regions_aref = $repeat_regions{$acc}) {
            my @seqchars = split(//, $sequence);
            foreach my $repeat_region_aref (@$repeat_regions_aref) {
                my ($lend, $rend) = @$repeat_region_aref;
                for (my $i = $lend; $i <= $rend; $i++) {
                    $seqchars[$i-1] = "N";
                }
            }
            $sequence = join("", @seqchars);
        }

        print $ofh ">$header\n$sequence\n";
    }
    
    close $ofh;
    
    exit(0);
    
}
    

####
sub parse_repeat_regions {
    my ($repeat_regions_file) = @_;

    my %repeat_regions;

    open(my $fh, $repeat_regions_file) or die "Error, cannot open file: $repeat_regions_file";
    while(<$fh>) {
        chomp;
        if (/^\#/) { next; }
        s/^\s+//; # trim leading ws.
        
        my @x = split(/\s+/);
        my ($lend, $rend) = sort {$a<=>$b} ($x[11], $x[12]); # using the envelope region
        my $acc = $x[2];

        push (@{$repeat_regions{$acc}}, [$lend, $rend]);

    }

    return(%repeat_regions);
}
        
        
        
        
        
