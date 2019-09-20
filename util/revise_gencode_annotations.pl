#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Pipeliner;
use File::Basename;


my $usage = <<__EOUSAGE__;

Creates IGH and IGL super-loci

######################################################
#
# --gencode_gtf <str>     gencode GTF annotation file
#
# --out_gtf <str>         output masked genome file
#
#######################################################


__EOUSAGE__

    ;


my $UTILDIR = $FindBin::Bin;

my $help_flag;

my $gencode_gtf;
my $out_gtf;


&GetOptions ( 'h' => \$help_flag,

              'gencode_gtf=s' => \$gencode_gtf,
              'out_gtf=s' => \$out_gtf,
    );

if ($help_flag) {
    die $usage;
}

unless ($gencode_gtf && $out_gtf) {
    die $usage;
}


main: {

    my $gtf_file = $gencode_gtf;

    print STDERR "-revising human gencode annotations\n";
    
    my $tmpdir = "__gencode_refinement_chkpts";

    my $pipeliner = new Pipeliner(-verbose => 2,
                                  -checkpoint_dir => $tmpdir);

    ## filter out the IG features that we'll add back later as super-loci
    my $feature_selected_gtf = "$tmpdir/" . basename($gtf_file) . ".feature_selected";
    my $cmd = "bash -c \"set -euxo pipefail; cat $gtf_file |  egrep -v 'IG_V_gene|IG_C_gene|IG_D_gene|IG_J_gene' > $feature_selected_gtf\" ";
    $pipeliner->add_commands(new Command($cmd, "feature_selection.ok"));
        
    ## remove readthru transcripts
    my $no_readthrus_gtf = "$feature_selected_gtf.noreadthrus";
    $cmd = "$UTILDIR/remove_long_intron_readthru_transcripts.pl $feature_selected_gtf 100000 > $no_readthrus_gtf";
    $pipeliner->add_commands(new Command($cmd, "remove_readthrus.ok"));
        
    ########################
    ## create IGH superlocus
    $cmd = "bash -c \"set -euxo pipefail; cat $gtf_file | egrep 'IG_V_gene|IG_C_gene|IG_D_gene|IG_J_gene' | awk '{ if (\\\$3 == \\\"exon\\\") { print } }' | egrep ^chr14  > $tmpdir/IGH_locus.gtf\"";
    $pipeliner->add_commands(new Command($cmd, "igh_locus.ok"));

    $cmd = "$UTILDIR/make_super_locus.pl $tmpdir/IGH_locus.gtf IGH\@ IGH.g\@ IGH.t\@ > $tmpdir/IGH.superlocus.gtf";
    $pipeliner->add_commands(new Command($cmd, "igh_superlocus.ok"));

    $cmd = "bash -c \"set -euxo pipefail; cat $tmpdir/IGH.superlocus.gtf | perl -lane 's/IGH/IGH-/g; print;'  | perl -lane '\@x = split(/\\t/); \\\$x[6] = \\\"+\\\"; print join(\\\"\\t\\\", \@x);' >  $tmpdir/IGH.superlocus.revcomp.gtf\" ";
    $pipeliner->add_commands(new Command($cmd, "igh_superlocus.revcomp.ok") );

    
    ########################
    ## create IGL superlocus
    
    ## extract the exon features for IGL
    $cmd = "bash -c \"set -euxo pipefail; cat $gtf_file | egrep 'IG_V_gene|IG_C_gene|IG_D_gene|IG_J_gene' | awk '{ if (\\\$3 == \\\"exon\\\") { print } }' | egrep ^chr22 > $tmpdir/IGL_locus.gtf\" ";
    $pipeliner->add_commands(new Command($cmd, "igl_locus.ok") );
    
    ## make IGL super-locus
    $cmd = "$UTILDIR/make_super_locus.pl $tmpdir/IGL_locus.gtf IGL\@ IGL.g\@ IGL.t\@ > $tmpdir/IGL.superlocus.gtf";
    $pipeliner->add_commands(new Command($cmd, "igl_superlocus.ok"));

## add reverse complement
    $cmd = "bash -c \"set -euxo pipefail; cat $tmpdir/IGL.superlocus.gtf  | perl -lane 's/IGL/IGL-/g; print;' | perl -lane '\@x = split(/\t/); \\\$x[6] = \\\"-\\\"; print join(\\\"\\t\\\", \@x);' >  $tmpdir/IGL.superlocus.revcomp.gtf \" ";
    $pipeliner->add_commands(new Command($cmd, "igl_superlocus.revcomp.ok") );


    ##########################
    ## Generate final annotation file:

    my $refined_gtf = "$no_readthrus_gtf.refined.gtf";
    $cmd = "cat $no_readthrus_gtf $tmpdir/IGH.superlocus.gtf $tmpdir/IGH.superlocus.revcomp.gtf $tmpdir/IGL.superlocus.gtf $tmpdir/IGL.superlocus.revcomp.gtf > $refined_gtf";
    
    $pipeliner->add_commands(new Command($cmd, "refined_gtf.ok"));
    
    $cmd = "cp $refined_gtf $out_gtf";
    $pipeliner->add_commands(new Command($cmd, "cp_refined_gtf_to_out.ok"));

    $pipeliner->run();
    
    
        
    exit(0);
}


