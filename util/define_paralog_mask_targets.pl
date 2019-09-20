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

my @MANUALLY_DEFINED_REFERENCE_GENES = qw( DUX4 );


my $usage = <<__EOUSAGE__;

###########################################################
#
# --gtf <str>          input gtf file
#
# --cdna_fasta <str>   output gtf file
#
# --tmpdir <str>       destination for intermediate outputs
#
# --output <str>       output list of targets as redundant.
#
###########################################################


__EOUSAGE__

    ;


my $UTILDIR = $FindBin::Bin;

my $help_flag;

my $in_gtf;
my $cdna_fasta;
my $output_file;

my $tmpdir = $ENV{TMPDIR} || "/tmp";


&GetOptions ( 'h' => \$help_flag,

              'gtf=s' => \$in_gtf,
              'cdna_fasta=s' => \$cdna_fasta,
              
              'tmpdir=s' => \$tmpdir,
              'output=s' => \$output_file,
    );

if ($help_flag) {
    die $usage;
}

unless ($in_gtf && $cdna_fasta && $tmpdir && $output_file) {
    die $usage;
}


main: {
    
    if (! -d $tmpdir) {
        mkdir($tmpdir) or die "Error, cannot mkdir $tmpdir";
    }
    my $chckpts_dir = "$tmpdir/__chckpts";
    
    my $pipeliner = new Pipeliner('-verbose' => 2,
                                  '-checkpoint_dir' => $chckpts_dir);
    
    # run cdhit
    my $cmd = "cd-hit-est -i $cdna_fasta -c 0.99 -aS 0.99 -s 0.9 -M 40000 -T 10 -r 0 -d 0 -G 0 -n 11 -o $tmpdir/cdhit.results";
    $pipeliner->add_commands(new Command($cmd, "cdhit.ok"));

    $pipeliner->run();
    
    my $clstr_file = "$tmpdir/cdhit.results.clstr";
    
    my @gene_clusters = &get_gene_clusters($clstr_file, $in_gtf);
    
    
    my @genes_isoforms_to_mask = &define_genes_n_isoforms_to_mask(\@gene_clusters, \@MANUALLY_DEFINED_REFERENCE_GENES);
    
    
    
    exit(0);
}


####
sub parse_gene_name_info {
    my ($gtf) = @_;
        
    my %trans_id_to_gene_name;
    
    open(my $fh, $gtf) or die "Error, cannot open file: $gtf";
    while(<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        my @x = split(/\t/);
        if (scalar(@x) < 9) { next; }
        
        unless ($x[2] eq "transcript") { next; }

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
        
        my $gene_name_use = ($gene_name) ? $gene_name : $gene_id;

        unless ($transcript_id && $gene_name_use) {
            print STDERR "-warning, missing transcript_id or (gene_name|gene_id) info in $info\n";
            next;
        }
        
        $trans_id_to_gene_name{$transcript_id} = $gene_name_use;
    }
    close $fh;
    
    return(%trans_id_to_gene_name);

}

####
sub get_gene_clusters {
    my ($clstr_file, $in_gtf) = @_;
    
    my %trans_id_to_gene_name = &parse_gene_name_info($in_gtf);
    
    my $gene_clusters_file = "$clstr_file.w_genes";
    open(my $ofh, ">$gene_clusters_file") or die "Error, cannot write to $gene_clusters_file";
    

    my @clusters;
    
    my $curr_cluster_aref = [];
    my %genes_seen;
    my @text_rows;
    my $ref_gene = "";

    open(my $fh, $clstr_file) or die "Error, cannot open file: $clstr_file";
    while(my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>/) {
            # process prev cluster.
            my $num_genes_seen = scalar(keys(%genes_seen));
            if ($num_genes_seen > 1) {
                ## keeper.
                push (@clusters, { num_genes => $num_genes_seen,
                                   text_blurb => join("\n", @text_rows) . "\n",
                                   trans_n_genes => $curr_cluster_aref,
                                   ref_gene => $ref_gene,
                      } );
            }
            
            # start a new cluster.
            $curr_cluster_aref = [];
            %genes_seen = ();
            @text_rows = ($line);
            $ref_gene = "";

            print $ofh "$line\n";
        }
        else {
            chomp $line;
            $line =~ s/\s+$//;
            my ($index, $len, $trans_id, $rest) = split(/\s+/, $line, 4);
            $trans_id =~ s/^>//;
            $trans_id =~ s/\.+$//;
            my $gene_name = $trans_id_to_gene_name{$trans_id} or confess "Error, no gene name found for $trans_id";
            
            $len =~ s/,$//;
            
            my $outline = join("\t", $index, $len, $trans_id, $gene_name, $rest);
            print $ofh "$outline\n";
            push (@$curr_cluster_aref, [$trans_id, $gene_name]);
            push (@text_rows, $outline);
            $genes_seen{$gene_name} = 1;
            if ($rest eq "*") {
                $ref_gene = $gene_name;
            }
        }
    }

    # get last one.
    my $num_genes_seen = scalar(keys(%genes_seen));
    if ($num_genes_seen > 1) {
        ## keeper.
        push (@clusters, { num_genes => $num_genes_seen,
                           text_blurb => join("\n", @text_rows) . "\n",
                           trans_n_genes => $curr_cluster_aref,
                           ref_gene => $ref_gene,
              } );
    }
    
    
    close $fh;
    close $ofh;
    

    # write sorted cluster output file.
    @clusters = reverse sort {$a->{num_genes} <=> $b->{num_genes}
                              ||
                                  length($a->{text_blurb}) <=> length($b->{text_blurb}) } @clusters;
    
    
    {
        my $ordered_multi_gene_clusters_file = "$clstr_file.w_genes.ordered_multi_gene_clusters";
        open(my $ofh, ">$ordered_multi_gene_clusters_file") or die "Error, cannot write to $ordered_multi_gene_clusters_file";    
        foreach my $cluster (@clusters) {
            print $ofh $cluster->{text_blurb};
        }
        close $ofh;
    }
    

    return (@clusters);
}

