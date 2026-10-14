#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use File::Basename;
use File::Temp qw(tempdir);
use Digest::MD5;
use Process_cmd;
use CTAT_supplement_gtf;


my $usage = <<__EOUSAGE__;

##################################################################################################
#
#  Adds supplemental transcript models to an annotation GTF, for building a CTAT genome lib that
#  includes them (used by prep_genome_lib.pl --supp_gtf).
#
#  Supplemental transcripts are assigned to genes of the input GTF exactly as by
#  supplement_ctat_genome_lib.pl (which adds them to an already built lib), and rewritten to
#  carry that gene's gene_id, gene_name and gene_type.
#
#  Required:
#
#   --in_gtf <string>             annotation GTF
#
#   --supp_gtf <string>           GTF of supplemental transcripts (can be given more than once)
#
#   --out_gtf <string>            output: the input GTF plus the normalized supplemental transcripts
#
#  Optional:
#
#   --source_tag <string>         source (column 2) given to the supplemental records (default: SUPPLEMENTAL)
#
#   --assign_by_overlap           assign transcripts lacking a recognizable gene_id/gene_name to the
#                                 single gene whose exons they overlap on the same strand
#
#   --allow_novel_genes           keep transcripts that cannot be assigned to a gene as new genes
#                                 (default: report them and stop)
#
#   --novel_gene_type <string>    gene_type given to novel genes (default: novel_gene)
#
#   --provenance_out <string>     write a table of the supplements added (file, md5, transcripts)
#
#  Unlike supplement_ctat_genome_lib.pl there is no --allow_nonoverlapping: a transcript assigned
#  to a gene whose span it doesn't overlap (eg. coordinates from a different genome build) is an error.
#
##################################################################################################


__EOUSAGE__

    ;


my $help_flag;
my $in_gtf;
my @supp_gtfs;
my $out_gtf;
my $source_tag = "SUPPLEMENTAL";
my $assign_by_overlap_flag = 0;
my $allow_novel_genes_flag = 0;
my $novel_gene_type = "novel_gene";
my $provenance_out;

&GetOptions ( 'h' => \$help_flag,
              'in_gtf=s' => \$in_gtf,
              'supp_gtf=s' => \@supp_gtfs,
              'out_gtf=s' => \$out_gtf,
              'source_tag=s' => \$source_tag,
              'assign_by_overlap' => \$assign_by_overlap_flag,
              'allow_novel_genes' => \$allow_novel_genes_flag,
              'novel_gene_type=s' => \$novel_gene_type,
              'provenance_out=s' => \$provenance_out,
    );

if ($help_flag) {
    die $usage;
}
unless ($in_gtf && @supp_gtfs && $out_gtf) {
    die $usage;
}
if (@ARGV) {
    die "Error, unrecognized parameters: @ARGV\n";
}
if ($source_tag =~ /\s/) {
    die "Error, --source_tag cannot contain whitespace\n";
}

my $UTILDIR = $FindBin::Bin;


main: {

    foreach my $file ($in_gtf, @supp_gtfs) {
        unless (-s $file) {
            die "Error, cannot locate file: $file";
        }
    }
    my %seen;
    foreach my $supp_gtf (@supp_gtfs) {
        if ($seen{$supp_gtf}++) {
            die "Error, supplemental gtf given more than once: $supp_gtf";
        }
    }

    my $workdir = tempdir("merge_supplemental_gtfs.XXXXXX", TMPDIR => 1, CLEANUP => 1);

    ## gene spans and exon records as the lib build derives them
    print STDERR "-deriving gene spans from $in_gtf\n";
    my $gene_spans_file = "$workdir/in.gene_spans";
    &process_cmd("$UTILDIR/gtf_to_gene_spans.pl $in_gtf > $gene_spans_file");

    my %existing_transcript_ids;
    &load_transcript_ids($in_gtf, \%existing_transcript_ids, "reference");

    my @normalized_gtfs;
    my @provenance;

    for (my $i = 0; $i <= $#supp_gtfs; $i++) {
        my $supp_gtf = $supp_gtfs[$i];

        ## genes and transcript ids of supplements added before this one count too
        my $genes_href = &load_gene_spans($gene_spans_file);
        &add_prior_supplement_genes($genes_href, \@normalized_gtfs);
        my $lookup = &build_gene_lookups($genes_href);

        my $ref_exons_href;
        if ($assign_by_overlap_flag) {
            my $exon_records = "$workdir/exon_records.$i";
            &process_cmd("cat $in_gtf @normalized_gtfs | $UTILDIR/gtf_to_exon_gene_records.pl /dev/stdin > $exon_records");
            $ref_exons_href = &load_ref_exons([$exon_records], $genes_href);
        }

        print STDERR "-parsing supplemental gtf: $supp_gtf\n";
        my ($transcripts_href, $parse_report) = &parse_supplemental_gtf($supp_gtf);

        my ($normalized_lines_aref, $report) = &assign_and_normalize($transcripts_href, $genes_href, $lookup,
                                                                     \%existing_transcript_ids, $ref_exons_href,
                                                                     { source_tag => $source_tag,
                                                                       allow_novel_genes => $allow_novel_genes_flag,
                                                                       novel_gene_type => $novel_gene_type,
                                                                       allow_nonoverlapping => 0 });
        $report->{skipped_unstranded} = $parse_report->{skipped_unstranded};
        $report->{dropped_gene_records} = $parse_report->{dropped_gene_records};
        unshift (@{$report->{errors}}, @{$parse_report->{errors}});

        &print_report($report);

        ## (this utility has no --allow_nonoverlapping)
        s/ \(see --allow_nonoverlapping\)// foreach @{$report->{errors}};

        if (my @errors = @{$report->{errors}}) {
            die "\nError, " . scalar(@errors) . " problem(s) found in supplemental gtf $supp_gtf:\n  "
                . join("\n  ", @errors) . "\n\n";
        }
        unless (@$normalized_lines_aref) {
            die "Error, no supplemental transcripts remain from $supp_gtf after normalization\n";
        }

        my $normalized_gtf = "$workdir/supp.$i.gtf";
        open(my $ofh, ">$normalized_gtf") or die "Error, cannot write $normalized_gtf";
        print $ofh join("\n", @$normalized_lines_aref) . "\n";
        close $ofh or die "Error, writing $normalized_gtf failed";

        &load_transcript_ids($normalized_gtf, \%existing_transcript_ids, "earlier supplement (" . basename($supp_gtf) . ")");
        push (@normalized_gtfs, $normalized_gtf);

        push (@provenance, join("\t", basename($supp_gtf), &md5_file($supp_gtf), $report->{num_transcripts_kept}, $source_tag));
    }

    ## input gtf, then the supplements, each starting on a new line
    my $tmp_out = "$out_gtf.tmp.$$";
    open(my $ofh, ">$tmp_out") or die "Error, cannot write $tmp_out";
    foreach my $file ($in_gtf, @normalized_gtfs) {
        open(my $fh, $file) or die "Error, cannot open $file";
        my $last = "";
        while (my $line = <$fh>) {
            print $ofh $line;
            $last = $line;
        }
        close $fh;
        if ($last ne "" && $last !~ /\n$/) {
            print $ofh "\n";
        }
    }
    close $ofh or die "Error, writing $tmp_out failed";
    rename($tmp_out, $out_gtf) or die "Error, cannot rename $tmp_out to $out_gtf";

    if ($provenance_out) {
        open(my $pfh, ">$provenance_out") or die "Error, cannot write $provenance_out";
        print $pfh join("\t", "#supplemental_gtf", "md5", "num_transcripts", "source_tag") . "\n";
        print $pfh join("\n", @provenance) . "\n";
        close $pfh;
    }

    print STDERR "-wrote $out_gtf (" . scalar(@supp_gtfs) . " supplement(s) added)\n";

    exit(0);
}


####
sub md5_file {
    my ($file) = @_;

    open(my $fh, $file) or die "Error, cannot open $file";
    binmode($fh);
    my $md5 = Digest::MD5->new->addfile($fh)->hexdigest;
    close $fh;

    return($md5);
}
