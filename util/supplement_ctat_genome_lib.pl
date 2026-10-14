#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use FindBin;
use lib ("$FindBin::Bin/../lib");
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use File::Basename;
use File::Copy;
use File::Path;
use File::Temp qw(tempdir);
use File::Spec;
use Cwd qw(abs_path);
use Sys::Hostname;
use Digest::MD5;
use JSON::PP;
use DB_File;
use Fcntl qw(O_WRONLY O_CREAT O_EXCL);
use Process_cmd;


my $usage = <<__EOUSAGE__;

##################################################################################################
#
#  Supplements an existing CTAT genome lib with additional transcript isoforms (eg. isoforms
#  missing from the reference annotation), modifying the lib in place.
#
#  Before any lib file is modified for the first time, the original is preserved as <file>.bkup
#  (a hard link to the original where possible, so the backup takes no extra disk space).
#  Each run regenerates the modified files from the .bkup originals plus ALL supplements added
#  so far (stored under <genome_lib_dir>/supplemental/), so reruns do not stack.
#
#  Revert to the original lib with:  revert_ctat_genome_lib_supplement.pl --genome_lib_dir <dir>
#
#  Required:
#
#   --genome_lib_dir <string>     path to ctat_genome_lib_build_dir
#
#   --supp_gtf <string>           GTF containing the supplemental transcripts
#
#  Optional:
#
#   --label <string>              short name for this supplement (default: basename of --supp_gtf)
#
#   --assign_by_overlap           assign transcripts lacking a recognizable reference gene_id/gene_name
#                                 to the single reference gene whose exons they overlap on the same strand
#
#   --allow_novel_genes           keep transcripts that cannot be assigned to a reference gene as new genes
#                                 (default: report them and stop)
#
#   --novel_gene_type <string>    gene_type given to novel genes (default: novel_gene)
#
#   --linked_copy <string>        leave --genome_lib_dir untouched: create a new lib at this path whose
#                                 files are hard links to the originals (no extra disk space), and
#                                 supplement that instead. Use for shared or read-only libs, or when
#                                 jobs are running against the lib. The copy is self-contained and can
#                                 be mounted on its own in docker/singularity. It must be on the same
#                                 filesystem as --genome_lib_dir, and you must own the lib files (or
#                                 have write access to them), otherwise symlinks are used instead and
#                                 the copy only works where the original lib is visible at its path.
#
#   --allow_nonoverlapping        allow a transcript to be assigned (by gene_id or gene_name) to a gene
#                                 whose span it does not overlap, extending that gene's span to include it
#                                 (default: report it and stop; this usually indicates a different genome
#                                 build, or a gene_id/gene_name reused for an unrelated locus)
#
#   --dry_run                     validate and normalize the supplement, report, and change nothing
#
#   --dry_run_out <string>        with --dry_run, write the normalized supplement (as it would be added to
#                                 the lib, with gene assignments applied) to this file
#
#  Supplemental transcripts are assigned to reference genes, in order of preference, by:
#     gene_id, ref_gene_id (eg. StringTie), version-stripped gene_id/ref_gene_id (reference genes
#     only, same chr and strand), gene_name/ref_gene_name (same chr and strand), then exon overlap
#     (--assign_by_overlap).
#  Assigned transcripts are rewritten to carry the reference gene's gene_id, gene_name and gene_type.
#
#  Files updated:
#     ref_annot.gtf, ref_annot.gtf.mini.sortu, ref_annot.gtf.gene_spans, fusion_annot_lib.idx,
#     ref_annot.gtf.mm2.splice.bed (if present in the lib), and
#     ref_annot.prot_info.dbm, ref_annot.cds, ref_annot.pep (only if supplements include CDS records)
#
#  Not updated: STAR index, blast_pairs, trans.blast, pfam_domains, cdna/cdsplus sequences.
#
##################################################################################################


__EOUSAGE__

    ;


my $help_flag;
my $genome_lib_dir;
my $supp_gtf;
my $label;
my $assign_by_overlap_flag = 0;
my $allow_novel_genes_flag = 0;
my $novel_gene_type = "novel_gene";
my $dry_run_flag = 0;
my $linked_copy;
my $allow_nonoverlapping_flag = 0;
my $dry_run_out;

&GetOptions ( 'h' => \$help_flag,
              'genome_lib_dir=s' => \$genome_lib_dir,
              'supp_gtf=s' => \$supp_gtf,
              'label=s' => \$label,
              'assign_by_overlap' => \$assign_by_overlap_flag,
              'allow_novel_genes' => \$allow_novel_genes_flag,
              'novel_gene_type=s' => \$novel_gene_type,
              'dry_run' => \$dry_run_flag,
              'linked_copy=s' => \$linked_copy,
              'allow_nonoverlapping' => \$allow_nonoverlapping_flag,
              'dry_run_out=s' => \$dry_run_out,
    );

if ($help_flag) {
    die $usage;
}

unless ($genome_lib_dir && $supp_gtf) {
    die $usage;
}

if (@ARGV) {
    die "Error, unrecognized parameters: @ARGV\n";
}
if ($dry_run_out && ! $dry_run_flag) {
    die "Error, --dry_run_out requires --dry_run\n";
}

## (not Process_cmd::ensure_full_path, which splits paths on commas)
$genome_lib_dir = File::Spec->rel2abs($genome_lib_dir);
$genome_lib_dir =~ s{/+$}{};
$supp_gtf = File::Spec->rel2abs($supp_gtf);

## run END cleanup (lock, temp files, incomplete linked copy) when interrupted
my $terminated_by_signal = 0;
$SIG{INT} = $SIG{TERM} = $SIG{HUP} = sub { $terminated_by_signal = 1; die "\nError, terminated by signal $_[0]\n"; };

my $HARDLINK_RULES = <<__EORULES__;
  Hard-linked files share storage: they are the same file under two names. So:
    - do NOT edit lib files in place (eg. in a text editor, or with tools that rewrite a file
      rather than replacing it): the change would appear under both names.
    - supplement_ctat_genome_lib.pl and revert_ctat_genome_lib_supplement.pl never edit files
      in place, so supplementing or reverting either lib is safe.
    - deleting a file, or a whole lib directory, does not affect the other name/lib.
__EORULES__

my $num_copy_hardlinks = 0;
my $num_copy_symlinks = 0;

my $source_lib_dir;
my $created_copy_dir; # set only once this run has created the linked copy; removed on failure
if ($linked_copy && $dry_run_flag) {
    print STDERR "-NOTE: --dry_run: checking the supplement against $genome_lib_dir; no linked copy is created.\n";
}
elsif ($linked_copy) {
    my $copy_dir = &create_linked_copy($genome_lib_dir, File::Spec->rel2abs($linked_copy));
    $created_copy_dir = $copy_dir;
    $source_lib_dir = $genome_lib_dir;
    $genome_lib_dir = $copy_dir;
}

unless ($label) {
    $label = basename($supp_gtf);
    $label =~ s/\.gtf$//i;
}
$label =~ s/[^\w\.\-]/_/g;

my $UTILDIR = $FindBin::Bin;

my $SUPP_DIR = "$genome_lib_dir/supplemental";
my $MANIFEST = "$SUPP_DIR/manifest.json";
my $LOCK_DIR = "$SUPP_DIR/.lock";

my $SOURCE_TAG = "SUPPLEMENTAL";

## lib files this utility may modify.
my @REQUIRED_MANAGED = ("ref_annot.gtf",
                        "ref_annot.gtf.mini.sortu",
                        "ref_annot.gtf.gene_spans",
                        "fusion_annot_lib.idx");

my @OPTIONAL_MANAGED = ("ref_annot.gtf.mm2.splice.bed");

my @PROT_INFO_MANAGED = ("ref_annot.prot_info.dbm",
                         "ref_annot.cds",
                         "ref_annot.pep");

my $have_lock = 0;
my $WORKDIR;
my $num_hardlinked_backups = 0;


main: {

    unless (-d $genome_lib_dir) {
        die "Error, cannot locate genome lib dir: $genome_lib_dir";
    }
    unless (-s $supp_gtf) {
        die "Error, cannot locate supplemental gtf file: $supp_gtf";
    }
    foreach my $file (@REQUIRED_MANAGED, "ref_genome.fa") {
        unless (-s "$genome_lib_dir/$file" || -s "$genome_lib_dir/$file.bkup") {
            die "Error, genome lib is missing required file: $genome_lib_dir/$file";
        }
    }

    unless ($dry_run_flag) {
        unless (-w $genome_lib_dir) {
            die "Error, genome lib dir is not writable: $genome_lib_dir\n"
                . "  (for a read-only install, use --linked_copy to supplement a copy)\n";
        }
        &acquire_lock();
        ## safe only while holding the lock: anything left by an earlier, interrupted run
        &remove_stale_files();
    }

    my $manifest = &read_manifest();
    if ($source_lib_dir) {
        if ($manifest->{source_lib}) {
            ## a linked copy of a linked copy: keep the chain
            $manifest->{earlier_source_libs} = [ @{$manifest->{earlier_source_libs} || []}, $manifest->{source_lib} ];
        }
        $manifest->{source_lib} = $source_lib_dir;
        $manifest->{source_lib_links} = { hard_links => $num_copy_hardlinks, symlinks => $num_copy_symlinks };
    }

    &check_backup_state($manifest);
    unless ($dry_run_flag) {
        &remove_orphan_supplement_files($manifest);
    }

    ## a supplement whose commit was interrupted is recorded with status 'committing'. Each run rebuilds
    ## the lib files from the originals plus every recorded supplement, so any run completes it;
    ## rerunning the interrupted command (same --label) completes it without adding anything new.
    my $resume_id;
    my @interrupted = grep { ($_->{status} || "complete") ne "complete" } @{$manifest->{supplements}};
    my $resume_flag = 0;
    if (@interrupted) {
        if (my ($entry) = grep { $_->{label} eq $label } @interrupted) {
            $resume_flag = 1;
            $resume_id = $entry->{id};
            print STDERR "-supplement '$label' was interrupted while being committed; completing it.\n";
        }
        else {
            print STDERR "-NOTE: supplement(s) " . join(", ", map { "'$_->{label}'" } @interrupted)
                . " were interrupted while being committed; this run completes them.\n";
        }
    }

    my $supp_id = $resume_id || sprintf("supp_%03d", scalar(@{$manifest->{supplements}}) + 1);
    unless ($resume_flag) {
        foreach my $prev (@{$manifest->{supplements}}) {
            if ($prev->{label} eq $label) {
                die "Error, a supplement labeled '$label' was already added ($prev->{normalized_gtf}). Use a different --label.\n";
            }
        }
    }

    if ($dry_run_flag) {
        $WORKDIR = tempdir("supplement_ctat_genome_lib.XXXXXX", TMPDIR => 1, CLEANUP => 1);
    }
    else {
        $WORKDIR = "$SUPP_DIR/.work.$$";
        mkpath($WORKDIR);
        unless (-d $WORKDIR) {
            die "Error, cannot mkdir $WORKDIR";
        }
    }

    my @prior_supp_gtfs = map { "$SUPP_DIR/" . $_->{normalized_gtf} } @{$manifest->{supplements}};

    my $report;
    my $normalized_gtf;

    unless ($resume_flag) {

        ######################################
        ## load reference gene and transcript info

        print STDERR "-loading reference gene spans\n";
        my $genes_href = &load_gene_spans(&src_file("ref_annot.gtf.gene_spans"));

        ## genes introduced by earlier supplements are valid assignment targets too.
        &add_prior_supplement_genes($genes_href, \@prior_supp_gtfs);

        my $lookup = &build_gene_lookups($genes_href);

        print STDERR "-loading reference transcript ids\n";
        my %existing_transcript_ids;
        &load_transcript_ids(&src_file("ref_annot.gtf"), \%existing_transcript_ids, "reference");
        foreach my $prior_gtf (@prior_supp_gtfs) {
            &load_transcript_ids($prior_gtf, \%existing_transcript_ids, "earlier supplement");
        }

        my $ref_exons_href;
        if ($assign_by_overlap_flag) {
            print STDERR "-loading reference exons for overlap-based assignment\n";
            my @exon_record_files = (&src_file("ref_annot.gtf.mini.sortu"));
            if (@prior_supp_gtfs) {
                ## include exons added by earlier supplements
                my $prior_exons = "$WORKDIR/prior_supplements.mini";
                &process_cmd("cat " . join(" ", map { &shq($_) } @prior_supp_gtfs)
                             . " | " . &shq("$UTILDIR/gtf_to_exon_gene_records.pl") . " /dev/stdin > " . &shq($prior_exons));
                push (@exon_record_files, $prior_exons);
            }
            $ref_exons_href = &load_ref_exons(\@exon_record_files, $genes_href);
        }

        ######################################
        ## parse, validate and normalize the supplement

        print STDERR "-parsing supplemental gtf: $supp_gtf\n";
        my ($transcripts_href, $parse_report) = &parse_supplemental_gtf($supp_gtf);

        my $normalized_lines_aref;
        ($normalized_lines_aref, $report) = &assign_and_normalize($transcripts_href, $genes_href, $lookup,
                                                                   \%existing_transcript_ids, $ref_exons_href);
        $report->{skipped_unstranded} = $parse_report->{skipped_unstranded};
        $report->{dropped_gene_records} = $parse_report->{dropped_gene_records};
        unshift (@{$report->{errors}}, @{$parse_report->{errors}});

        &print_report($report);

        if (@{$report->{errors}}) {
            my $num_errors = scalar(@{$report->{errors}});
            my @show = @{$report->{errors}};
            if (@show > 25) { @show = @show[0..24]; }
            die "\nError, $num_errors problem(s) found in the supplemental gtf; no changes made to the genome lib:\n  "
                . join("\n  ", @show) . "\n" . ($num_errors > 25 ? "  ...\n" : "") . "\n";
        }
        unless (@$normalized_lines_aref) {
            die "Error, no supplemental transcripts remain after normalization; no changes made.\n";
        }

        $normalized_gtf = "$WORKDIR/$supp_id.$label.gtf";
        &write_lines($normalized_gtf, $normalized_lines_aref);
    }

    if ($dry_run_flag) {
        if ($dry_run_out && $normalized_gtf) {
            copy($normalized_gtf, $dry_run_out) or die "Error, cannot write $dry_run_out: $!";
            print STDERR "-dry run: normalized supplement written to $dry_run_out\n";
        }
        print STDERR "\n-dry run: genome lib unchanged.\n\n";
        exit(0);
    }

    ######################################
    ## combine all supplements and build replacement lib files (nothing in the lib changes yet)

    my $combined_gtf = "$WORKDIR/all_supplements.gtf";
    &process_cmd("cat " . join(" ", map { &shq($_) } @prior_supp_gtfs, ($normalized_gtf || ())) . " > " . &shq($combined_gtf));

    my $combined_has_CDS = &gtf_has_feature($combined_gtf, "CDS");

    my @managed = (@REQUIRED_MANAGED);
    foreach my $file (@OPTIONAL_MANAGED) {
        if (-e "$genome_lib_dir/$file" || -e "$genome_lib_dir/$file.bkup") {
            push (@managed, $file);
        }
    }
    if ($combined_has_CDS) {
        unless (-e "$genome_lib_dir/ref_annot.prot_info.dbm" || -e "$genome_lib_dir/ref_annot.prot_info.dbm.bkup") {
            die "Error, supplements include CDS records but the lib has no ref_annot.prot_info.dbm";
        }
        foreach my $file (@PROT_INFO_MANAGED) {
            if (-e "$genome_lib_dir/$file" || -e "$genome_lib_dir/$file.bkup") {
                push (@managed, $file);
            }
        }
    }
    else {
        print STDERR "-supplements contain no CDS records; ref_annot.prot_info.dbm left unchanged.\n";
    }

    ## check required tools before any of the (large) copies are made
    if (grep { $_ eq "ref_annot.gtf.mm2.splice.bed" } @managed) {
        &require_tool("k8", "to update ref_annot.gtf.mm2.splice.bed (paftools)");
    }
    if ($combined_has_CDS) {
        &require_tool("samtools", "to update ref_annot.prot_info.dbm for supplements with CDS records");
    }

    &check_disk_space(\@managed, ($combined_has_CDS ? &prot_info_workspace_estimate($combined_gtf) : 0));

    my %tmp_outputs;
    foreach my $file (@managed) {
        my $tmp = "$genome_lib_dir/$file.tmp.$$";
        ## never write into an existing file: it could be a hard link to a lib file
        unlink($tmp);
        if (-e $tmp || -l $tmp) {
            die "Error, cannot remove $tmp";
        }
        $tmp_outputs{$file} = $tmp;
    }

    print STDERR "-building ref_annot.gtf\n";
    &append_files(&src_file("ref_annot.gtf"), $combined_gtf, $tmp_outputs{"ref_annot.gtf"});

    print STDERR "-building ref_annot.gtf.mini.sortu\n";
    &build_mini_sortu($combined_gtf, $tmp_outputs{"ref_annot.gtf.mini.sortu"});

    print STDERR "-building ref_annot.gtf.gene_spans\n";
    my $affected_gene_keys_href = &build_gene_spans($combined_gtf, $tmp_outputs{"ref_annot.gtf.gene_spans"});

    print STDERR "-building fusion_annot_lib.idx\n";
    &build_fusion_annot_idx($tmp_outputs{"ref_annot.gtf.gene_spans"}, $affected_gene_keys_href, $tmp_outputs{"fusion_annot_lib.idx"});

    if ($tmp_outputs{"ref_annot.gtf.mm2.splice.bed"}) {
        print STDERR "-building ref_annot.gtf.mm2.splice.bed\n";
        &build_mm2_splice_bed($combined_gtf, $tmp_outputs{"ref_annot.gtf.mm2.splice.bed"});
    }

    if ($combined_has_CDS) {
        print STDERR "-building ref_annot.prot_info.dbm\n";
        &build_prot_info($combined_gtf, \%tmp_outputs);
    }

    ## guard, before anything in the lib changes: each replacement must be a new file, sharing
    ## storage with neither the current lib file nor its original
    foreach my $file (@managed) {
        my $tmp_inode = (stat($tmp_outputs{$file}))[1];
        foreach my $existing ("$genome_lib_dir/$file", "$genome_lib_dir/$file.bkup") {
            my $inode = (stat($existing))[1];
            if (defined $inode && $tmp_inode == $inode) {
                die "Error, replacement for $file shares storage with $existing; refusing to commit";
            }
        }
    }

    ######################################
    ## commit: back up originals (first time only), record the supplement, then move the new files into place

    print STDERR "-committing changes to $genome_lib_dir\n";

    foreach my $file (@managed) {
        my $lib_file = "$genome_lib_dir/$file";
        my $bkup = "$lib_file.bkup";
        unless (-e $bkup) {
            my $bkup_tmp = "$bkup.tmp.$$";
            unlink($bkup_tmp);
            ## A hard link is safe here because lib files are never written in place: each is
            ## replaced by renaming a separately built file over it, leaving .bkup the only link
            ## to the original. Otherwise copy-on-write where supported, or a regular copy.
            if (-l $lib_file) {
                ## a symlinked file in a linked copy: the backup points at the same original
                symlink(readlink($lib_file), $bkup_tmp) or die "Error, cannot create symlink $bkup_tmp: $!";
                print STDERR "\t-backing up original $file -> $file.bkup (SYMLINK to the original)\n";
            }
            elsif (link($lib_file, $bkup_tmp)) {
                print STDERR "\t-backing up original $file -> $file.bkup (HARD LINK to the original; no extra disk space)\n";
                $num_hardlinked_backups++;
            }
            else {
                print STDERR "\t-backing up original $file -> $file.bkup (copy)\n";
                my $ret = system("cp -p --reflink=auto " . &shq($lib_file) . " " . &shq($bkup_tmp) . " 2>/dev/null");
                if ($ret) {
                    &process_cmd("cp -p " . &shq($lib_file) . " " . &shq($bkup_tmp));
                }
            }
            rename($bkup_tmp, $bkup) or die "Error, cannot rename $bkup_tmp to $bkup";
            $manifest->{backups}->{$file} = &md5_file($bkup);
            ## record each backup as soon as it exists
            &write_manifest($manifest);
        }
        unless ($manifest->{backups}->{$file}) {
            die "Error, $bkup exists but is not recorded in $MANIFEST";
        }
    }

    unless ($resume_flag) {
        copy($supp_gtf, "$SUPP_DIR/$supp_id.$label.input.gtf") or die "Error, cannot copy $supp_gtf to $SUPP_DIR: $!";
        rename($normalized_gtf, "$SUPP_DIR/$supp_id.$label.gtf") or die "Error, cannot move $normalized_gtf to $SUPP_DIR";

        ## recorded before the lib files are replaced, so that an interrupted commit is completed by
        ## the next run rather than silently dropped
        push (@{$manifest->{supplements}}, { id => $supp_id,
                                             label => $label,
                                             status => "committing",
                                             added => scalar(localtime()),
                                             input_gtf => $supp_gtf,
                                             input_gtf_copy => "$supp_id.$label.input.gtf",
                                             normalized_gtf => "$supp_id.$label.gtf",
                                             options => { assign_by_overlap => $assign_by_overlap_flag,
                                                          allow_novel_genes => $allow_novel_genes_flag,
                                                          allow_nonoverlapping => $allow_nonoverlapping_flag,
                                                          novel_gene_type => $novel_gene_type },
                                             num_transcripts => $report->{num_transcripts_kept},
                                             num_existing_genes => scalar(keys %{$report->{existing_genes}}),
                                             num_earlier_supplement_genes => scalar(keys %{$report->{earlier_supplement_genes}}),
                                             num_novel_genes => scalar(keys %{$report->{novel_genes}}),
                                           } );
    }
    foreach my $supp (@{$manifest->{supplements}}) {
        $supp->{status} = "committing";
    }
    &write_manifest($manifest);

    foreach my $file (@managed) {
        rename($tmp_outputs{$file}, "$genome_lib_dir/$file") or die "Error, cannot rename $tmp_outputs{$file} to $genome_lib_dir/$file";
    }

    foreach my $supp (@{$manifest->{supplements}}) {
        $supp->{status} = "complete";
    }
    $manifest->{modified_files} = [ sort keys %{$manifest->{backups}} ];
    $manifest->{last_updated} = scalar(localtime());
    &write_manifest($manifest);

    rmtree($WORKDIR);

    print STDERR "\n-done. Genome lib supplemented with '$label' ($supp_id).\n"
        . "  Modified files (originals preserved as .bkup): " . join(", ", @managed) . "\n"
        . "  Supplements are recorded in $MANIFEST\n"
        . ($source_lib_dir ? "  Linked copy of $source_lib_dir (unchanged)\n" : "")
        . "  To revert:  $UTILDIR/revert_ctat_genome_lib_supplement.pl --genome_lib_dir $genome_lib_dir\n\n";

    if ($source_lib_dir) {
        if ($num_copy_hardlinks) {
            print STDERR "-WARNING: unmodified files in $genome_lib_dir are HARD LINKS to the files in\n"
                . "  $source_lib_dir" . ($num_copy_symlinks ? " (except $num_copy_symlinks, which are symlinks)" : "") . "\n"
                . $HARDLINK_RULES . "\n";
        }
        if ($num_copy_symlinks) {
            print STDERR "-WARNING: $num_copy_symlinks files in $genome_lib_dir are SYMLINKS to files in\n"
                . "  $source_lib_dir: the copy only works where that path is visible\n"
                . "  (eg. in a container, only if it is also mounted at that same path).\n\n";
        }
    }
    elsif ($num_hardlinked_backups) {
        print STDERR "-WARNING: the new .bkup files were made as HARD LINKS to the original lib files (no extra\n"
            . "  disk space), so each .bkup is now the original file itself, and also shares storage with\n"
            . "  any linked copies made from this lib. Do not edit or replace .bkup files: revert restores them.\n\n";
    }

    exit(0);
}


END {
    if ($created_copy_dir && $? != 0 && -d $created_copy_dir) {
        ## don't leave behind an unsupplemented (or partly supplemented) linked copy
        print STDERR "-removing incomplete linked copy $created_copy_dir\n";
        rmtree($created_copy_dir);
        $have_lock = 0;
    }
    elsif ($WORKDIR && -d $WORKDIR && $? != 0) {
        if ($terminated_by_signal) {
            rmtree($WORKDIR);
        }
        else {
            print STDERR "-work files retained for inspection: $WORKDIR\n";
        }
    }
    if ($have_lock) {
        if ($? != 0) {
            ## remove any partially-built replacement files.
            if (opendir(my $dh, $genome_lib_dir)) {
                foreach my $entry (readdir($dh)) {
                    if ($entry =~ /\.tmp\.(\d+)$/ && $1 == $$) {
                        unlink("$genome_lib_dir/$entry");
                    }
                }
                closedir($dh);
            }
        }
        &release_lock();
    }
    ## consistent failure status (die otherwise exits with whatever errno was last set)
    if ($? != 0) {
        $? = 1;
    }
}


####
sub create_linked_copy {
    my ($source_dir, $copy_dir) = @_;

    ## resolve symlinks, so the checks below compare real locations
    unless (-d $source_dir) {
        die "Error, cannot locate genome lib dir: $source_dir";
    }
    $source_dir = abs_path($source_dir);
    $copy_dir =~ s{/+$}{};
    if (-e $copy_dir || -l $copy_dir) {
        unless (-d $copy_dir) {
            die "Error, --linked_copy target $copy_dir exists and is not a directory\n";
        }
        $copy_dir = abs_path($copy_dir);
    }
    else {
        my $parent = dirname($copy_dir);
        unless (-d $parent) {
            die "Error, parent directory of --linked_copy target does not exist: $parent\n";
        }
        $copy_dir = abs_path($parent) . "/" . basename($copy_dir);
    }

    ## a new lib whose files are hard links to those of the source lib; supplementing it
    ## replaces only its own links, leaving the source lib unchanged.

    if ($copy_dir eq $source_dir || index("$copy_dir/", "$source_dir/") == 0) {
        die "Error, --linked_copy target must be outside $source_dir\n";
    }
    if (-d "$source_dir/supplemental/.lock") {
        die "Error, $source_dir/supplemental/.lock exists: a supplement or revert may be running on it.\n";
    }

    ## claim the target: create it, or use it if it is an existing empty directory, then create a
    ## marker exclusively so that two runs can never both build (or later remove) the same copy.
    unless (-d $copy_dir) {
        mkdir($copy_dir) or die "Error, cannot mkdir $copy_dir: $!";
    }
    opendir(my $dh, $copy_dir) or die "Error, cannot read $copy_dir: $!";
    my @existing = grep { !/^\.\.?$/ } readdir($dh);
    closedir($dh);
    if (@existing) {
        die "Error, --linked_copy target $copy_dir exists and is not empty\n";
    }
    my $claim = "$copy_dir/.linked_copy_in_progress";
    sysopen(my $claim_fh, $claim, O_WRONLY | O_CREAT | O_EXCL)
        or die "Error, cannot claim $copy_dir (another run may be creating a copy there): $!\n";
    print $claim_fh hostname() . " $$\n";
    close $claim_fh;

    my $probe = "$copy_dir/.hardlink_probe.$$";
    my $can_hardlink = link("$source_dir/ref_annot.gtf", $probe);
    unlink($probe);

    if ($can_hardlink) {
        print STDERR "\n-WARNING: creating a linked copy of $source_dir\n"
            . "  at $copy_dir using HARD LINKS: the copy's files will be the same files as the source's\n"
            . "  (no extra disk space) until the supplement replaces them.\n"
            . $HARDLINK_RULES . "\n";
    }
    else {
        print STDERR "\n-WARNING: hard links to $source_dir are not possible here (different filesystem,\n"
            . "  or files not owned by you), so the linked copy at $copy_dir will use SYMLINKS.\n"
            . "  It will only work where $source_dir is visible at that same path\n"
            . "  (eg. in a container, only if it is also mounted at that path).\n\n";
    }

    my @symlinked;
    my $num_linked = 0;

    my $copy_entry;
    $copy_entry = sub {
        my ($rel) = @_;
        my $src = ($rel eq "") ? $source_dir : "$source_dir/$rel";
        my $dest = ($rel eq "") ? $copy_dir : "$copy_dir/$rel";

        if (-l $src) {
            ## keep the lib's own symlinks as they are (eg. relative links within the lib)
            symlink(readlink($src), $dest) or die "Error, cannot create symlink $dest: $!";
        }
        elsif (-d $src) {
            unless (-d $dest) {
                mkdir($dest) or die "Error, cannot mkdir $dest: $!";
            }
            opendir(my $dh, $src) or die "Error, cannot read dir $src: $!";
            my @entries = sort grep { !/^\.\.?$/ } readdir($dh);
            closedir($dh);
            foreach my $entry (@entries) {
                my $entry_rel = ($rel eq "") ? $entry : "$rel/$entry";
                ## skip transient state of an in-progress supplement of the source lib
                if ($entry_rel =~ m{^supplemental/\.(work|lock)} || $entry_rel =~ /\.tmp\.[\d\.]+$/
                    || $entry_rel eq ".linked_copy_in_progress") { next; }
                $copy_entry->($entry_rel);
            }
        }
        elsif (link($src, $dest)) {
            $num_linked++;
        }
        else {
            symlink($src, $dest) or die "Error, cannot link or symlink $src to $dest: $!";
            push (@symlinked, $rel);
        }
    };

    ## from here on, this run owns $copy_dir: END removes it if the run fails
    $created_copy_dir = $copy_dir;

    $copy_entry->("");
    unlink($claim);

    $num_copy_hardlinks = $num_linked;
    $num_copy_symlinks = scalar(@symlinked);

    print STDERR "-linked copy created: $num_linked files hard linked, " . scalar(@symlinked) . " symlinked\n";
    if (@symlinked) {
        print STDERR "\n-WARNING: " . scalar(@symlinked) . " files could not be hard linked (different filesystem, or\n"
            . "  files not owned by you) and are symlinks to $source_dir instead, eg. $symlinked[0]\n"
            . "  The copy will only work where that path is visible, eg. in a container only if\n"
            . "  $source_dir is also mounted at that same path.\n\n";
    }

    return($copy_dir);
}


####
sub acquire_lock {

    unless (-d $SUPP_DIR) {
        mkpath($SUPP_DIR);
        unless (-d $SUPP_DIR) {
            die "Error, cannot mkdir $SUPP_DIR";
        }
    }
    unless (mkdir($LOCK_DIR)) {
        my $owner = "";
        if (open(my $fh, "$LOCK_DIR/owner")) {
            $owner = <$fh>;
            close $fh;
            chomp $owner if defined $owner;
        }
        die "Error, $LOCK_DIR exists: another supplement or revert may be running"
            . ($owner ? " ($owner)" : "") . ".\n"
            . "  If not (eg. that process no longer exists), remove that directory and try again.\n";
    }
    $have_lock = 1;
    if (open(my $ofh, ">$LOCK_DIR/owner")) {
        print $ofh "supplement_ctat_genome_lib.pl, host " . hostname() . ", pid $$, started " . scalar(localtime()) . "\n";
        close $ofh;
    }
    return;
}


####
sub release_lock {
    unlink("$LOCK_DIR/owner");
    rmdir($LOCK_DIR);
    $have_lock = 0;
    return;
}


####
sub remove_stale_files {

    ## Left by an earlier run that was killed. Removed only while holding the lock, so no
    ## running supplement or revert can own them. Replacement files are never written into an
    ## existing file, so this also guarantees a stale one can't be a link to a lib file.

    my %managed = map { $_ => 1 } (@REQUIRED_MANAGED, @OPTIONAL_MANAGED, @PROT_INFO_MANAGED);

    opendir(my $dh, $genome_lib_dir) or die "Error, cannot read $genome_lib_dir: $!";
    foreach my $entry (readdir($dh)) {
        if ($entry =~ /^(.+?)(\.bkup)?\.tmp\.\d+$/ && $managed{$1}) {
            print STDERR "-removing stale temporary file from an interrupted run: $entry\n";
            unlink("$genome_lib_dir/$entry") or die "Error, cannot remove $genome_lib_dir/$entry: $!";
        }
    }
    closedir($dh);

    opendir($dh, $SUPP_DIR) or die "Error, cannot read $SUPP_DIR: $!";
    foreach my $entry (readdir($dh)) {
        if ($entry =~ /^\.work\.\d+$/) {
            print STDERR "-removing stale work directory from an interrupted run: supplemental/$entry\n";
            rmtree("$SUPP_DIR/$entry");
        }
        elsif ($entry =~ /^manifest\.json\.tmp\.\d+$/) {
            unlink("$SUPP_DIR/$entry");
        }
    }
    closedir($dh);

    return;
}


####
sub remove_orphan_supplement_files {
    my ($manifest) = @_;

    ## supplement gtfs not recorded in the manifest: left by a run interrupted just before
    ## recording its supplement (the lib files were not yet changed)
    my %recorded;
    foreach my $supp (@{$manifest->{supplements}}) {
        $recorded{$supp->{normalized_gtf}} = 1;
        $recorded{$supp->{input_gtf_copy}} = 1;
    }
    opendir(my $dh, $SUPP_DIR) or die "Error, cannot read $SUPP_DIR: $!";
    foreach my $entry (readdir($dh)) {
        if ($entry =~ /^supp_\d+\..*\.gtf$/ && ! $recorded{$entry}) {
            print STDERR "-removing unrecorded supplement file from an interrupted run: supplemental/$entry\n";
            unlink("$SUPP_DIR/$entry");
        }
    }
    closedir($dh);

    return;
}


####
sub require_tool {
    my ($tool, $purpose) = @_;

    foreach my $dir (split(/:/, $ENV{PATH} || "")) {
        if (-x "$dir/$tool" && ! -d "$dir/$tool") {
            return;
        }
    }
    die "Error, '$tool' is required $purpose, but was not found in PATH; no changes made.\n";
}


####
sub shq {
    my ($str) = @_;
    ## single-quote for the shell
    $str =~ s/'/'\\''/g;
    return("'$str'");
}


####
sub src_file {
    my ($file) = @_;
    ## always build from the original version of a lib file.
    my $lib_file = "$genome_lib_dir/$file";
    if (-e "$lib_file.bkup") {
        return("$lib_file.bkup");
    }
    return($lib_file);
}


####
sub read_manifest {

    if (-e $MANIFEST) {
        open(my $fh, $MANIFEST) or die "Error, cannot read $MANIFEST";
        my $json = do { local $/; <$fh> };
        close $fh;
        return(JSON::PP->new->decode($json));
    }

    return( { supplements => [],
              backups => {},
              modified_files => [],
            } );
}


####
sub write_manifest {
    my ($manifest) = @_;

    my $json = JSON::PP->new->pretty->canonical->encode($manifest);
    open(my $ofh, ">$MANIFEST.tmp.$$") or die "Error, cannot write $MANIFEST.tmp.$$";
    print $ofh $json;
    close $ofh;
    rename("$MANIFEST.tmp.$$", $MANIFEST) or die "Error, cannot rename $MANIFEST.tmp.$$ to $MANIFEST";

    return;
}


####
sub check_backup_state {
    my ($manifest) = @_;

    ## every .bkup must be one this utility created and recorded.
    foreach my $file (@REQUIRED_MANAGED, @OPTIONAL_MANAGED, @PROT_INFO_MANAGED) {
        my $lib_file = "$genome_lib_dir/$file";
        my $bkup = "$lib_file.bkup";
        if (-e $bkup && ! $manifest->{backups}->{$file}) {
            ## A run interrupted right after making the backup, before recording it: the lib file
            ## had not been replaced yet, so the two are identical and the backup is the original.
            if (-e $lib_file && &same_content($lib_file, $bkup)) {
                print STDERR "-recording backup $file.bkup made by an interrupted run (identical to $file)\n";
                $manifest->{backups}->{$file} = &md5_file($bkup);
                &write_manifest($manifest) unless $dry_run_flag;
                next;
            }
            die "Error, found $bkup but it is not recorded in $MANIFEST, and it differs from $file.\n"
                . "  It was not made by this utility (or the lib was changed by other means); resolve this before supplementing.\n";
        }
        if ($manifest->{backups}->{$file} && ! -e $bkup) {
            die "Error, $MANIFEST records a backup of $file but $bkup is missing.\n"
                . "  If a revert was interrupted, run revert_ctat_genome_lib_supplement.pl again to complete it.\n";
        }
    }
    foreach my $supp (@{$manifest->{supplements}}) {
        my $gtf = "$SUPP_DIR/$supp->{normalized_gtf}";
        unless (-s $gtf) {
            die "Error, $MANIFEST lists supplement $supp->{id} but $gtf is missing";
        }
    }

    return;
}


####
sub same_content {
    my ($file_a, $file_b) = @_;

    my @stat_a = stat($file_a);
    my @stat_b = stat($file_b);
    if ($stat_a[0] == $stat_b[0] && $stat_a[1] == $stat_b[1]) {
        return(1); # the same file
    }
    if ($stat_a[7] != $stat_b[7]) {
        return(0);
    }
    return(&md5_file($file_a) eq &md5_file($file_b));
}


####
sub check_disk_space {
    my ($managed_aref, $workspace_bytes) = @_;

    my $needed = $workspace_bytes || 0;
    foreach my $file (@$managed_aref) {
        my $size = -s &src_file($file);
        $needed += $size; # replacement file
        unless (-e "$genome_lib_dir/$file.bkup") {
            $needed += $size; # backup, in case it can't be a hard link
        }
    }

    my $df = `df -Pk @{[ &shq($genome_lib_dir) ]} | tail -1`;
    my @x = split(/\s+/, $df);
    unless (@x >= 4 && $x[3] =~ /^\d+$/) {
        print STDERR "-WARNING: could not determine free disk space in $genome_lib_dir\n";
        return;
    }
    my $avail = $x[3] * 1024;
    if ($avail < $needed * 1.1) {
        die sprintf("Error, insufficient disk space in %s: need ~%.2f GB, have %.2f GB\n",
                    $genome_lib_dir, $needed / 1e9, $avail / 1e9);
    }

    return;
}


####
sub prot_info_workspace_estimate {
    my ($combined_gtf) = @_;

    ## build_prot_info writes a genome subset of every chromosome with supplemental CDS records
    my %chrs;
    open(my $fh, $combined_gtf) or die "Error, cannot open $combined_gtf";
    while (<$fh>) {
        my @x = split(/\t/);
        if ($x[2] && $x[2] eq "CDS") {
            $chrs{$x[0]} = 1;
        }
    }
    close $fh;

    my $bytes = 0;
    if (open(my $fai, "$genome_lib_dir/ref_genome.fa.fai")) {
        while (<$fai>) {
            my ($chr, $len) = split(/\t/);
            $bytes += $len if $chrs{$chr};
        }
        close $fai;
    }

    return($bytes * 1.02);
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


####
sub parse_attributes {
    my ($info) = @_;

    my @atts;
    ## the final attribute's ';' is optional
    while ($info =~ /\s*([^\s;]+)\s+(\"[^\"]*\"|[^;]+?)\s*(?:;|$)/g) {
        my ($key, $val) = ($1, $2);
        $val =~ s/^\"|\"$//g;
        push (@atts, [$key, $val]);
    }

    return(@atts);
}


####
sub load_gene_spans {
    my ($gene_spans_file) = @_;

    my %genes;
    open(my $fh, $gene_spans_file) or die "Error, cannot open $gene_spans_file";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        my ($gene_id, $chr, $lend, $rend, $orient, $gene_name, $gene_type) = split(/\t/);
        ## gtf_to_gene_spans.pl lists a gene lacking a gene_name under its gene_id: it has no name
        if (! defined $gene_name || $gene_name eq "" || $gene_name eq $gene_id) {
            $gene_name = ".";
        }
        my $key = "$gene_id\t$chr";
        $genes{$key} = { gene_id => $gene_id,
                         chr => $chr,
                         lend => $lend,
                         rend => $rend,
                         orient => $orient,
                         gene_name => $gene_name,
                         gene_type => $gene_type,
                         novel => 0,
        };
    }
    close $fh;

    return(\%genes);
}


####
sub add_prior_supplement_genes {
    my ($genes_href, $prior_gtfs_aref) = @_;

    foreach my $gtf (@$prior_gtfs_aref) {
        open(my $fh, $gtf) or die "Error, cannot open $gtf";
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            unless ($x[2] eq "exon") { next; }
            my %att = map { $_->[0] => $_->[1] } &parse_attributes($x[8]);
            my $key = "$att{gene_id}\t$x[0]";
            if (my $gene = $genes_href->{$key}) {
                if ($gene->{novel}) {
                    $gene->{lend} = $x[3] if $x[3] < $gene->{lend};
                    $gene->{rend} = $x[4] if $x[4] > $gene->{rend};
                }
                next;
            }
            $genes_href->{$key} = { gene_id => $att{gene_id},
                                    chr => $x[0],
                                    lend => $x[3],
                                    rend => $x[4],
                                    orient => $x[6],
                                    gene_name => ($att{gene_name} && $att{gene_name} ne $att{gene_id}) ? $att{gene_name} : ".",
                                    gene_type => $att{gene_type} || ".",
                                    novel => 1,
            };
        }
        close $fh;
    }

    return;
}


####
sub build_gene_lookups {
    my ($genes_href) = @_;

    my %by_id;
    my %by_id_nover;
    my %by_name;

    foreach my $key (keys %$genes_href) {
        my $gene = $genes_href->{$key};
        push (@{$by_id{$gene->{gene_id}}}, $key);

        ## version-stripped ids only for reference genes: ids of genes added by supplements (eg.
        ## StringTie's STRG.100) are not versioned, and stripping them would conflate unrelated genes
        unless ($gene->{novel}) {
            my $id_nover = $gene->{gene_id};
            if ($id_nover =~ s/\.\d+$//) {
                push (@{$by_id_nover{$id_nover}}, $key);
            }
        }

        if ($gene->{gene_name} ne ".") {
            push (@{$by_name{$gene->{gene_name}}}, $key);
        }
    }

    return( { by_id => \%by_id,
              by_id_nover => \%by_id_nover,
              by_name => \%by_name } );
}


####
sub load_transcript_ids {
    my ($gtf, $tids_href, $source) = @_;

    open(my $fh, $gtf) or die "Error, cannot open $gtf";
    while (<$fh>) {
        if (/transcript_id \"([^\"]+)\"/) {
            $tids_href->{$1} = $source;
        }
    }
    close $fh;

    return;
}


####
sub load_ref_exons {
    my ($exon_record_files_aref, $genes_href) = @_;

    ## map NAME^ID tokens back to gene keys
    my %token_to_key;
    foreach my $key (keys %$genes_href) {
        my $gene = $genes_href->{$key};
        if ($gene->{gene_name} ne ".") {
            $token_to_key{"$gene->{gene_name}^$gene->{gene_id}\t$gene->{chr}"} = $key;
        }
        else {
            $token_to_key{"$gene->{gene_id}\t$gene->{chr}"} = $key;
            ## as written for novel genes by earlier versions of this utility
            $token_to_key{"$gene->{gene_id}^$gene->{gene_id}\t$gene->{chr}"} = $key;
        }
    }

    my %exons; # chr -> orient -> [ [lend, rend, gene_key], ... ] sorted by lend
    foreach my $exon_records_file (@$exon_record_files_aref) {
        open(my $fh, $exon_records_file) or die "Error, cannot open $exon_records_file";
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            unless ($x[8] =~ /gene_id \"([^\"]+)\"/) { next; }
            my $gene_key = $token_to_key{"$1\t$x[0]"} or next;
            push (@{$exons{$x[0]}->{$x[6]}}, [$x[3], $x[4], $gene_key]);
        }
        close $fh;
    }

    my %sorted_exons; # chr -> orient -> { exons => [ sorted by lend ], max_len => longest exon }
    foreach my $chr (keys %exons) {
        foreach my $orient (keys %{$exons{$chr}}) {
            my @sorted = sort { $a->[0] <=> $b->[0] } @{$exons{$chr}->{$orient}};
            my $max_len = 0;
            foreach my $exon (@sorted) {
                my $len = $exon->[1] - $exon->[0] + 1;
                $max_len = $len if $len > $max_len;
            }
            $sorted_exons{$chr}->{$orient} = { exons => \@sorted, max_len => $max_len };
        }
    }

    return(\%sorted_exons);
}


####
sub overlapping_gene_keys {
    my ($exon_set, $lend, $rend) = @_;

    my $exons_aref = $exon_set->{exons};

    ## binary search for the first exon that could reach $lend
    my $min_start = $lend - $exon_set->{max_len};
    my ($lo, $hi) = (0, scalar(@$exons_aref));
    while ($lo < $hi) {
        my $mid = int(($lo + $hi) / 2);
        if ($exons_aref->[$mid]->[0] < $min_start) {
            $lo = $mid + 1;
        }
        else {
            $hi = $mid;
        }
    }

    my @gene_keys;
    for (my $i = $lo; $i <= $#$exons_aref; $i++) {
        my $ref_exon = $exons_aref->[$i];
        if ($ref_exon->[0] > $rend) { last; }
        if ($ref_exon->[1] >= $lend) {
            push (@gene_keys, $ref_exon->[2]);
        }
    }

    return(@gene_keys);
}


####
sub parse_supplemental_gtf {
    my ($gtf) = @_;

    my %transcripts;
    my %report = ( skipped_unstranded => [],
                   dropped_gene_records => 0,
                   errors => [] );
    my $errors = $report{errors};
    my %gene_id_conflict;

    open(my $fh, $gtf) or die "Error, cannot open $gtf";
    read($fh, my $magic, 2);
    if (defined $magic && $magic eq "\x1f\x8b") {
        die "Error, $gtf is gzip-compressed; please provide it uncompressed.\n";
    }
    seek($fh, 0, 0);
    while (<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        chomp;
        s/\r$//;
        my @x = split(/\t/);
        unless (scalar(@x) == 9) {
            push (@$errors, "line $. does not have 9 tab-delimited fields: $_");
            next;
        }
        if ($x[2] eq "gene") {
            ## gene records are taken from the reference (or derived for novel genes).
            $report{dropped_gene_records}++;
            next;
        }
        unless ($x[3] =~ /^\d+$/ && $x[4] =~ /^\d+$/ && $x[3] <= $x[4]) {
            push (@$errors, "line $. has invalid coordinates: $_");
            next;
        }
        my @atts = &parse_attributes($x[8]);
        my %att = map { $_->[0] => $_->[1] } @atts;
        my $tid = $att{transcript_id};
        unless (defined $tid && $tid ne "") {
            push (@$errors, "line $. lacks a transcript_id: $_");
            next;
        }
        my $transcript = $transcripts{$tid};
        unless ($transcript) {
            $transcript = $transcripts{$tid} = { transcript_id => $tid,
                                                 chr => $x[0],
                                                 orient => $x[6],
                                                 atts => \%att,
                                                 lines => [],
                                                 exons => [],
                                                 order => scalar(keys %transcripts),
            };
        }
        if ($x[0] ne $transcript->{chr} || $x[6] ne $transcript->{orient}) {
            push (@$errors, "transcript $tid has records on more than one chromosome/strand (line $.)");
            $transcript->{bad} = 1;
            next;
        }
        ## the same transcript_id under different genes (eg. merged GTFs) must not be combined
        my $gene_id = $att{gene_id} // "";
        my $first_gene_id = $transcript->{atts}->{gene_id} // "";
        if ($gene_id ne "" && $first_gene_id ne "" && $gene_id ne $first_gene_id && ! $gene_id_conflict{$tid}) {
            push (@$errors, "transcript $tid has records with different gene_ids ('$first_gene_id', '$gene_id'; line $.)");
            $gene_id_conflict{$tid} = 1;
            $transcript->{bad} = 1;
        }
        foreach my $key (keys %att) {
            unless (exists $transcript->{atts}->{$key}) {
                $transcript->{atts}->{$key} = $att{$key};
            }
        }
        push (@{$transcript->{lines}}, [ @x[0..7], \@atts ]);
        if ($x[2] eq "exon") {
            push (@{$transcript->{exons}}, [$x[3], $x[4]]);
        }
    }
    close $fh;

    foreach my $tid (sort { $transcripts{$a}->{order} <=> $transcripts{$b}->{order} } keys %transcripts) {
        my $transcript = $transcripts{$tid};
        if ($transcript->{bad}) {
            delete $transcripts{$tid};
            next;
        }
        unless (@{$transcript->{exons}}) {
            push (@$errors, "transcript $tid has no exon records");
            delete $transcripts{$tid};
            next;
        }
        if ($transcript->{orient} ne '+' && $transcript->{orient} ne '-') {
            push (@{$report{skipped_unstranded}}, $tid);
            delete $transcripts{$tid};
        }
    }

    return(\%transcripts, \%report);
}


####
sub assign_and_normalize {
    my ($transcripts_href, $genes_href, $lookup, $existing_tids_href, $ref_exons_href) = @_;

    my %report = ( errors => [],
                   warnings => [],
                   num_transcripts_kept => 0,
                   existing_genes => {},
                   earlier_supplement_genes => {},
                   novel_genes => {},
                   assigned_by => {},
        );

    my %novel_gene_info; # gene_id -> { chr, orient }
    my %gene_to_transcripts;

    foreach my $tid (sort { $transcripts_href->{$a}->{order} <=> $transcripts_href->{$b}->{order} } keys %$transcripts_href) {
        my $transcript = $transcripts_href->{$tid};
        my $chr = $transcript->{chr};
        my $orient = $transcript->{orient};
        my $atts = $transcript->{atts};

        if ($tid =~ /[\|\^\s]/) {
            push (@{$report{errors}}, "transcript_id '$tid' contains a disallowed character ('|', '^' or whitespace)");
            next;
        }
        if (my $source = $existing_tids_href->{$tid}) {
            push (@{$report{errors}}, "transcript_id $tid already exists in the $source annotation");
            next;
        }

        my ($gene_key, $method, $problem) = &assign_gene($transcript, $genes_href, $lookup, $ref_exons_href);

        if ($problem) {
            push (@{$report{errors}}, "transcript $tid: $problem");
            next;
        }

        if ($gene_key) {
            my $gene = $genes_href->{$gene_key};
            if ($gene->{chr} ne $chr || $gene->{orient} ne $orient) {
                push (@{$report{errors}}, "transcript $tid ($chr:$orient) assigned by $method to gene "
                      . &gene_label($gene) . " which is on $gene->{chr}:$gene->{orient}");
                next;
            }
            ## catch likely mismatches (eg. a different genome build, or a reused id or name)
            my @coords = sort { $a <=> $b } map { @$_ } @{$transcript->{exons}};
            if (defined $gene->{lend} && ($coords[-1] < $gene->{lend} || $coords[0] > $gene->{rend})) {
                my $msg = "transcript $tid ($chr:$coords[0]-$coords[-1]) does not overlap the gene it was assigned to by $method, "
                    . &gene_label($gene) . " ($gene->{chr}:$gene->{lend}-$gene->{rend})";
                if ($allow_nonoverlapping_flag) {
                    push (@{$report{warnings}}, "$msg; the gene span will be extended to include it");
                }
                else {
                    push (@{$report{errors}}, "$msg (see --allow_nonoverlapping)");
                    next;
                }
            }

            $report{assigned_by}->{$method}++;
            if ($gene->{novel}) {
                $report{earlier_supplement_genes}->{$gene_key} = 1;
            }
            else {
                $report{existing_genes}->{$gene_key} = 1;
            }
        }
        else {
            unless ($allow_novel_genes_flag) {
                push (@{$report{errors}}, "transcript $tid (gene_id: " . ($atts->{gene_id} || "none") . ") "
                      . "could not be assigned to a reference gene (see --assign_by_overlap, --allow_novel_genes)");
                next;
            }
            my $gene_id = $atts->{gene_id};
            unless ($gene_id) {
                push (@{$report{errors}}, "transcript $tid has no gene_id to use for a novel gene");
                next;
            }
            if ($gene_id =~ /[\|\^\s]/ || $gene_id =~ /--/) {
                push (@{$report{errors}}, "novel gene_id '$gene_id' (transcript $tid) contains a disallowed character ('|', '^', '--' or whitespace)");
                next;
            }
            if ($lookup->{by_name}->{$gene_id}) {
                push (@{$report{errors}}, "novel gene_id '$gene_id' (transcript $tid) collides with an existing gene_name");
                next;
            }
            if (my $info = $novel_gene_info{$gene_id}) {
                if ($info->{chr} ne $chr || $info->{orient} ne $orient) {
                    push (@{$report{errors}}, "novel gene $gene_id has transcripts on more than one chromosome/strand");
                    next;
                }
            }
            $novel_gene_info{$gene_id} = { chr => $chr, orient => $orient };
            $gene_key = "$gene_id\t$chr";
            $genes_href->{$gene_key} = { gene_id => $gene_id,
                                         chr => $chr,
                                         orient => $orient,
                                         gene_name => ".", # listed under its gene_id, as for any gene without a name
                                         gene_type => $novel_gene_type,
                                         novel => 1,
            };
            $report{novel_genes}->{$gene_id} = 1;
            $report{assigned_by}->{"novel gene"}++;
        }

        push (@{$gene_to_transcripts{$gene_key}}, $tid);
    }

    ## write normalized records, grouped by gene and transcript
    my @lines;
    foreach my $gene_key (sort keys %gene_to_transcripts) {
        my $gene = $genes_href->{$gene_key};
        foreach my $tid (@{$gene_to_transcripts{$gene_key}}) {
            my $transcript = $transcripts_href->{$tid};
            $report{num_transcripts_kept}++;

            my @recs = sort { $a->[3] <=> $b->[3]
                                  || ($a->[2] eq "transcript" ? 0 : 1) <=> ($b->[2] eq "transcript" ? 0 : 1)
                                  || $a->[2] cmp $b->[2] } @{$transcript->{lines}};

            foreach my $rec (@recs) {
                my @x = @$rec[0..7];
                $x[1] = $SOURCE_TAG;
                $x[8] = &normalized_attributes($gene, $tid, $rec->[8]);
                push (@lines, join("\t", @x));
            }
        }
    }

    return(\@lines, \%report);
}


####
sub assign_gene {
    my ($transcript, $genes_href, $lookup, $ref_exons_href) = @_;

    my $chr = $transcript->{chr};
    my $orient = $transcript->{orient};
    my $atts = $transcript->{atts};

    my $pick = sub {
        my ($keys_aref) = @_;
        ## prefer candidates on the transcript's chromosome
        my @on_chr = grep { $genes_href->{$_}->{chr} eq $chr } @$keys_aref;
        if (scalar(@on_chr) == 1) {
            return($on_chr[0]);
        }
        elsif (scalar(@$keys_aref) == 1) {
            return($keys_aref->[0]);
        }
        return(undef);
    };

    ## by gene_id, then ref_gene_id
    foreach my $att ("gene_id", "ref_gene_id") {
        my $val = $atts->{$att} or next;
        if (my $keys_aref = $lookup->{by_id}->{$val}) {
            my $key = &$pick($keys_aref);
            return($key, $att) if $key;
            my @chrs = sort map { $genes_href->{$_}->{chr} } @$keys_aref;
            return(undef, undef, "$att '$val' matches genes on more than one other chromosome (@chrs), none on $chr");
        }
    }

    ## by version-stripped gene_id / ref_gene_id (reference genes only), requiring chr and strand agreement
    foreach my $att ("gene_id", "ref_gene_id") {
        my $val = $atts->{$att} or next;
        $val =~ s/\.\d+$//;
        if (my $keys_aref = $lookup->{by_id_nover}->{$val}) {
            my @candidates = grep { $genes_href->{$_}->{chr} eq $chr
                                        && $genes_href->{$_}->{orient} eq $orient } @$keys_aref;
            if (scalar(@candidates) == 1) {
                return($candidates[0], "$att (version-stripped)");
            }
            elsif (scalar(@candidates) > 1) {
                return(undef, undef, "$att '$atts->{$att}' (version-stripped) matches more than one reference gene on $chr:$orient");
            }
            ## otherwise fall through to gene_name matching
        }
    }

    ## by gene name, requiring chr and strand agreement
    foreach my $att ("gene_name", "ref_gene_name") {
        my $val = $atts->{$att} or next;
        if (my $keys_aref = $lookup->{by_name}->{$val}) {
            my @candidates = grep { $genes_href->{$_}->{chr} eq $chr
                                        && $genes_href->{$_}->{orient} eq $orient } @$keys_aref;
            if (scalar(@candidates) == 1) {
                return($candidates[0], $att);
            }
            elsif (scalar(@candidates) > 1) {
                return(undef, undef, "$att '$val' matches more than one reference gene on $chr:$orient");
            }
            else {
                return(undef, undef, "$att '$val' matches a reference gene, but not on $chr:$orient");
            }
        }
    }

    ## by exon overlap
    if ($ref_exons_href) {
        my %overlapping_genes;
        if (my $exon_set = $ref_exons_href->{$chr}->{$orient}) {
            foreach my $exon (@{$transcript->{exons}}) {
                foreach my $gene_key (&overlapping_gene_keys($exon_set, @$exon)) {
                    $overlapping_genes{$gene_key} = 1;
                }
            }
        }
        my @candidates = keys %overlapping_genes;
        if (scalar(@candidates) == 1) {
            return($candidates[0], "exon overlap");
        }
        elsif (scalar(@candidates) > 1) {
            my @names = sort map { &gene_label($genes_href->{$_}) } @candidates;
            return(undef, undef, "exons overlap more than one reference gene (@names)");
        }
    }

    return(undef);
}


####
sub gene_label {
    my ($gene) = @_;
    return( ($gene->{gene_name} ne ".") ? "$gene->{gene_name} ($gene->{gene_id})" : $gene->{gene_id} );
}


####
sub normalized_attributes {
    my ($gene, $tid, $atts_aref) = @_;

    my @parts = ("gene_id \"$gene->{gene_id}\"",
                 "transcript_id \"$tid\"");
    if ($gene->{gene_name} ne ".") {
        push (@parts, "gene_name \"$gene->{gene_name}\"");
    }
    if ($gene->{gene_type} ne ".") {
        push (@parts, "gene_type \"$gene->{gene_type}\"");
    }

    my %skip = map { $_ => 1 } qw(gene_id transcript_id gene_name gene_type gene_biotype
                                  ref_gene_id ref_gene_name);
    my $have_transcript_type = 0;
    foreach my $att (@$atts_aref) {
        my ($key, $val) = @$att;
        if ($skip{$key}) { next; }
        if ($key eq "transcript_type") { $have_transcript_type = 1; }
        push (@parts, "$key \"$val\"");
    }
    unless ($have_transcript_type) {
        push (@parts, "transcript_type \"supplemental\"");
    }

    return(join("; ", @parts) . ";");
}


####
sub print_report {
    my ($report) = @_;

    print STDERR "\n#### supplemental annotation summary ####\n";
    print STDERR "  transcripts kept: $report->{num_transcripts_kept}\n";
    print STDERR "  existing genes receiving isoforms: " . scalar(keys %{$report->{existing_genes}}) . "\n";
    if (my $num = scalar(keys %{$report->{earlier_supplement_genes}})) {
        print STDERR "  genes from earlier supplements receiving isoforms: $num\n";
    }
    print STDERR "  novel genes: " . scalar(keys %{$report->{novel_genes}}) . "\n";
    foreach my $method (sort keys %{$report->{assigned_by}}) {
        print STDERR "    assigned by $method: $report->{assigned_by}->{$method}\n";
    }
    if (my $num = scalar(@{$report->{skipped_unstranded}})) {
        print STDERR "  skipped (no strand): $num\n";
    }
    if ($report->{dropped_gene_records}) {
        print STDERR "  'gene' records ignored: $report->{dropped_gene_records}\n";
    }
    if (my @warnings = @{$report->{warnings}}) {
        print STDERR "  WARNINGS: " . scalar(@warnings) . "\n";
        foreach my $warning (@warnings[0 .. ($#warnings < 24 ? $#warnings : 24)]) {
            print STDERR "    $warning\n";
        }
        print STDERR "    ...\n" if (@warnings > 25);
    }
    print STDERR "  problems: " . scalar(@{$report->{errors}}) . "\n\n";

    return;
}


####
sub write_lines {
    my ($file, $lines_aref) = @_;

    open(my $ofh, ">$file") or die "Error, cannot write $file";
    foreach my $line (@$lines_aref) {
        print $ofh "$line\n";
    }
    close $ofh;

    return;
}


####
sub gtf_has_feature {
    my ($gtf, $feature) = @_;

    open(my $fh, $gtf) or die "Error, cannot open $gtf";
    while (<$fh>) {
        my @x = split(/\t/);
        if ($x[2] && $x[2] eq $feature) {
            close $fh;
            return(1);
        }
    }
    close $fh;

    return(0);
}


####
sub append_files {
    my ($orig_file, $append_file, $out_file) = @_;

    &copy_with_trailing_newline($orig_file, $out_file);
    &process_cmd("cat " . &shq($append_file) . " >> " . &shq($out_file));

    return;
}


####
sub copy_with_trailing_newline {
    my ($orig_file, $out_file) = @_;

    &process_cmd("cp " . &shq($orig_file) . " " . &shq($out_file));

    ## so appended records don't run onto the last original line
    if (-s $out_file) {
        open(my $fh, $out_file) or die "Error, cannot open $out_file";
        seek($fh, -1, 2);
        read($fh, my $last_char, 1);
        close $fh;
        if ($last_char ne "\n") {
            open(my $ofh, ">>$out_file") or die "Error, cannot append to $out_file";
            print $ofh "\n";
            close $ofh;
        }
    }

    return;
}


####
sub build_mini_sortu {
    my ($combined_gtf, $out_file) = @_;

    my $supp_records = "$WORKDIR/all_supplements.mini";
    &process_cmd(&shq("$UTILDIR/gtf_to_exon_gene_records.pl") . " " . &shq($combined_gtf)
                 . " | sort -k 1,1 -k4,4g -k5,5g | uniq > " . &shq($supp_records));

    ## an exon is identified as the reader does: by gene, coordinates and strand (source and score ignored)
    my $exon_key = sub {
        my @x = split(/\t/, $_[0]);
        return(join("\t", @x[0,3,4,6,8]));
    };

    my %supp_rows; # exon key -> 1 while still to be appended
    my @supp_rows_ordered;
    open(my $fh, $supp_records) or die "Error, cannot open $supp_records";
    while (<$fh>) {
        my $key = &$exon_key($_);
        unless (exists $supp_rows{$key}) {
            push (@supp_rows_ordered, [$key, $_]);
            $supp_rows{$key} = 1;
        }
    }
    close $fh;

    ## only append exon records not already in the reference.
    my $src = &src_file("ref_annot.gtf.mini.sortu");
    open($fh, $src) or die "Error, cannot open $src";
    while (<$fh>) {
        my $key = &$exon_key($_);
        if (exists $supp_rows{$key}) {
            $supp_rows{$key} = 0;
        }
    }
    close $fh;

    &copy_with_trailing_newline($src, $out_file);
    open(my $ofh, ">>$out_file") or die "Error, cannot append to $out_file";
    foreach my $row (@supp_rows_ordered) {
        my ($key, $line) = @$row;
        if ($supp_rows{$key}) {
            print $ofh $line;
        }
    }
    close $ofh;

    return;
}


####
sub build_gene_spans {
    my ($combined_gtf, $out_file) = @_;

    ## gather supplemental exon extents per gene
    my %supp_spans; # gene_id \t chr -> { lend, rend, orient, gene_name, gene_type }
    my @novel_order;
    open(my $fh, $combined_gtf) or die "Error, cannot open $combined_gtf";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        unless ($x[2] eq "exon") { next; }
        my %att = map { $_->[0] => $_->[1] } &parse_attributes($x[8]);
        my $key = "$att{gene_id}\t$x[0]";
        my $span = $supp_spans{$key};
        unless ($span) {
            $span = $supp_spans{$key} = { gene_id => $att{gene_id},
                                          chr => $x[0],
                                          lend => $x[3],
                                          rend => $x[4],
                                          orient => $x[6],
                                          ## as gtf_to_gene_spans.pl: a gene without a gene_name is listed under its gene_id
                                          gene_name => $att{gene_name} || $att{gene_id},
                                          gene_type => $att{gene_type} || ".",
            };
            push (@novel_order, $key);
        }
        $span->{lend} = $x[3] if $x[3] < $span->{lend};
        $span->{rend} = $x[4] if $x[4] > $span->{rend};
    }
    close $fh;

    my %affected_gene_keys; # fusion_annot_lib.idx keys whose COORDS may have changed

    my $src = &src_file("ref_annot.gtf.gene_spans");
    open($fh, $src) or die "Error, cannot open $src";
    open(my $ofh, ">$out_file") or die "Error, cannot write $out_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $key = "$x[0]\t$x[1]";
        if (my $span = $supp_spans{$key}) {
            $x[2] = $span->{lend} if $span->{lend} < $x[2];
            $x[3] = $span->{rend} if $span->{rend} > $x[3];
            $span->{in_reference} = 1;
            $affected_gene_keys{ ($x[5] && $x[5] ne ".") ? $x[5] : $x[0] } = 1;
        }
        print $ofh join("\t", @x) . "\n";
    }
    close $fh;

    ## novel genes
    foreach my $key (@novel_order) {
        my $span = $supp_spans{$key};
        if ($span->{in_reference}) { next; }
        print $ofh join("\t", $span->{gene_id}, $span->{chr}, $span->{lend}, $span->{rend},
                        $span->{orient}, $span->{gene_name}, $span->{gene_type}) . "\n";
        $affected_gene_keys{ ($span->{gene_name} ne ".") ? $span->{gene_name} : $span->{gene_id} } = 1;
    }
    close $ofh;

    return(\%affected_gene_keys);
}


####
sub build_fusion_annot_idx {
    my ($new_gene_spans, $affected_gene_keys_href, $out_file) = @_;

    &process_cmd("cp " . &shq(&src_file("fusion_annot_lib.idx")) . " " . &shq($out_file));

    my %idx;
    tie (%idx, 'DB_File', $out_file, O_RDWR, 0666, $DB_BTREE) or die "Error, cannot open $out_file for update: $!";

    ## mirrors build_fusion_annot_db_index.pl, restricted to the affected genes
    ## (last record wins for a gene name listed more than once, as in a full rebuild).
    open(my $fh, $new_gene_spans) or die "Error, cannot open $new_gene_spans";
    while (<$fh>) {
        chomp;
        my ($gene_id, $chr, $lend, $rend, $orient, $gene_name, $gene_type) = split(/\t/);
        if ($gene_name && $gene_name ne ".") {
            $gene_id = $gene_name;
        }
        unless ($affected_gene_keys_href->{$gene_id}) { next; }
        $idx{"$gene_id$;COORDS"} = "$chr:$lend-$rend:$orient";
        if ($gene_type) {
            unless ($idx{$gene_id}) {
                $idx{"$gene_id$;GENE_TYPE"} = $gene_type;
            }
        }
    }
    close $fh;

    untie(%idx);

    return;
}


####
sub build_mm2_splice_bed {
    my ($combined_gtf, $out_file) = @_;

    my $supp_bed = "$WORKDIR/all_supplements.mm2.splice.bed";
    &process_cmd(&shq("$UTILDIR/paftools.ctat.js") . " gff2bed " . &shq($combined_gtf) . " > " . &shq($supp_bed));
    unless (-s $supp_bed) {
        die "Error, paftools gff2bed produced no records from $combined_gtf";
    }

    &append_files(&src_file("ref_annot.gtf.mm2.splice.bed"), $supp_bed, $out_file);

    return;
}


####
sub build_prot_info {
    my ($combined_gtf, $tmp_outputs_href) = @_;

    ## genes with supplemental coding isoforms
    my %coding_gene_ids;
    my %chrs;
    my %supp_tids;
    open(my $fh, $combined_gtf) or die "Error, cannot open $combined_gtf";
    while (<$fh>) {
        my @x = split(/\t/);
        if ($x[8] =~ /transcript_id \"([^\"]+)\"/) {
            $supp_tids{$1} = 1;
        }
        unless ($x[2] eq "CDS") { next; }
        $x[8] =~ /gene_id \"([^\"]+)\"/ or die "Error, no gene_id in $_";
        $coding_gene_ids{$1} = 1;
        $chrs{$x[0]} = 1;
    }
    close $fh;

    ## subset gtf: all reference records for those genes, plus their supplemental records
    my $subset_gtf = "$WORKDIR/prot_info_subset.gtf";
    open(my $ofh, ">$subset_gtf") or die "Error, cannot write $subset_gtf";
    foreach my $gtf (&src_file("ref_annot.gtf"), $combined_gtf) {
        open($fh, $gtf) or die "Error, cannot open $gtf";
        while (<$fh>) {
            if (/^\#/) { next; }
            if (/gene_id \"([^\"]+)\"/ && $coding_gene_ids{$1}) {
                print $ofh $_;
            }
        }
        close $fh;
    }
    close $ofh;

    my $subset_genome = "$WORKDIR/prot_info_subset.genome.fa";
    &process_cmd("samtools faidx " . &shq("$genome_lib_dir/ref_genome.fa") . " " . join(" ", map { &shq($_) } sort keys %chrs)
                 . " > " . &shq($subset_genome));

    my $subset_prefix = "$WORKDIR/prot_info_subset";
    &process_cmd(&shq("$UTILDIR/build_prot_info_db.pl") . " --gtf " . &shq($subset_gtf) . " --genome_fa " . &shq($subset_genome)
                 . " --out_prefix " . &shq($subset_prefix) . " 2> " . &shq("$subset_prefix.log"));

    ## update a copy of the original db with the rebuilt records for these genes
    my $out_dbm = $tmp_outputs_href->{"ref_annot.prot_info.dbm"};
    &process_cmd("cp " . &shq(&src_file("ref_annot.prot_info.dbm")) . " " . &shq($out_dbm));

    my (%subset_idx, %idx);
    tie (%subset_idx, 'DB_File', "$subset_prefix.prot_info.dbm", O_RDONLY, 0, $DB_BTREE) or die "Error, cannot open $subset_prefix.prot_info.dbm: $!";
    tie (%idx, 'DB_File', $out_dbm, O_RDWR, 0666, $DB_BTREE) or die "Error, cannot open $out_dbm for update: $!";
    foreach my $key (keys %subset_idx) {
        $idx{$key} = $subset_idx{$key};
    }
    untie(%idx);
    untie(%subset_idx);

    ## append the supplemental cds and protein sequences
    foreach my $type ("cds", "pep") {
        my $out_file = $tmp_outputs_href->{"ref_annot.$type"} or next;
        &process_cmd("cp " . &shq(&src_file("ref_annot.$type")) . " " . &shq($out_file));
        open($fh, "$subset_prefix.$type") or die "Error, cannot open $subset_prefix.$type";
        open($ofh, ">>$out_file") or die "Error, cannot append to $out_file";
        my $keep = 0;
        while (<$fh>) {
            if (/^>(\S+)/) {
                $keep = $supp_tids{$1};
            }
            print $ofh $_ if $keep;
        }
        close $fh;
        close $ofh;
    }

    return;
}
