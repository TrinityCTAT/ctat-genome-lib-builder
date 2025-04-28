#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/lib");
use Pipeliner;
use Cwd;
use File::Path;


## Note, ideas related to IGH/IGL super-locus creation and masking out pseudogenes derive from Daniel Nicorici and FusionCatcher (personal comm. w/ Daniel).

my $CPU = 4;
my $outTmpDir=undef;

my $max_readlength = 150;

my $output_dir = "ctat_genome_lib_build_dir";

my $annot_filter_rule = "$FindBin::Bin/AnnotFilterRuleLib/AnnotFilterRule.pm";



my $usage = <<__EOUSAGE__;

##################################################################################
#
#  Required:
#
#  --genome_fa  <string>           genome fasta file
#
#  --gtf <string>                  transcript structure annotation
#                                     Note: can restrict to coding genes and lncRNAs
#  --dfam_db <string>              DNA transposable element database (Dfam.hmm), required for repeat masking. (** highly recommended **)
#                                  (use organism-specific library if possible. Indicating 'human' or 'mouse' will automatically pull the resources from dfam directly for convenience.
#
#  Misc options:
#
#  --output_dir <string>           output directory (default: $output_dir)
#
#  --fusion_annot_lib <string>     fusion annotation library (key/val pairs, tab-delimited)
#
#  --annot_filter_rule <string>    target AnnotFilterRule.pm (default: $annot_filter_rule)
#
#  --pfam_db <string>              /path/to/Pfam-A.hmm  
#                                  (get it from here: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz)
#                                  Note, if keyword 'current' is used, this is retrieved at runtime.
#
#  --max_readlength <int>          max length for an individual RNA-Seq read (ie. default: $max_readlength)
#
#  --CPU <int>                     number of threads (defalt: $CPU)
#
#  --outTmpDir	<string>	   passed to STAR (very useful if local disks are faster than network disks)
#
#  --human_gencode_filter      customized prep operations for human/gencode genome and annotation data.
#
##################################################################################


__EOUSAGE__

    ;



my $help_flag;
my $genome_fa_file;
my $gtf_file;
my $fusion_annot_lib;

my $gmap_build_flag = 0;
my $pfam_db = "";
my $dfam_db = "";

my $SKIP_STAR_FLAG = 0;
my $STAR_ONLY_FLAG = 0;

my $HUMAN_GENCODE_FILTER = 0;

&GetOptions ( 'h' => \$help_flag,

              # required for STAR-Fusion, FusionInspector, GMAP-fusion
              'genome_fa=s' => \$genome_fa_file,
              'gtf=s' => \$gtf_file,        
              
              # required for star
              'max_readlength=i' => \$max_readlength,
              
              # optional
              'output_dir=s' => \$output_dir,
              'CPU=i' => \$CPU,
              'outTmpDir=s' => \$outTmpDir,

              'human_gencode_filter' => \$HUMAN_GENCODE_FILTER,
              
              'fusion_annot_lib=s' => \$fusion_annot_lib,
   
              'pfam_db=s' => \$pfam_db,
              'dfam_db=s' => \$dfam_db,
   
              'annot_filter_rule=s' => \$annot_filter_rule,
              
              # misc opts
              'skip_star' => \$SKIP_STAR_FLAG,
              'STAR_ONLY' => \$STAR_ONLY_FLAG,
              
 );


if ($help_flag) {
    die $usage;
}

unless ($genome_fa_file && $gtf_file && $max_readlength && $dfam_db) {
    die $usage;
}

if (@ARGV) {
    die "Error, didn't understand parameters: @ARGV";
}
        

my @required_tools = ("STAR", "makeblastdb", "blastn", "minimap2");


# for dfam
push (@required_tools, "dfamscan.pl", "nhmmscan");

if ($pfam_db) {
    push (@required_tools, "hmmscan");
}

if ($STAR_ONLY_FLAG) {
    @required_tools = ("STAR");
}

&check_for_required_tools(@required_tools);


$genome_fa_file = &Pipeliner::ensure_full_path($genome_fa_file) if $genome_fa_file;
$gtf_file = &Pipeliner::ensure_full_path($gtf_file) if $gtf_file;
$output_dir = &Pipeliner::ensure_full_path($output_dir) if $output_dir;
$fusion_annot_lib = &Pipeliner::ensure_full_path($fusion_annot_lib) if $fusion_annot_lib;

my $UTILDIR = $FindBin::Bin . "/util";

unless ($output_dir) {
    $output_dir = cwd();
}




main: {

    my $output_dir_checkpoints_dir = "$output_dir/__chkpts";
    unless (-d $output_dir_checkpoints_dir) {
        mkpath($output_dir_checkpoints_dir) or die "Error, cannot mkpath $output_dir_checkpoints_dir";
    }
    
    #################
    # Prep the genome

    unless (-d $output_dir) {
        mkpath($output_dir) or die "Error, cannot mkpath $output_dir";
    }
    
    my $local_checkpoints_dir = "__loc_chkpts";
    unless (-d $local_checkpoints_dir) {
        mkpath($local_checkpoints_dir) or die "Error, cannot mkpath $local_checkpoints_dir";
    }
    
    my $pipeliner = new Pipeliner(-verbose => 2); ## going to need precise control over the checkpoints dir.
    
    if ($dfam_db eq "human" || $dfam_db eq "mouse") {
        $dfam_db = &prep_dfam_db($dfam_db);
    }
    if ($pfam_db eq "current") {
        $pfam_db = &prep_pfam_db();
    }
    
    my $cmd = "cp $genome_fa_file $output_dir/ref_genome.fa";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/ref_genome.fa.ok"));
    
    $cmd = "samtools faidx $output_dir/ref_genome.fa";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/ref_genome_fai.ok"));


    #################
    ## make blastable
    $cmd = "makeblastdb -in $output_dir/ref_genome.fa -dbtype nucl";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/makeblastdb.ok"));
    
    my $genome_fa_file_for_blast = "$output_dir/ref_genome.fa";  ## reset so used below.
    
    my $genome_fa_for_STAR_index = $genome_fa_file;


    my $is_grch38 = ($genome_fa_file =~ /GRCh38/i) ? 1 : 0;
    
    if ($HUMAN_GENCODE_FILTER) {
        
        my $new_gtf_file = "$gtf_file.revised.gtf";
        my $cmd = "$UTILDIR/revise_gencode_annotations.pl --gencode_gtf $gtf_file --out_gtf $new_gtf_file";
        $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/revised_gencode_annots.ok"));

        my $customized_gtf = "$gtf_file.revised.custom.gtf";
        $cmd = "$UTILDIR/add_gtf_custom_targets.pl --in_gtf $new_gtf_file --out_gtf $customized_gtf";
        $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/customized_gtf.ok"));
        
        
        $gtf_file = $customized_gtf;
        
        my $masked_genome_fa_file = "$genome_fa_file.pseudo_masked.fa";
        $cmd = "$UTILDIR/mask_confounding_features_from_genome.pl --gencode_gtf $gtf_file --genome_fa $genome_fa_file_for_blast --out_masked $masked_genome_fa_file --is_grch38 $is_grch38 ";
        $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/pseudo_mask_genome.ok"));
        
        $genome_fa_for_STAR_index = $masked_genome_fa_file;
        
    }
    
        
    ###############################
    ## and copy the annotation file
    
    $cmd = "cp $gtf_file $output_dir/ref_annot.gtf";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/ref_annot.gtf.ok"));
    
    # extract exon records    
    $cmd = "bash -c \" set -euxo pipefail; $UTILDIR/gtf_to_exon_gene_records.pl $output_dir/ref_annot.gtf  | sort -k 1,1 -k4,4g -k5,5g | uniq  > $output_dir/ref_annot.gtf.mini.sortu \" ";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/ref_annot.gtf.mini.sortu.ok"));
    

    ##
    # build minimap2 index
    ## 
    
    $cmd = "minimap2 -d $output_dir/ref_genome.fa.mm2 $output_dir/ref_genome.fa";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/mm2_genome_idx.ok"));

    $cmd = "$UTILDIR/paftools.ctat.js gff2bed $output_dir/ref_annot.gtf > $output_dir/ref_annot.gtf.mm2.splice.bed";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/mm2.splice_bed.ok"));
    

    ##
    # build star index
    ##

    unless ($SKIP_STAR_FLAG) {
    
        my $star_index = "$output_dir/ref_genome.fa.star.idx";
        unless (-d $star_index) {
            mkpath $star_index or die "Error, cannot mkdir $star_index";
        }
        
        my $maybe_tmpdir= defined($outTmpDir)? " --outTmpDir $outTmpDir " : "";
        
        $cmd = "STAR --runThreadN $CPU --runMode genomeGenerate --genomeDir $star_index $maybe_tmpdir "
            . " --genomeFastaFiles $genome_fa_for_STAR_index " ## using the pseudogene and paralog-masked genome here.
            . " --limitGenomeGenerateRAM 40419136213 "
            . " --limitSjdbInsertNsj 10000000 "
            #. " --genomeChrBinNbits 16 " # needed for >4k contigs w/ FI, default is 18 now, so commenting it out, Oct 2023 bhaas.
            . " --sjdbGTFfile $gtf_file "
            . " --sjdbOverhang $max_readlength ";
        
        $pipeliner->add_commands(new Command($cmd, "$star_index/build.ok"));
    
    }
    
    if ($STAR_ONLY_FLAG) {
        $pipeliner->run();
        print STDERR "** --STAR_ONLY set, stopping now.\n";
        exit(0);
    }

    
    $cmd = "$UTILDIR/gtf_to_gene_spans.pl $output_dir/ref_annot.gtf > $output_dir/ref_annot.gtf.gene_spans";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/ref_annot.gtf.gene_spans.ok") );
    
    
    #############################
    ## Begin the BLAST compendium

    # CDS and ncRNA blastn for quick homology identification
    
    my $ref_annot_cdsplus_fa = "ref_annot.cdsplus.fa";
    
    $cmd = "$UTILDIR/gtf_file_to_feature_seqs.pl --gtf_file $output_dir/ref_annot.gtf --genome_fa $output_dir/ref_genome.fa --seqType CDSplus > $ref_annot_cdsplus_fa";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/ref_annot.cdsplus.fa.ok"));
    
    if ($dfam_db) {
        # repeat masking of transposable elements.
        my $dfam_masked_cdsplus = "ref_annot.cdsplus.dfam_masked.fa";
        
        $cmd = "$UTILDIR/dfam_repeat_masker.pl --dfam_hmm $dfam_db --target_fa $ref_annot_cdsplus_fa --out_masked $dfam_masked_cdsplus --CPU $CPU";
        $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/$dfam_masked_cdsplus.ok"));
        
        $ref_annot_cdsplus_fa = $dfam_masked_cdsplus;
    }
    
    $cmd = "cp $ref_annot_cdsplus_fa $output_dir/ref_annot.cdsplus.fa";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/$ref_annot_cdsplus_fa.cp.ok"));

    $cmd = "$UTILDIR/index_cdna_seqs.pl $output_dir/ref_annot.cdsplus.fa"; 
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/$ref_annot_cdsplus_fa.idx.ok"));
        

    ## blasts for homology detection (seq-similar pair definitions)

    $cmd = "makeblastdb -in $ref_annot_cdsplus_fa -dbtype nucl";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/$ref_annot_cdsplus_fa.blidx.ok"));

    $cmd = "blastn -query $ref_annot_cdsplus_fa -db $ref_annot_cdsplus_fa -max_target_seqs 10000 -outfmt 6 -evalue 1e-10 -num_threads $CPU -dust no -lcase_masking > $ref_annot_cdsplus_fa.allvsall.outfmt6";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/$ref_annot_cdsplus_fa.allvsall.outfmt6.ok"));
    
    $cmd = "bash -c \" set -euxo pipefail; $UTILDIR/blast_outfmt6_replace_trans_id_w_gene_symbol.pl $ref_annot_cdsplus_fa $ref_annot_cdsplus_fa.allvsall.outfmt6 | gzip > $ref_annot_cdsplus_fa.allvsall.outfmt6.genesym.gz\" ";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/$ref_annot_cdsplus_fa.allvsall.outfmt6.genesym.gz.ok"));

    
    $cmd = "bash -c \" set -euxo pipefail; $UTILDIR/blast_select_single_per_gene_pair.pl $ref_annot_cdsplus_fa.allvsall.outfmt6.genesym.gz | gzip > $ref_annot_cdsplus_fa.allvsall.outfmt6.genesym.best.gz\" ";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/$ref_annot_cdsplus_fa.allvsall.outfmt6.genesym.best.gz.ok"));


    $cmd = "bash -c \" set -euxo pipefail; $UTILDIR/filter_overlapping_blast_hits.pl $ref_annot_cdsplus_fa.allvsall.outfmt6.genesym.best.gz $output_dir/ref_annot.gtf.gene_spans | gzip > $ref_annot_cdsplus_fa.allvsall.outfmt6.genesym.best.overlaps_filt.gz\" ";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/$ref_annot_cdsplus_fa.allvsall.outfmt6.genesym.best.overlaps_filt.ok"));

    $cmd = "cp $ref_annot_cdsplus_fa.allvsall.outfmt6.genesym.best.overlaps_filt.gz $output_dir/blast_pairs.dat.gz";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/cp_gene_blast_pairs.ok"));

    # index the blast hits:
    $cmd = "$UTILDIR/index_blast_pairs.pl $output_dir/blast_pairs.idx $output_dir/blast_pairs.dat.gz";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/blast_pairs.idx.ok"));
    
    ##################################
    # blast across full cDNA sequences (for read filtering)

    # extract the cDNA sequences
    my $ref_annot_cdna_fa = "ref_annot.cdna.fa";
    $cmd = "$UTILDIR/gtf_file_to_feature_seqs.pl --gtf_file $gtf_file --genome_fa $genome_fa_file --seqType cDNA > $ref_annot_cdna_fa";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/ref_annot.cdna.fa.ok"));
    
    $cmd = "cp $ref_annot_cdna_fa $output_dir/$ref_annot_cdna_fa";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/cp_ref_annot_cdna.ok"));
    
    $cmd = "$UTILDIR/index_cdna_seqs.pl $output_dir/$ref_annot_cdna_fa";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/index_ref_annot_cdna.ok"));
    

    $cmd = "makeblastdb -in $ref_annot_cdna_fa -dbtype nucl";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/$ref_annot_cdna_fa.blastidx.ok"));
    
    $cmd = "blastn -query $ref_annot_cdna_fa -db $ref_annot_cdna_fa -max_target_seqs 10000 -outfmt 6 -evalue 1e-10 -num_threads $CPU -dust no > $ref_annot_cdna_fa.allvsall.outfmt6";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/$ref_annot_cdna_fa.blastn.allvsall.outfmt6.ok"));
    
    $cmd = "$UTILDIR/isoform_blast_gene_chr_conversion.pl --blast_outfmt6 $ref_annot_cdna_fa.allvsall.outfmt6 --gtf $gtf_file > $ref_annot_cdna_fa.allvsall.outfmt6.toGenes";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/$ref_annot_cdna_fa.allvsall.blastn.outfmt6.toGenes.ok"));

    $cmd = "sort -k2,2 -k7,7 $ref_annot_cdna_fa.allvsall.outfmt6.toGenes > $ref_annot_cdna_fa.allvsall.outfmt6.toGenes.sorted";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/$ref_annot_cdna_fa.allvsall.blastn.outfmt6.toGenes.sorted.ok"));

    $cmd = "gzip -f $ref_annot_cdna_fa.allvsall.outfmt6.toGenes.sorted";
    $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/$ref_annot_cdna_fa.allvsall.blastn.outfmt6.toGenes.sorted.gzip.ok"));
    
    $cmd = "cp $ref_annot_cdna_fa.allvsall.outfmt6.toGenes.sorted.gz $output_dir/trans.blast.dat.gz";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/trans.blast.dat.cp.ok"));

    $cmd = "$UTILDIR/build_chr_gene_alignment_index.pl --blast_genes $output_dir/trans.blast.dat.gz  --out_prefix $output_dir/trans.blast.align_coords";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/trans.blast.dat.index.ok"));
    
    
    ####################################
    ## integrate protein structure info
    
    $cmd = "$UTILDIR/build_prot_info_db.pl --gtf $gtf_file --genome_fa $genome_fa_file --out_prefix $output_dir/ref_annot";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/_prot_info_db.ok"));
    
    
    ######################
    # build the fusion annotation database
    # 

    my $build_time = time();  ## so we always rerun the annotation build step, no skipping here...
    ## TODO: include --force to rebuild
    # copy over the AnnotFilterRule:
    $cmd = "cp $annot_filter_rule $output_dir/.";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/annotfiltrule_cp.ok"));
    
    if ($fusion_annot_lib) {
        unless ($fusion_annot_lib =~ /\.gz$/) {
            die "Error, fusion annot lib: $fusion_annot_lib doesn't appear to be gzipped";
        }
        $cmd = "cp $fusion_annot_lib $output_dir/fusion_annot_lib.gz";
        $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/fusion_annot_lib.cp.ok"));
        $fusion_annot_lib = "$output_dir/fusion_annot_lib.gz";
    }
    
    $cmd = "$UTILDIR/build_fusion_annot_db_index.pl --gene_spans $output_dir/ref_annot.gtf.gene_spans --out_db_file $output_dir/fusion_annot_lib.idx";
    if ($fusion_annot_lib) {
        $cmd .= " --key_pairs $fusion_annot_lib";
    }
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/_fusion_annot_lib.idx.ok"));
    

    ############################################################
    ## Sections below are optional depending on optional params
    ############################################################
    
    
    if ($pfam_db) {

        # extract the protein sequences:
        my $cmd = "$UTILDIR/gtf_file_to_feature_seqs.pl --gtf_file $gtf_file --genome_fa $genome_fa_file --seqType prot > ref_annot.pep.fa"; 
        $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/make_pep_file.ok"));

        # run pfam
        $cmd = "hmmscan --cpu 4 --domtblout PFAM.domtblout.dat $pfam_db ref_annot.pep.fa";
        $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/hmmscan.ok"));

        # gzip pfam results
        $cmd= "gzip PFAM.domtblout.dat";
        $pipeliner->add_commands(new Command($cmd, "$local_checkpoints_dir/gzip_pfam.ok"));

        $cmd = "cp PFAM.domtblout.dat.gz $output_dir/PFAM.domtblout.dat.gz";
        $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/cp_pfam_dat.ok"));

        # index the pfam hits:
        $cmd = "$UTILDIR/index_pfam_domain_info.pl --pfam_domains $output_dir/PFAM.domtblout.dat.gz --genome_lib_dir $output_dir";
        $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/index_pfam_hits.ok"));
        
    }


    ## perform final validation:
    $cmd = "$UTILDIR/validate_ctat_genome_lib.pl $output_dir";
    $pipeliner->add_commands(new Command($cmd, "$output_dir_checkpoints_dir/validate_ctat_genome_lib.ok"));
    
    
    $pipeliner->run();

    exit(0);
}


####
sub check_for_required_tools {
    my (@tools) = @_;

    my $missing_flag = 0;

    foreach my $tool (@tools) {
        my $path = `which $tool`;
        if ($path =~ /\w/) {
            print STDERR "-found $tool at $path\n";
        }
        else {
            print STDERR "- *** ERROR, cannot locate required tool in PATH setting: $tool ***\n";
            $missing_flag = 1;
        }
    }

    if ($missing_flag) {
        
        die "Error, missing at least one required tool. See error messages and perform required software installations before running.";
    }

}

####
sub prep_dfam_db {
    my ($dfam_db_type) = @_;

    my $pipeliner = new Pipeliner(-verbose => 2,
                                  -checkpoint_dir => "_dfam_db_prep_chckpts");

    my $db_base_url = "http://dfam.org/releases/Dfam_3.1/infrastructure/dfamscan";
    
    my %db_filenames = ( 'human' => "homo_sapiens_dfam.hmm",
                         'mouse' => "mus_musculus_dfam.hmm" );

    my $db_filename = $db_filenames{$dfam_db_type} or confess "Error, cannot find url for type: $dfam_db_type"; 
    my $url = $db_base_url . "/" . $db_filename;
    
    ## pull down database and indexes.
    foreach my $index_type ("", "h3f", "h3i", "h3m", "h3p") {
        my $wget_url = "$url";
        if ($index_type) {
            $wget_url .= ".$index_type";
        }
        my $cmd = "wget $wget_url";
        $pipeliner->add_commands(new Command($cmd, "dfam.$index_type.ok"));
    }

    $pipeliner->run();
    
    return($db_filename);
}

####
sub prep_pfam_db {
    
    my $pipeliner = new Pipeliner(-verbose => 2,
                                  -checkpoint_dir => "_pfam_db_prep_chckpts");
    
    my $url = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz";
    
    $pipeliner->add_commands(new Command("wget $url", "pfam_wget.ok"));

    $pipeliner->add_commands(new Command("gunzip Pfam-A.hmm.gz", "pfam_gunzip.ok"));

    $pipeliner->add_commands(new Command("hmmpress Pfam-A.hmm", "pfam_hmmpress.ok"));
    

    $pipeliner->run();
    
    return("Pfam-A.hmm");
}
