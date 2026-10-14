#!/bin/bash

# Tests for util/supplement_ctat_genome_lib.pl and util/revert_ctat_genome_lib_supplement.pl
#
# usage: ./run_tests.sh /path/to/small/ctat_genome_lib_build_dir
#
#   The lib is not modified; tests run on copies under ./_test_out
#   (build a small lib with:  cd ../test_ctat_build && make test)

set -u

if [ $# -ne 1 ] || [ ! -d "$1" ]; then
    echo "usage: $0 /path/to/small/ctat_genome_lib_build_dir" >&2
    exit 1
fi

TESTDIR=$(cd $(dirname $0) && pwd)
UTILDIR=$TESTDIR/../../util
LIB=$(cd $1 && pwd)
OUT=$TESTDIR/_test_out

NUM_FAILED=0

ok() { echo "  ok   $1"; }
fail() { echo "  FAIL $1"; NUM_FAILED=$((NUM_FAILED+1)); }
check() { if eval "$2"; then ok "$1"; else fail "$1"; fi; }

# run a command, logging to last.log; returns its exit status
run() { "$@" > $OUT/last.log 2>&1; }

# md5 of every lib file outside supplemental/
snapshot() { (cd $1 && find . -type f ! -path './supplemental/*' -exec md5sum {} + | sort -k2); }

dump_db() {
    perl -MDB_File -e 'my %h; tie(%h, "DB_File", $ARGV[0], O_RDONLY, 0, $DB_BTREE) or die "cannot open $ARGV[0]";
                       foreach my $k (sort keys %h) { (my $p = $k) =~ s/$;/::/g; print "$p\t$h{$k}\n"; }' $1
}

rm -rf $OUT && mkdir -p $OUT
cp -a $LIB $OUT/pristine
WORK=$OUT/work

# supplement fixtures derived from the lib itself
# supp2: a coding isoform (exon + CDS records) under a new transcript_id
awk -F'\t' -v OFS='\t' '/transcript_id "ENST00000263741.7"/ && ($3=="exon" || $3=="CDS" || $3=="transcript") {
        $2="Custom"; sub(/transcript_id "ENST00000263741.7"/, "transcript_id \"SDF4_CODING_ALT.1\"", $9); print }' \
    $LIB/ref_annot.gtf > $OUT/supp2_cds.gtf
# a reference transcript_id (must be rejected)
awk -F'\t' '/transcript_id "ENST00000263741.7"/ && $3=="exon"' $LIB/ref_annot.gtf | head -1 > $OUT/supp_collide.gtf


echo "## dry run"
cp -a $OUT/pristine $WORK
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --dry_run
check "unassignable transcripts rejected (exit 1)" "[ $? -eq 1 ]"
check "both unassignable transcripts reported" "[ \$(grep -c 'could not be assigned' $OUT/last.log) -eq 2 ]"
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --dry_run --assign_by_overlap --allow_novel_genes
check "dry run with flags succeeds" "[ $? -eq 0 ]"
check "dry run leaves lib unchanged" "diff -r $OUT/pristine $WORK > /dev/null"
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --dry_run --assign_by_overlap --allow_novel_genes --dry_run_out $OUT/dry_run.gtf
check "--dry_run_out writes the normalized supplement" "[ $? -eq 0 ] && grep -q 'gene_id \"ENSG00000078808.12\"; transcript_id \"STRG.1.1\"' $OUT/dry_run.gtf && diff -r $OUT/pristine $WORK > /dev/null"


echo "## supplement 1 (isoforms of an existing gene, a novel gene)"
orig_gtf_inode=$(stat -c %i $WORK/ref_annot.gtf)
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
check "supplement succeeds" "[ $? -eq 0 ]"
check "transcripts assigned by ref_gene_id, gene_name, exon overlap" \
    "grep -q 'assigned by ref_gene_id: 1' $OUT/last.log && grep -q 'assigned by gene_name: 1' $OUT/last.log && grep -q 'assigned by exon overlap: 1' $OUT/last.log"
check "unstranded transcript skipped" "grep -q 'skipped (no strand): 1' $OUT/last.log"
for f in ref_annot.gtf ref_annot.gtf.mini.sortu ref_annot.gtf.gene_spans fusion_annot_lib.idx ref_annot.gtf.mm2.splice.bed; do
    check ".bkup identical to original: $f" "cmp -s $OUT/pristine/$f $WORK/$f.bkup"
done
for f in ref_annot.gtf ref_annot.gtf.mini.sortu ref_annot.gtf.mm2.splice.bed; do
    n=$(wc -l < $OUT/pristine/$f)
    check "original records intact at start of $f" "head -n $n $WORK/$f | cmp -s - $OUT/pristine/$f"
done
check "backup is a hard link to the original (no copy)" "[ \$(stat -c %i $WORK/ref_annot.gtf.bkup) -eq $orig_gtf_inode ]"
check "user warned that backups are hard links" "grep -q 'WARNING: the new .bkup files were made as HARD LINKS' $OUT/last.log"
check "mini.sortu gains only exons not already in the reference" \
    "[ \$(( \$(wc -l < $WORK/ref_annot.gtf.mini.sortu) - \$(wc -l < $OUT/pristine/ref_annot.gtf.mini.sortu) )) -eq 6 ]"
check "SDF4 gene span extended" "grep -qP '^ENSG00000078808.12\tchr1\t1152288\t1170200\t-\tSDF4' $WORK/ref_annot.gtf.gene_spans"
check "novel gene added to gene_spans" "grep -qP '^NOVEL1\tchr1\t2800000\t2801400\t\+\tNOVEL1\tnovel_gene' $WORK/ref_annot.gtf.gene_spans"
$UTILDIR/build_fusion_annot_db_index.pl --gene_spans $WORK/ref_annot.gtf.gene_spans --out_db_file $OUT/full_rebuild.idx \
    $( [ -s $LIB/fusion_annot_lib.gz ] && echo "--key_pairs $LIB/fusion_annot_lib.gz" ) 2> /dev/null
check "updated fusion_annot_lib.idx matches a full rebuild" "cmp -s <(dump_db $OUT/full_rebuild.idx) <(dump_db $WORK/fusion_annot_lib.idx)"
check "prot_info untouched (no CDS in supplement)" "[ ! -e $WORK/ref_annot.prot_info.dbm.bkup ] && cmp -s $OUT/pristine/ref_annot.prot_info.dbm $WORK/ref_annot.prot_info.dbm"
check "validator reports supplemented lib" "perl $UTILDIR/validate_ctat_genome_lib.pl $WORK 2>&1 | grep -q 'has been supplemented'"
snapshot $WORK > $OUT/after1.md5


echo "## rejected supplements leave the lib unchanged"
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
check "same label rejected (exit 1)" "[ $? -eq 1 ] && grep -q 'already added' $OUT/last.log"
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --label again --assign_by_overlap --allow_novel_genes
check "transcript_ids from an earlier supplement rejected" "[ $? -eq 1 ] && [ \$(grep -c 'already exists in the earlier supplement' $OUT/last.log) -eq 4 ]"
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $OUT/supp_collide.gtf
check "reference transcript_id rejected" "[ $? -eq 1 ] && grep -q 'already exists in the reference' $OUT/last.log"
check "lib unchanged" "snapshot $WORK | cmp -s - $OUT/after1.md5"
check "no lock or temp files left" "[ ! -d $WORK/supplemental/.lock ] && ! ls $WORK/*.tmp.* > /dev/null 2>&1"


echo "## supplement 2 (coding isoform)"
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $OUT/supp2_cds.gtf
check "supplement succeeds" "[ $? -eq 0 ]"
check "ref_annot.gtf == original + supplement 1 + supplement 2 (no stacking)" \
    "cat $WORK/ref_annot.gtf.bkup $WORK/supplemental/supp_001.supp1.gtf $WORK/supplemental/supp_002.supp2_cds.gtf | cmp -s - $WORK/ref_annot.gtf"
check "prot_info backed up" "cmp -s $OUT/pristine/ref_annot.prot_info.dbm $WORK/ref_annot.prot_info.dbm.bkup"
check "only SDF4 prot_info records changed, adding the new isoform" \
    "perl -MDB_File -MJSON::PP -e '
        my (%a, %b); tie(%a, \"DB_File\", \$ARGV[0], O_RDONLY, 0, \$DB_BTREE) or die; tie(%b, \"DB_File\", \$ARGV[1], O_RDONLY, 0, \$DB_BTREE) or die;
        my @changed = sort grep { !exists \$a{\$_} || \$a{\$_} ne \$b{\$_} } keys %b;
        my @lost = grep { !exists \$b{\$_} } keys %a;
        my @before = sort map { \$_->{cds_id} } @{decode_json(\$a{SDF4})};
        my @after = sort map { \$_->{cds_id} } @{decode_json(\$b{SDF4})};
        exit !(\"@changed\" eq \"ENSG00000078808.12 SDF4\" && !@lost && \"@after\" eq join(\" \", sort(@before, \"SDF4_CODING_ALT.1\")));
     ' $OUT/pristine/ref_annot.prot_info.dbm $WORK/ref_annot.prot_info.dbm"
check "new isoform appended to ref_annot.cds and .pep" \
    "grep -q '>SDF4_CODING_ALT.1' $WORK/ref_annot.cds && grep -q '>SDF4_CODING_ALT.1' $WORK/ref_annot.pep"


echo "## supplement 3 (exon overlap with a gene added by an earlier supplement)"
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp3_overlaps_novel.gtf --assign_by_overlap
check "assigned to earlier novel gene by exon overlap" \
    "[ $? -eq 0 ] && grep -q 'gene_id \"NOVEL1\"; transcript_id \"STRG.50.1\"' $WORK/supplemental/supp_003.supp3_overlaps_novel.gtf"


echo "## revert"
run $UTILDIR/revert_ctat_genome_lib_supplement.pl --genome_lib_dir $WORK
check "revert succeeds" "[ $? -eq 0 ]"
check "lib identical to original" "diff -r $OUT/pristine $WORK > /dev/null"
run $UTILDIR/revert_ctat_genome_lib_supplement.pl --genome_lib_dir $WORK
check "second revert reports nothing to revert (exit 1)" "[ $? -eq 1 ] && grep -q 'has not been supplemented' $OUT/last.log"

echo "## revert refuses a modified backup"
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
echo "# tampered" >> $WORK/ref_annot.gtf.gene_spans.bkup
snapshot $WORK > $OUT/before_bad_revert.md5
run $UTILDIR/revert_ctat_genome_lib_supplement.pl --genome_lib_dir $WORK
check "revert refused (exit 1)" "[ $? -eq 1 ] && grep -q 'differs from the one recorded' $OUT/last.log"
check "no files changed by refused revert" "snapshot $WORK | cmp -s - $OUT/before_bad_revert.md5"
rm -rf $WORK


echo "## linked copy (source lib left unchanged)"
SRC=$OUT/src
COPY=$OUT/linked
cp -a $OUT/pristine $SRC
snapshot $SRC > $OUT/src.md5
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $SRC --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes --linked_copy $COPY
check "linked copy supplement succeeds" "[ $? -eq 0 ]"
check "user warned before and after that the copy uses hard links" \
    "grep -q 'using HARD LINKS' $OUT/last.log && grep -q 'are HARD LINKS to the files in' $OUT/last.log"
check "source lib unchanged" "snapshot $SRC | cmp -s - $OUT/src.md5 && [ ! -e $SRC/supplemental ]"
check "unmodified files are hard links to the source" \
    "[ \$(stat -c %i $COPY/ref_genome.fa) -eq \$(stat -c %i $SRC/ref_genome.fa) ] && [ \$(stat -c %i $COPY/ref_annot.prot_info.dbm) -eq \$(stat -c %i $SRC/ref_annot.prot_info.dbm) ]"
check "no symlinks created" "[ \$(find $COPY -type l | wc -l) -eq \$(find $SRC -type l | wc -l) ]"
check "modified files are new; .bkup links the source original" \
    "[ \$(stat -c %i $COPY/ref_annot.gtf) -ne \$(stat -c %i $SRC/ref_annot.gtf) ] && [ \$(stat -c %i $COPY/ref_annot.gtf.bkup) -eq \$(stat -c %i $SRC/ref_annot.gtf) ]"
check "copy's ref_annot.gtf == source + supplement" "cat $SRC/ref_annot.gtf $COPY/supplemental/supp_001.supp1.gtf | cmp -s - $COPY/ref_annot.gtf"
check "manifest records source lib" "grep -q '\"source_lib\" : \"$SRC\"' $COPY/supplemental/manifest.json"
check "validator reports linked copy" "perl $UTILDIR/validate_ctat_genome_lib.pl $COPY 2>&1 | grep -q 'linked copy of \\[$SRC\\]'"
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $COPY --supp_gtf $OUT/supp2_cds.gtf
check "second supplement to the copy succeeds; source unchanged" "[ $? -eq 0 ] && snapshot $SRC | cmp -s - $OUT/src.md5"
run $UTILDIR/revert_ctat_genome_lib_supplement.pl --genome_lib_dir $COPY
check "reverted copy identical to source" "[ $? -eq 0 ] && diff -r $SRC $COPY > /dev/null"
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $SRC --supp_gtf $TESTDIR/supp1.gtf --linked_copy $COPY
check "non-empty target rejected; source and target untouched" \
    "[ $? -eq 1 ] && grep -q 'not empty' $OUT/last.log && snapshot $SRC | cmp -s - $OUT/src.md5 && diff -r $SRC $COPY > /dev/null"
rm -rf $COPY
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $SRC --supp_gtf $TESTDIR/supp1.gtf --linked_copy $COPY
check "failed supplement removes the linked copy" "[ $? -eq 1 ] && [ ! -e $COPY ] && snapshot $SRC | cmp -s - $OUT/src.md5"
rm -rf $SRC


fresh() { rm -rf "${WORK:?}" && cp -a "${OUT:?}/pristine" "$WORK"; }
supp() { run $UTILDIR/supplement_ctat_genome_lib.pl "$@"; }
SDF4_EXONS() { awk -F'\t' '/transcript_id "ENST00000263741.7"/ && $3=="exon"' $OUT/pristine/ref_annot.gtf; }


echo "## a second StringTie supplement: STRG ids are not conflated with earlier novel STRG genes"
fresh
printf 'chr1\tStringTie\texon\t2950000\t2950400\t1000\t+\t.\tgene_id "STRG.100"; transcript_id "STRG.100.1";\nchr1\tStringTie\texon\t2951000\t2951300\t1000\t+\t.\tgene_id "STRG.100"; transcript_id "STRG.100.2";\nchr1\tStringTie\texon\t2970000\t2970400\t1000\t+\t.\tgene_id "STRG.101"; transcript_id "STRG.101.1";\n' > $OUT/strg_run1.gtf
SDF4_EXONS | awk -F'\t' -v OFS='\t' '{$2="StringTie"; $9="gene_id \"STRG.200\"; transcript_id \"STRG.200.1\"; ref_gene_name \"SDF4\";"; print}' > $OUT/strg_run2.gtf
supp --genome_lib_dir $WORK --supp_gtf $OUT/strg_run1.gtf --allow_novel_genes
check "run 1 adds novel genes STRG.100, STRG.101" "[ $? -eq 0 ] && grep -qP '^STRG.100\t' $WORK/ref_annot.gtf.gene_spans && grep -qP '^STRG.101\t' $WORK/ref_annot.gtf.gene_spans"
supp --genome_lib_dir $WORK --supp_gtf $OUT/strg_run2.gtf
check "run 2 transcript assigned to SDF4 by ref_gene_name, not to an earlier STRG gene" \
    "[ $? -eq 0 ] && grep -q 'assigned by ref_gene_name: 1' $OUT/last.log && grep -q 'gene_id \"ENSG00000078808.12\"; transcript_id \"STRG.200.1\"' $WORK/supplemental/supp_002.strg_run2.gtf"
check "earlier novel gene spans unchanged" "grep -qP '^STRG.100\tchr1\t2950000\t2951300\t' $WORK/ref_annot.gtf.gene_spans"


echo "## a transcript not overlapping the gene it is assigned to"
fresh
printf 'chr1\tx\texon\t2990000\t2990400\t.\t-\t.\tgene_id "G1"; transcript_id "FAR.1"; gene_name "SDF4";\n' > $OUT/far.gtf
supp --genome_lib_dir $WORK --supp_gtf $OUT/far.gtf
check "rejected by default" "[ $? -eq 1 ] && grep -q 'does not overlap the gene it was assigned to by gene_name' $OUT/last.log && snapshot $WORK | cmp -s - <(snapshot $OUT/pristine)"
supp --genome_lib_dir $WORK --supp_gtf $OUT/far.gtf --allow_nonoverlapping
check "allowed with --allow_nonoverlapping (span extended, with a warning)" \
    "[ $? -eq 0 ] && grep -q 'the gene span will be extended' $OUT/last.log && grep -qP '^ENSG00000078808.12\tchr1\t1152288\t2990400\t' $WORK/ref_annot.gtf.gene_spans"


echo "## GTF parsing"
fresh
printf 'chr1\tx\texon\t1158624\t1158737\t.\t-\t.\tgene_id "ENSG00000078808.12"; transcript_id "NOSEMI.1"\n' > $OUT/nosemi.gtf
supp --genome_lib_dir $WORK --supp_gtf $OUT/nosemi.gtf --dry_run
check "attribute column without a final ';' is parsed" "[ $? -eq 0 ] && grep -q 'transcripts kept: 1' $OUT/last.log"
printf 'chr1\tx\texon\t1158624\t1158737\t.\t-\t.\tgene_id "G1"; transcript_id "T1";\nchr1\tx\texon\t1163848\t1164356\t.\t-\t.\tgene_id "G2"; transcript_id "T1";\n' > $OUT/conflict.gtf
supp --genome_lib_dir $WORK --supp_gtf $OUT/conflict.gtf --dry_run --allow_novel_genes
check "a transcript_id under two gene_ids is rejected" "[ $? -eq 1 ] && grep -q 'has records with different gene_ids' $OUT/last.log"
printf 'chr1\tx\texon\t1158624\n' > $OUT/bad.gtf
printf 'chr1\tx\texon\t1158624\t1158737\t.\t-\t.\tgene_id "G1";\n' >> $OUT/bad.gtf
supp --genome_lib_dir $WORK --supp_gtf $OUT/bad.gtf --dry_run
check "all structural problems reported together" "[ $? -eq 1 ] && grep -q 'does not have 9 tab-delimited fields' $OUT/last.log && grep -q 'lacks a transcript_id' $OUT/last.log"
gzip -c $OUT/nosemi.gtf > $OUT/nosemi.gtf.gz
supp --genome_lib_dir $WORK --supp_gtf $OUT/nosemi.gtf.gz --dry_run
check "gzipped input reported as such" "[ $? -eq 1 ] && grep -q 'gzip-compressed' $OUT/last.log"


echo "## reference gene without a gene_name (eg. Ensembl-style libs)"
fresh
sed -i '/ENSG00000078808.12/ s/ gene_name "SDF4";//' $WORK/ref_annot.gtf
$UTILDIR/gtf_to_exon_gene_records.pl $WORK/ref_annot.gtf | sort -k 1,1 -k4,4g -k5,5g | uniq > $WORK/ref_annot.gtf.mini.sortu
$UTILDIR/gtf_to_gene_spans.pl $WORK/ref_annot.gtf > $WORK/ref_annot.gtf.gene_spans 2> /dev/null
rm -f "${WORK:?}/fusion_annot_lib.idx"
$UTILDIR/build_fusion_annot_db_index.pl --gene_spans $WORK/ref_annot.gtf.gene_spans --out_db_file $WORK/fusion_annot_lib.idx 2> /dev/null
printf 'chr1\tx\texon\t1152400\t1153184\t.\t-\t.\tgene_id "XLOC_9"; transcript_id "OVL.1";\nchr1\tx\texon\t1158624\t1158737\t.\t-\t.\tgene_id "XLOC_9"; transcript_id "OVL.1";\n' > $OUT/ovl.gtf
supp --genome_lib_dir $WORK --supp_gtf $OUT/ovl.gtf --assign_by_overlap
check "assigned to the unnamed gene by exon overlap" "[ $? -eq 0 ] && grep -q 'assigned by exon overlap: 1' $OUT/last.log"
check "no synthetic gene_name in the supplement" "! grep -q 'gene_name' $WORK/supplemental/supp_001.ovl.gtf"
check "mini.sortu lists the gene under a single token" \
    "[ \$(grep -o 'gene_id \"[^\"]*ENSG00000078808.12[^\"]*\"' $WORK/ref_annot.gtf.mini.sortu | sort -u | wc -l) -eq 1 ]"


echo "## stale temporary files from an interrupted run"
fresh
ln $WORK/ref_annot.gtf.gene_spans $WORK/ref_annot.gtf.gene_spans.tmp.99999
mkdir -p $WORK/supplemental/.work.99999
supp --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
check "stale files removed; original intact" \
    "[ $? -eq 0 ] && [ ! -e $WORK/ref_annot.gtf.gene_spans.tmp.99999 ] && [ ! -e $WORK/supplemental/.work.99999 ] && cmp -s $OUT/pristine/ref_annot.gtf.gene_spans $WORK/ref_annot.gtf.gene_spans.bkup"


echo "## interrupted by a signal"
fresh
mkdir -p $OUT/slowbin
printf '#!/bin/bash\ntouch %s/k8.started\nsleep 60\n' $OUT > $OUT/slowbin/k8 && chmod +x $OUT/slowbin/k8
rm -f "${OUT:?}/k8.started"
PATH=$OUT/slowbin:$PATH $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes > $OUT/last.log 2>&1 &
pid=$!
for i in $(seq 1 600); do [ -e $OUT/k8.started ] && break; sleep 0.1; done
kill -TERM $pid; wait $pid; status=$?
pkill -f "$OUT/slowbin/k8" 2> /dev/null
check "terminated run exits 1 and cleans up (no lock, temp or work files); lib unchanged" \
    "[ $status -eq 1 ] && [ ! -e $WORK/supplemental/.lock ] && ! ls $WORK/*.tmp.* > /dev/null 2>&1 && ! ls -d $WORK/supplemental/.work.* > /dev/null 2>&1 && snapshot $WORK | cmp -s - <(snapshot $OUT/pristine)"
supp --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
check "next run succeeds" "[ $? -eq 0 ]"


echo "## interrupted commit is completed by rerunning"
fresh
supp --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
cp $WORK/ref_annot.gtf $OUT/expected_gtf
## simulate a kill during the renames: supplement recorded as committing, one file not yet replaced
perl -MJSON::PP -i -0pe '$m = decode_json($_); $_->{status} = "committing" for @{$m->{supplements}}; $_ = JSON::PP->new->pretty->canonical->encode($m)' $WORK/supplemental/manifest.json
rm -f "${WORK:?}/ref_annot.gtf" && ln $WORK/ref_annot.gtf.bkup $WORK/ref_annot.gtf
supp --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
check "rerun completes it" "[ $? -eq 0 ] && grep -q 'was interrupted while being committed; completing it' $OUT/last.log && cmp -s $OUT/expected_gtf $WORK/ref_annot.gtf && grep -q '\"status\" : \"complete\"' $WORK/supplemental/manifest.json"
check "no duplicate supplement recorded" "[ \$(grep -c '\"label\"' $WORK/supplemental/manifest.json) -eq 1 ]"


echo "## backup made but not recorded (interrupted run)"
fresh
mkdir $WORK/supplemental && ln $WORK/ref_annot.gtf $WORK/ref_annot.gtf.bkup
supp --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
check "identical unrecorded backup is adopted" "[ $? -eq 0 ] && grep -q 'recording backup ref_annot.gtf.bkup made by an interrupted run' $OUT/last.log"
fresh
echo "# other" > $WORK/ref_annot.gtf.gene_spans.bkup
supp --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
check "a different unrecorded backup is refused" "[ $? -eq 1 ] && grep -q 'not recorded in .* and it differs' $OUT/last.log"


echo "## interrupted revert is completed by rerunning"
fresh
supp --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
mv $WORK/ref_annot.gtf.bkup $WORK/ref_annot.gtf      # killed after restoring one file
run $UTILDIR/revert_ctat_genome_lib_supplement.pl --genome_lib_dir $WORK
check "revert after a partial revert succeeds; lib identical to original" "[ $? -eq 0 ] && diff -r $OUT/pristine $WORK > /dev/null"
supp --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
for f in ref_annot.gtf ref_annot.gtf.mini.sortu ref_annot.gtf.gene_spans fusion_annot_lib.idx ref_annot.gtf.mm2.splice.bed; do mv $WORK/$f.bkup $WORK/$f; done  # killed before removing supplemental/
run $UTILDIR/revert_ctat_genome_lib_supplement.pl --genome_lib_dir $WORK
check "revert after all files were restored succeeds; lib identical to original" "[ $? -eq 0 ] && diff -r $OUT/pristine $WORK > /dev/null"


echo "## paths with spaces; symlinked lib dir with --linked_copy"
SPACED="${OUT:?}/with space/ctat lib"
mkdir -p "${OUT:?}/with space" && cp -a $OUT/pristine "$SPACED"
supp --genome_lib_dir "$SPACED" --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
check "supplement of a lib path with spaces" "[ $? -eq 0 ]"
run $UTILDIR/revert_ctat_genome_lib_supplement.pl --genome_lib_dir "$SPACED"
check "revert of a lib path with spaces" "[ $? -eq 0 ] && diff -r $OUT/pristine \"$SPACED\" > /dev/null"
rm -rf "${OUT:?}/with space"
ln -sfn $OUT/pristine $OUT/lib_link
supp --genome_lib_dir $OUT/lib_link --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes --linked_copy $OUT/linked2
check "--linked_copy of a symlinked lib dir" "[ $? -eq 0 ] && [ \$(stat -c %i $OUT/linked2/ref_genome.fa) -eq \$(stat -c %i $OUT/pristine/ref_genome.fa) ] && [ ! -e $OUT/pristine/supplemental ]"
rm -rf "${OUT:?}/linked2" "${OUT:?}/lib_link"
supp --genome_lib_dir $OUT/pristine --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes --linked_copy $OUT/linked3 --dry_run
check "--dry_run with --linked_copy creates no copy" "[ $? -eq 0 ] && [ ! -e $OUT/linked3 ] && grep -q 'no linked copy is created' $OUT/last.log"
rm -rf "${WORK:?}"


echo "## original file lacking a final newline"
cp -a $OUT/pristine $WORK
head -c -1 $OUT/pristine/ref_annot.gtf > $WORK/ref_annot.gtf
run $UTILDIR/supplement_ctat_genome_lib.pl --genome_lib_dir $WORK --supp_gtf $TESTDIR/supp1.gtf --assign_by_overlap --allow_novel_genes
check "supplemental records start on a new line" \
    "[ $? -eq 0 ] && cat $OUT/pristine/ref_annot.gtf $WORK/supplemental/supp_001.supp1.gtf | cmp -s - $WORK/ref_annot.gtf"


echo
if [ $NUM_FAILED -eq 0 ]; then
    echo "All tests passed."
    rm -rf $OUT
    exit 0
else
    echo "$NUM_FAILED test(s) failed; see $OUT"
    exit 1
fi
