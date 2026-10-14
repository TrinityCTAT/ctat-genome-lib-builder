
test: clean
	../../prep_genome_lib.pl  --genome_fa data/small.fa --gtf data/small.gtf --dfam_db data/dfam.mini.hmm --pfam_db data/pfam.mini.hmm 


test_human_filter: clean
	../../prep_genome_lib.pl  --genome_fa data/small.fa --gtf data/small.gtf --dfam_db data/dfam.mini.hmm --pfam_db data/pfam.mini.hmm 	--human_gencode_filter

## build with supplemental transcript models (--supp_gtf), and check they reach the lib's resources
test_supp_gtf: clean
	../../prep_genome_lib.pl  --genome_fa data/small.fa --gtf data/small.gtf --dfam_db data/dfam.mini.hmm --pfam_db data/pfam.mini.hmm --supp_gtf data/small.supp.gtf
	grep -q 'SUPPLEMENTAL.*gene_id "ENSG00000188157.9"; transcript_id "AGRN_skip5.1"; gene_name "AGRN"' ctat_genome_lib_build_dir/ref_annot.gtf
	grep -qP '^chr1\t976261\t976857\t\+' ctat_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbList.out.tab
	grep -q '^>AGRN_skip5.1' ctat_genome_lib_build_dir/ref_annot.cdna.fa
	grep -qP '^small.supp.gtf\t' ctat_genome_lib_build_dir/ref_annot.build_supplements.tsv
	test ! -e ctat_genome_lib_build_dir/supplemental
	@echo "test_supp_gtf: passed"

## the curated GRCh38 supplements must be rejected for an annotation they don't fit (this one is GRCh37)
test_curated_wrong_build:
	! ../../util/merge_supplemental_gtfs.pl --in_gtf data/small.gtf --supp_gtf ../../resources/GRCh38.ctat_curated_supplements.gtf --out_gtf curated_wrong_build.gtf
	test ! -e curated_wrong_build.gtf
	@echo "test_curated_wrong_build: passed"

test_docker:
	docker run --rm -it -v `pwd`:/workspace --rm trinityctat/starfusion bash -c 'cd /workspace && /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl --genome_fa data/small.fa --gtf data/small.gtf --dfam_db human'


clean:
	rm -f ./Log.out
	rm -f ./ref_annot.long_readthru_introns_removed.gtf
	rm -f ./pipeliner.*.cmds
	rm -f ./ref_annot.*
	rm -rf ./ctat_genome_lib_build_dir
	rm -rf ./__loc_chkpts
	rm -rf ./__dfam_ref*
	rm -rf ./__gencode_refinement_chkpts
	rm -f ./PFAM.domtblout.dat.gz
	rm -rf ./__paralogs
	rm -f ./para.mask.list
	rm -f ./data/small.gtf.*
	rm -f ./curated_wrong_build.gtf


