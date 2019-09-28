
test: clean
	../../prep_genome_lib.pl  --genome_fa data/small.fa --gtf data/small.gtf --dfam_db data/dfam.mini.hmm --pfam_db data/pfam.mini.hmm 


test_human_filter: clean
	../../prep_genome_lib.pl  --genome_fa data/small.fa --gtf data/small.gtf --dfam_db data/dfam.mini.hmm --pfam_db data/pfam.mini.hmm 	--human_gencode_filter

test_docker:
	docker run --rm -it -v `pwd`:/workspace --rm trinityctat/ctatfusion bash -c 'cd /workspace && /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl --genome_fa data/small.fa --gtf data/small.gtf'


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


