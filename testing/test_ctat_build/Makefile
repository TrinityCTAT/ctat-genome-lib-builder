

test:
	../../prep_genome_lib.pl  --genome_fa data/small.fa --gtf data/small.gtf --dfam_db /seq/RNASEQ/TOOLS/DFAM/homo_sapiens_dfam.hmm


test_docker:
	docker run --rm -it -v `pwd`:/workspace --rm trinityctat/ctatfusion bash -c 'cd /workspace && /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl --genome_fa data/small.fa --gtf data/small.gtf'


clean:
	rm -f ./Log.out
	rm -f ./ref_annot.long_readthru_introns_removed.gtf
	rm -f ./pipeliner.*.cmds
	rm -f ./ref_annot.*
	rm -rf ./ctat_genome_lib_build_dir
	rm -rf ./__loc_chkpts
