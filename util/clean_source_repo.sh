#!/bin/bash

set -ex

if [ $# -ne 1 ]; then
    echo -e "\n\tusage: outdirname\n\n"
    exit 1
fi

dirname=$1

if [ ! -e $dirname ];then
    mkdir $dirname
fi

set +e

mv Pfam-A.* \
    homo_sapiens_dfam.* \
    *cmds \
    ref_annot.cdsplus.dfam_masked.fa.n* ref_annot.cdsplus.dfam_masked.fa.allvsall.outfmt6.genesym.gz ref_annot.cdsplus.dfam_masked.fa.allvsall.outfmt6.genesym.best.gz \
    *outfmt6 ref_annot.cdna.fa.n* \
    *outfmt6.toGenes \
    *ok \
    *log \
    __gencode_refinement_chkpts/gencode.*.annotation.gtf* \
    $dirname/


