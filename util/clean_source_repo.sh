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

mv IGH* IGL* *.gtf.feature_selected *.gtf.feature_selected.noreadthrus Log.out *.nhr *.nin *nsq *.outfmt6 *.outfmt6.toGenes    $dirname/

