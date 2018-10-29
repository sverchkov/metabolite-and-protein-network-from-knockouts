#!/bin/bash

dir=raw-bgk

mkdir -p ${dir}
cd ${dir}

echo "Created by scripts/download-bgk.sh" > README

curl https://stringdb-static.org/download/protein.links.v10.5.txt.gz -o protein.links.v10.5.txt.gz
gunzip protein.links.v10.5.txt.gz

curl https://string-db.org/mapping_files/uniprot_mappings/full_uniprot_2_string.04_2015.tsv.gz -o full_uniprot_2_string.04_2015.tsv.gz
gunzip full_uniprot_2_string.04_2015.tsv.gz