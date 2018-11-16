#!/bin/bash

dir=raw-bgk

mkdir -p ${dir}
cd ${dir}

echo "Created by scripts/download-bgk.sh" > README

# protein-proteins links with association scores. Not informative enough for our purposes.
#curl https://stringdb-static.org/download/protein.links.v10.5.txt.gz -o protein.links.v10.5.txt.gz
#gunzip protein.links.v10.5.txt.gz

# download mappings from uniprot to string-db ids
curl https://string-db.org/mapping_files/uniprot_mappings/full_uniprot_2_string.04_2015.tsv.gz -o full_uniprot_2_string.04_2015.tsv.gz
gunzip full_uniprot_2_string.04_2015.tsv.gz

# download protein-protein interactions with interaction type annotation
curl https://stringdb-static.org/download/protein.actions.v10.5.txt.gz -o protein.actions.v10.5.txt.gz
gunzip protein.actions.v10.5.txt.gz

# filter protein-proteing interactions down to human genes only
#  sed command explanation:
#   -n - don't print to stdout unless p command is used
#   1p - print the first line (table header)
#   ; - separates commands
#   /626523\./ - apply command to lines containing "626523." (the . is escaped since this is regex)
#   p - print
# NOTE: 626523 designates the species (human) in string-db
sed -n '1p;/9606\./p' protein.actions.v10.5.txt > human.protein.actions.txt