#!/bin/bash

dir=raw-bgk

mkdir -p ${dir}
cd ${dir}

echo "Created by scripts/download-bgk.sh" > README

curl https://stringdb-static.org/download/protein.links.v10.5.txt.gz -o protein.links.v10.5.txt.gz
gunzip protein.links.v10.5.txt.gz
