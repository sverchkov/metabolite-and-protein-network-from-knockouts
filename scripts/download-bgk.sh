#!/bin/bash

dir=raw-bgk

mkdir -p ${dir}
echo "Created by scripts/download-bgk.sh" > ${dir}/README

curl https://stringdb-static.org/download/protein.links.v10.5.txt.gz -o ${dir}/protein.links.v10.5.txt.gz