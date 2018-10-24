#!/bin/bash

for n in `seq 1 13080`; do
	echo JOB DOT-$n dot.condor
	echo VARS DOT-$n n=\"$n\"
	echo VARS DOT-$n result=\"dots/dots-$n.rds\"
done;
