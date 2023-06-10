#!/bin/bash

##e.g. ./prep_graph_reps.sh N250_R3_popsize.coal 10

PREFIX=${1?Error: no input example file (e.g. N250_R1_popsize.coal)}
NREPS=${2?Error: number of replicates not specified}
for i in $(eval echo {1..$2} ) 
do
cat $(echo ${1} | \
	sed -e "s/_R[^_]*_/_R${i}_/") | \
	tail -n +2 | \
	sed 's/^0 0 //' | \
	sed 's/ $//g' | \
	sed 's/ /,/g' \
	> $(echo ${1} | \
	sed -e "s/_R[^_]*_.*/_R${i}.recoal/")
done
