#!/bin/bash

chr=${1?Error: no chromosome name specified}

grep chr${1} muts.bed > chr${1}.muts

echo -e "start\tend\trate" > chr${1}_mutation.map

rm temp.map 

bedtools getfasta -fi GRCh37_latest_genomic.fna -bed chr${1}.muts > chr${1}_bins.fa

sed 's/a/A/g' chr${1}_bins.fa | sed 's/c/C/g' | sed 's/t/T/g' | sed 's/g/G/g' | sed 's/Chr/chr/g' > allCAPS_chr${1}.fa
echo ">end" >> allCAPS_chr${1}.fa

grep '>' allCAPS_chr${1}.fa | sed 's/.*-//g' | while read line; do 

grep -m 1 -A 10000000 ${line} allCAPS_chr${1}.fa | grep -m 2 -B 10000000 '>' | grep -v '>' > temp.txt

for i in {A,C,G,T}; do grep -o ${i} temp.txt | uniq -c; done > temp_counts.txt

for i in {A,C,G,T}; do grep ${i} temp_counts.txt | awk '{print $1}' |  xargs -I [] awk "BEGIN {print []/$(awk '{sum += $1} END {print sum}' temp_counts.txt)}" | sed "s/^/${i},/g"; done  > temp_props.txt

grep A temp_props.txt | sed 's/.,//g' | xargs -I [] awk -v var="${line}" '$3 == var {sum = $9+$14+$16; result = sum * [] ; print result}' chr${1}.muts | sed "s/^/A,/" > temp_mprops.txt
grep T temp_props.txt | sed 's/.,//g' | xargs -I [] awk -v var="${line}" '$3 == var {sum = $8+$15+$17; result = sum * [] ; print result}' chr${1}.muts | sed "s/^/T,/" >> temp_mprops.txt
grep G temp_props.txt | sed 's/.,//g' | xargs -I [] awk -v var="${line}" '$3 == var {sum = $7+$10+$13; result = sum * [] ; print result}' chr${1}.muts | sed "s/^/G,/" >> temp_mprops.txt
grep C temp_props.txt | sed 's/.,//g' | xargs -I [] awk -v var="${line}" '$3 == var {sum = $6+$11+$12; result = sum * [] ; print result}' chr${1}.muts | sed "s/^/C,/" >> temp_mprops.txt

awk -F ',' '{sum += $2} END {print sum}' temp_mprops.txt | sed "s/^/$line\t/" >> temp.map

done

#echo -e "end\trate" > chr${1}_mutation.map
head -n-1 temp.map >> chr${1}_mutation.map
