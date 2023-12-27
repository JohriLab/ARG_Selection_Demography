#!/bin/bash

# ./build_mask.sh mask_BGS_nopos_constant_R9 1040001 exon_UTR_pos.txt 

VCF=${1?Error: no input VCF (e.g. mask_BGS_nopos_constant_R9)}
END=${2?Error: end of genomic interval not specified (e.g. 1040001 assumes starts at position 1 and ends +1)}
REGIONS=${3?Error: regions to be masked not specified (e.g. exon_UTR_pos.txt see hashed comments below)} 

#grep -ae "initializeGenomicElement(g4" mask_BGS_nopos_constant.o | sort | uniq > UTR.temp
#grep -ae "initializeGenomicElement(g5" mask_BGS_nopos_constant.o | sort | uniq >> UTR.temp
#sed 's/.*g...//' UTR.temp | sed 's/,//' | sed 's/).*//' | sort -nk1 > UTR.txt
## Already 1-base for exons (2col_positions.txt)
#./print_nums.sh 2col_positions.txt > exon_positions.txt
#./print_nums.sh UTR.txt | awk '$1 = $1 + 1 {print $0}' > UTR_positions.txt

grep MT= $1.vcf | grep -v MT=5 | awk '{print $2}' > mno5_$1.sites
python -c "for i in range(1, $2): print(i)" | cat - $3 mno5_$1.sites | sort | uniq -c | sort -k2 -n > $1_counts.txt
echo ">1:1-${2} mask fasta for Relate" > $1.mask
sed 's/.*1 .*/P/' $1_counts.txt | sed 's/.* .*/R/' | tr -d '\n' >> $1.mask
