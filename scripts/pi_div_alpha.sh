
##Alpha calculation

for y in {low,mid,high}
do
for x in {low,high}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "huber_BGS_${x}S_${y}G_${q}"
grep -ae p1 /nas/longleaf/home/jimarsh/LOGFILES/production_huber_BGS_${x}S_${y}G_${q}.o | awk '{if ($NF >= 12000) print}' | grep "m." | sed 's/ [^ ]* [^ ]* p1.*//' | sed 's/.*m/m/' | grep -wf exon_positions.txt - | sed 's/ .*//' | sort | uniq -c | awk '{if ($2 == "m1") $1 = $1*.4; print}' > mcount_production_huber_BGS_${x}S_${y}G_${q}.txt
grep m6 mcount_production_huber_BGS_${x}S_${y}G_${q}.txt | sed "s/ //g" | sed "s/m6//g" | xargs -I [] awk '{sum+=$1;} END {print []/sum}' mcount_production_huber_BGS_${x}S_${y}G_${q}.txt
done
done
done

##Divergence in neutral intergenic/intronic regions

"Theta pi for neutral interg/intron across whole genome = 0.000647246"
0.000647246*138600000 = 129819787

for y in {low,mid,high}
do
for x in {low,high}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "huber_BGS_${x}S_${y}G_${q}"
grep -ae m5 /nas/longleaf/home/jimarsh/LOGFILES/production_huber_BGS_${x}S_${y}G_${q}.o | awk '{if ($NF >= 12000) print}' | grep p1 | wc -l | xargs -I [] echo [] "/1298197870" | bc -l
done
done
done

for x in {huber,johri}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "${x}_${q}"
grep -ae m5 /nas/longleaf/home/jimarsh/LOGFILES/production_${x}_BGS_nopos_${q}.o | awk '{if ($NF >= 12000) print}' | grep p1 | wc -l | xargs -I [] echo [] "/1298197870" | bc -l
done
done
done

for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "neutral_${q}"
grep -ae m5 /nas/longleaf/home/jimarsh/LOGFILES/production_neutral_${q}.o | awk '{if ($NF >= 12000) print}' | grep p1 | wc -l | xargs -I [] echo [] "/1298197870" | bc -l
done
done

#Then need to divide each demographic scenario by nGenerations after 6N

#e.g. "1879832"
#div = (1879832/10)/129819787

##Divergence in exonic regions


for y in {low,mid,high}
do
for x in {low,high}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "huber_BGS_${x}S_${y}G_${q}"

grep -ae p1 /nas/longleaf/home/jimarsh/LOGFILES/production_huber_BGS_${x}S_${y}G_${q}.o | awk '{if ($NF >= 12000) print}' | sed 's/.*m. //' | sed 's/ .*//' > subs_production_huber_BGS_${x}S_${y}G_${q}.txt

grep -wf exon_positions.txt subs_production_huber_BGS_${x}S_${y}G_${q}.txt | wc -l | xargs -I [] echo "[]/(11*797*311*10)" | bc -l
done
done
done

for x in {huber,johri}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "${x}_${q}"
grep -ae p1 /nas/longleaf/home/jimarsh/LOGFILES/production_${x}_BGS_nopos_${q}.o | awk '{if ($NF >= 12000) print}' | sed 's/.*m. //' | sed 's/ .*//' > subs_production_${x}_BGS_nopos_${q}.txt

grep -wf exon_positions.txt subs_production_${x}_BGS_nopos_${q}.txt | wc -l | xargs -I [] echo "[]/(11*797*311*10)" | bc -l
done
done


for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "neutral_${q}"

grep -ae p1 /nas/longleaf/home/jimarsh/LOGFILES/production_neutral_${q}.o | awk '{if ($NF >= 12000) print}' | sed 's/.*m. //' | sed 's/ .*//' > subs_production_neutral_${q}.txt

grep -wf exon_positions.txt subs_production_neutral_${q}.txt | wc -l | xargs -I [] echo "[]/(11*797*311*10)" | bc -l

done
done


##Pi in neutral intergenic/intronic regions

grep -ae "Rep." production_johri_BGS_nopos_constant.o | sed 's/.*= //' | sed 's/"//' | awk '{ tot+=$1;cnt++ } END { print tot/cnt }'

##Pi in exonic regions
for y in {low,mid,high}
do
for x in {low,high}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo ${q}
for i in {1..10}
do
bgzip ../production_huber_BGS_${x}S_${y}G_${q}_R${i}.vcf
bcftools index ../production_huber_BGS_${x}S_${y}G_${q}_R${i}.vcf.gz
bedtools intersect -a ../production_huber_BGS_${x}S_${y}G_${q}_R${i}.vcf.gz -b ../../slim_exons.bed -header > g3_${i}.vcf
vcftools --vcf g3_${i}.vcf --site-pi --out g3_${i}
awk '{ sum += $NF } END { print sum / (11*797*311)}' g3_${i}.sites.pi >> production_huber_BGS_${x}S_${y}G_${q}.exon.pi
done
done
done
done

for z in {johri,huber}
do
for i in {1..10}
do
bgzip ../production_${z}_BGS_nopos_constant_R${i}.vcf
bcftools index ../production_${z}_BGS_nopos_constant_R${i}.vcf.gz
bedtools intersect -a ../production_${z}_BGS_nopos_constant_R${i}.vcf.gz -b ../../slim_exons.bed -header > g3_${i}.vcf
vcftools --vcf g3_${i}.vcf --site-pi --out g3_${i}
awk '{ sum += $NF } END { print sum / (11*797*311)}' g3_${i}.sites.pi >> production_${z}_BGS_nopos_constant.exon.pi
done
done

