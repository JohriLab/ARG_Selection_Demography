
##Alpha calculation

for y in {low,mid,mhi,high}
do
for x in {low,high}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "huber_BGS_${x}S_${y}G_${q}"
grep -ae p1 /nas/longleaf/home/jimarsh/LOGFILES/re_production_huber_BGS_${x}S_${y}G_${q}.o | awk '{if ($NF >= 20000) print}' | grep "m." | sed 's/ [^ ]* [^ ]* p1.*//' | sed 's/.*m/m/' | grep -wf exon_positions.txt - | sed 's/ .*//' | sort | uniq -c | awk '{if ($2 == "m1") $1 = $1*.4; print}' > mcount_re_production_huber_BGS_${x}S_${y}G_${q}.txt
grep m6 mcount_re_production_huber_BGS_${x}S_${y}G_${q}.txt | sed "s/ //g" | sed "s/m6//g" | xargs -I [] awk '{sum+=$1;} END {print []/sum}' mcount_re_production_huber_BGS_${x}S_${y}G_${q}.txt
done
done
done

##Divergence in neutral intergenic/intronic regions

"Theta pi for neutral interg/intron across whole genome = 0.000647246"
0.000647246*138600000 = 129819787

for y in {low,mid,mhi,high}
do
for x in {low,high}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "huber_BGS_${x}S_${y}G_${q}"
grep -ae m5 /nas/longleaf/home/jimarsh/LOGFILES/re_production_huber_BGS_${x}S_${y}G_${q}.o | awk '{if ($NF >= 20000) print}' | grep p1 | wc -l | xargs -I [] echo [] "/1298197870" | bc -l
done
done
done

for x in {huber,johri}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "${x}_${q}"
grep -ae m5 /nas/longleaf/home/jimarsh/LOGFILES/re_production_${x}_BGS_nopos_${q}.o | awk '{if ($NF >= 20000) print}' | grep p1 | wc -l | xargs -I [] echo [] "/1298197870" | bc -l
done
done
done

for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "neutral_${q}"
grep -ae m5 /nas/longleaf/home/jimarsh/LOGFILES/re_production_neutral_${q}.o | awk '{if ($NF >= 20000) print}' | grep p1 | wc -l | xargs -I [] echo [] "/1298197870" | bc -l
done
done

#Then need to divide each demographic scenario by nGenerations after 6N

#e.g. "1879832"
#div = (1879832/10)/129819787

##Divergence in exonic regions


for y in {low,mid,mhi,high}
do
for x in {low,high}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "huber_BGS_${x}S_${y}G_${q}"

grep -ae p1 /nas/longleaf/home/jimarsh/LOGFILES/re_production_huber_BGS_${x}S_${y}G_${q}.o | awk '{if ($NF >= 20000) print}' | sed 's/.*m. //' | sed 's/ .*//' > subs_re_production_huber_BGS_${x}S_${y}G_${q}.txt

grep -wf exon_positions.txt subs_re_production_huber_BGS_${x}S_${y}G_${q}.txt | wc -l | xargs -I [] echo "[]/(11*797*311*10)" | bc -l
done
done
done

for x in {huber,johri}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "${x}_${q}"
grep -ae p1 /nas/longleaf/home/jimarsh/LOGFILES/re_production_${x}_BGS_nopos_${q}.o | awk '{if ($NF >= 20000) print}' | sed 's/.*m. //' | sed 's/ .*//' > subs_re_production_${x}_BGS_nopos_${q}.txt

grep -wf exon_positions.txt subs_re_production_${x}_BGS_nopos_${q}.txt | wc -l | xargs -I [] echo "[]/(11*797*311*10)" | bc -l
done
done


for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "neutral_${q}"

grep -ae p1 /nas/longleaf/home/jimarsh/LOGFILES/re_production_neutral_${q}.o | awk '{if ($NF >= 20000) print}' | sed 's/.*m. //' | sed 's/ .*//' > subs_re_production_neutral_${q}.txt

grep -wf exon_positions.txt subs_re_production_neutral_${q}.txt | wc -l | xargs -I [] echo "[]/(11*797*311*10)" | bc -l

done


##Pi in neutral intergenic/intronic regions

for y in {low,mid,mhi,high}
do
for x in {low,high}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "huber_BGS_${x}S_${y}G_${q}"
grep -ae "Rep." /nas/longleaf/home/jimarsh/LOGFILES/re_production_huber_BGS_${x}S_${y}G_${q}.o | sed 's/.*= //' | sed 's/"//' | awk '{ tot+=$1;cnt++ } END { print tot/cnt }'
done
done
done

for x in {huber,johri}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "${x}_nopos_${q}"
grep -ae "Rep." /nas/longleaf/home/jimarsh/LOGFILES/re_production_${x}_BGS_nopos_${q}.o | sed 's/.*= //' | sed 's/"//' | awk '{ tot+=$1;cnt++ } END { print tot/cnt }'
done
done

for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "neutral_${q}"
grep -ae "Rep." /nas/longleaf/home/jimarsh/LOGFILES/re_production_neutral_${q}.o | sed 's/.*= //' | sed 's/"//' | awk '{ tot+=$1;cnt++ } END { print tot/cnt }'
done

##Pi in exonic regions
for y in {low,mid,mhi,high}
do
for x in {low,high}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo ${q}
for i in {1..10}
do
bgzip ../re_production_huber_BGS_${x}S_${y}G_${q}_R${i}.vcf
bcftools index ../re_production_huber_BGS_${x}S_${y}G_${q}_R${i}.vcf.gz
bedtools intersect -a ../re_production_huber_BGS_${x}S_${y}G_${q}_R${i}.vcf.gz -b ../../slim_exons.bed -header > g3_${i}.vcf
vcftools --vcf g3_${i}.vcf --site-pi --out g3_${i}
awk '{ sum += $NF } END { print sum / (11*797*311)}' g3_${i}.sites.pi >> re_production_huber_BGS_${x}S_${y}G_${q}.exon.pi
done
done
done
done

for z in {johri,huber}
do
for q in {expansion,contract,bottleneck,repbneck}
do
for i in {1..10}
do
bgzip ../re_production_${z}_BGS_nopos_${q}_R${i}.vcf
bcftools index ../re_production_${z}_BGS_nopos_${q}_R${i}.vcf.gz
bedtools intersect -a ../re_production_${z}_BGS_nopos_${q}_R${i}.vcf.gz -b ../../slim_exons.bed -header > g3_${i}.vcf
vcftools --vcf g3_${i}.vcf --site-pi --out g3_${i}
awk '{ sum += $NF } END { print sum / (11*797*311)}' g3_${i}.sites.pi >> re_production_${z}_BGS_nopos_${q}.exon.pi
done
done
done

for q in {constant,expansion,contract,bottleneck,repbneck}
do
for i in {1..10} 
do
bgzip ../re_production_neutral_${q}_R${i}.vcf
bcftools index ../re_production_neutral_${q}_R${i}.vcf.gz
bedtools intersect -a ../re_production_neutral_${q}_R${i}.vcf.gz -b ../../slim_exons.bed -header > g3_${i}.vcf
vcftools --vcf g3_${i}.vcf --site-pi --out g3_${i}
awk '{ sum += $NF } END { print sum / (11*797*311)}' g3_${i}.sites.pi >> re_production_neutral_${q}.exon.pi
done
done

lh | grep re_production.*pi | sed 's/.* //' | while read line; do echo ${line}; awk '{ tot+=$1;cnt++ } END { print tot/cnt }' ${line}; done
