##Alpha calculation

for y in {low,mid}
do
for x in {mid,high}
do
for q in {constant,expansion,contract,bottleneck,repbneck}
do
echo -e "production_BGS_${x}S_${y}G_${q}"
grep -ae p1 /nas/longleaf/home/jimarsh/LOGFILES/production_BGS_${x}S_${y}G_${q}.o | awk '{if ($NF >= 20000) print}' | grep "m." | sed 's/ [^ ]* [^ ]* p1.*//' | sed 's/.*m/m/' | grep -wf /work/users/j/i/jimarsh/drosoph_sims/exon_positions.txt - | sed 's/ .*//' | sort | uniq -c | awk '{if ($2 == "m1") $1 = $1*.4; print}' > mcount_production_BGS_${x}S_${y}G_${q}.txt
grep m6 mcount_production_BGS_${x}S_${y}G_${q}.txt | sed "s/ //g" | sed "s/m6//g" | xargs -I [] awk '{sum+=$1;} END {print []/sum}' mcount_production_BGS_${x}S_${y}G_${q}.txt
done
done
done

##Div in neutral interg/intron regions

#3732050 is length of non-selected interg/intron region * 10 replicates

for y in {low,mid}
do
	for x in {mid,high}
	do
		for q in {constant,expansion,contract,bottleneck,repbneck}
		do
echo -e "production_BGS_${x}S_${y}G_${q}"
grep -ae m5 /nas/longleaf/home/jimarsh/LOGFILES/production_BGS_${x}S_${y}G_${q}.o | awk '{if ($NF >= 100000) print}' | grep p1 | wc -l | xargs -I [] echo [] "/3732050" | bc -l
		done
	done
done

for y in {BGS_nopos,neutral,neutral_nomaps}
do
	for q in {constant,expansion,contract,bottleneck,repbneck}
	do
echo -e "production_${y}_${q}"
grep -ae m5 /nas/longleaf/home/jimarsh/LOGFILES/production_${y}_${q}.o | awk '{if ($NF >= 100000) print}' | grep p1 | wc -l | xargs -I [] echo [] "/3732050" | bc -l
	done
done

##Div in coding regions

for y in {low,mid}
do
        for x in {mid,high}
        do
                for q in {constant,expansion,contract,bottleneck,repbneck}
                do
echo -e "production_BGS_${x}S_${y}G_${q}"
grep -ae p1 /nas/longleaf/home/jimarsh/LOGFILES/production_BGS_${x}S_${y}G_${q}.o | awk '{if ($NF >= 100000) print}' | sed 's/.*m. //' | sed 's/ .*//' > subs_production_BGS_${x}S_${y}G_${q}.txt
grep -wf exon_positions.txt subs_production_BGS_${x}S_${y}G_${q}.txt | wc -l | xargs -I [] echo "[]/(5*314*10*113)" | bc -l
		done
	done
done

for y in {BGS_nopos,neutral,neutral_nomaps}
do
        for q in {constant,expansion,contract,bottleneck,repbneck}
        do
echo -e "production_${y}_${q}"
grep -ae p1 /nas/longleaf/home/jimarsh/LOGFILES/production_${y}_${q}.o | awk '{if ($NF >= 100000) print}' | sed 's/.*m. //' | sed 's/ .*//' > subs_production_${y}_${q}.txt
grep -wf exon_positions.txt subs_production_${y}_${q}.txt | wc -l | xargs -I [] echo "[]/(5*314*10*113)" | bc -l
	done
done

##Pi in neutral intergenic/intronic regions


for y in {low,mid}
do
        for x in {mid,high}
        do
                for q in {constant,expansion,contract,bottleneck,repbneck}
                do
echo -e "production_BGS_${x}S_${y}G_${q}"
grep -ae "Rep." /nas/longleaf/home/jimarsh/LOGFILES/production_BGS_${x}S_${y}G_${q}.o | sed 's/.*= //' | sed 's/"//' | awk '{ tot+=$1;cnt++ } END { print tot/cnt }'
		done
	done
done

for y in {BGS_nopos,neutral,neutral_nomaps}
do
        for q in {constant,expansion,contract,bottleneck,repbneck}
        do
echo -e "production_${y}_${q}"
grep -ae "Rep." /nas/longleaf/home/jimarsh/LOGFILES/production_${y}_${q}.o | sed 's/.*= //' | sed 's/"//' | awk '{ tot+=$1;cnt++ } END { print tot/cnt }'
	done
done

##Pi in exon regions

#grep -ae g3 production_neutral_nomaps_bottleneck.o | sort | uniq | sed 's/.*g3, /1\t/' | sed 's/, /\t/' | sed 's/).*//' | sort -nk2 > slim_exons.bed

#module load samtools
#module load vcftools
#module load bedtools
#
#for y in {low,mid}
#do
#	for x in {mid,high}
#	do
#		for q in {constant,expansion,contract,bottleneck,repbneck}
#		do
#echo ${q}
#			for i in {1..10}
#			do
#bgzip ../production_BGS_${x}S_${y}G_${q}_R${i}.vcf
#bcftools index ../production_BGS_${x}S_${y}G_${q}_R${i}.vcf.gz
#bedtools intersect -a ../production_BGS_${x}S_${y}G_${q}_R${i}.vcf.gz -b ../../slim_exons.bed -header > g3_${i}.vcf
#vcftools --vcf g3_${i}.vcf --site-pi --out g3_${i}
#awk '{ sum += $NF } END { print sum / (5*314*10*113)}' g3_${i}.sites.pi >> production_BGS_${x}S_${y}G_${q}.exon.pi
#			done
#		done
#	done
#done
#
#for y in {BGS_nopos,neutral,neutral_nomaps}
#do
#	for q in {constant,expansion,contract,bottleneck,repbneck}
#        do
#echo ${q}
#                for i in {1..10}
#                do
#bgzip ../production_${y}_${q}_R${i}.vcf
#bcftools index ../production_${y}_${q}_R${i}.vcf.gz
#bedtools intersect -a ../production_${y}_${q}_R${i}.vcf.gz -b ../../slim_exons.bed -header > g3_${i}.vcf
#vcftools --vcf g3_${i}.vcf --site-pi --out g3_${i}
#awk '{ sum += $NF } END { print sum / (5*313*113)}' g3_${i}.sites.pi >> production_${y}_${q}.exon.pi
#                done
#        done
#done
#
#lh | grep production.*pi | sed 's/.* //' | while read line; do echo ${line}; awk '{ tot+=$1;cnt++ } END { print tot/cnt }' ${line}; done
