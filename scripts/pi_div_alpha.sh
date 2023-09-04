
##Alpha calculation

alpha calculate from uniq -c of m1-6

grep -m 1 -A 1000000 Mutations production_huber_BGS_highS_midG_constant.o | grep -m 1 -B 10000000 '"' | head -n -1 | awk '{print $3}' | sort | uniq -c > temp.txt
cat temp.txt | grep m6 | sed "s/ //g" | sed "s/m6//g" | xargs -I [] awk '{sum+=$1;} END {print []/sum}' temp.txt

##Divergence in neutral intergenic/intronic regions

"Theta pi for neutral interg/intron across whole genome = 0.000647246"
0.000647246*138600000 = 129819787

grep m5 /nas/longleaf/home/jimarsh/LOGFILES/production_huber_BGS_highS_midG_constant.o | grep p1 | wc -l
"1879832"

div = (1879832/10)/129819787

##Divergence in exonic regions

# Read the content of file 1 and extract the two columns as ranges
with open('file1.txt', 'r') as file1:
    lines = file1.readlines()
    ranges = [(int(line.split()[1]), int(line.split()[2])) for line in lines]

# Read the content of file 2 and extract the integers
with open('file2.txt', 'r') as file2:
    integers = [int(line.strip()) for line in file2.readlines()]

# Iterate through the integers and check if they are within any of the ranges
for num in integers:
    for range_start, range_end in ranges:
        if range_start <= num < range_end:
            print(num)
            break  # No need to check other ranges if the number is already found within one

grep p1 /nas/longleaf/home/jimarsh/LOGFILES/production_huber_BGS_highS_midG_constant.o | sed 's/[^ ]* //' | sed 's/ m.*//' > test_subs.txt

div = ans/(11*979*311)

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

