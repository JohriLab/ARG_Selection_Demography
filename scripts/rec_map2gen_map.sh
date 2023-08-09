#track type=wiggle_0 name="HapMap YRI" description="HapMap Release 24 YRI recombination map"


cat chr7_trimmed_recs_4col.txt | awk '{scaledstart =$2*(138600000/146042804)} {scaledstop = $3*(138600000/146042804)} {print $1" "int(scaledstart)" "int(scaledstop)" "$4}' > scaled_chr7.txt
awk '$2 != $3' scaled_chr7.txt | awk '{cM = ($3 - $2)*($4*1e-8)*100} {cMMb = cM/(($3-$2)*1e-6)} {print $2" "cMMb" "cM}' | sed "$ s/.*/138600000 0 0/" > 138Mb_r_rate.map
awk '{sum += $3} { OFMT = "%6.8f" ;print sum}' 138Mb_r_rate.map | paste 138Mb_r_rate.map - | awk '{print $1,$2,$4}' > 138Mb.map
