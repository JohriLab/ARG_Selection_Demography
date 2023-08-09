grep chr7 recs.txt | awk '$3 <= 152054332 {print $0}' | awk '$2 >= 6010001 {print $0}' | awk '$2 = $2-6010000 {print $0}' | awk '$3 = $3-6010000 {print $0}' > chr7_trimmed_recs_4col.txt 

awk '{print $3, $4}' chr7_trimmed_recs_4col.txt  > chr7_trimmed_recs.txt
