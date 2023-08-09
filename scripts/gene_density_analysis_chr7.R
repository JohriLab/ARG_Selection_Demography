library(tidyverse)

# Read in gene search results from NCBI
# Searched "10[Chr] AND "Homo sapiens"[Organism] AND ("has ccds"[Properties] AND alive[prop])"
qed <- read_tsv("gene_result-3.txt") %>% 
  select("start_position_on_the_genomic_accession","end_position_on_the_genomic_accession", "exon_count", "orientation") %>% 
  rename(start = "start_position_on_the_genomic_accession", end = "end_position_on_the_genomic_accession") %>% 
# Calculate midpoint of each gene, and length of each gene
  mutate("mid" = (as.numeric(start) + as.numeric(end)) / 2) %>% 
  mutate("length" = abs(as.numeric(start) - as.numeric(end))) %>% 
  arrange(., end)

# % of chromosome that's genic
message(paste0(round(
  sum(qed$length, na.rm = T)/250000000*100, 
  digits = 2), "% of chromosomal positions are within genes"))

# Can simply plot gene density as histogram across chromosome using midpoints
ggplot2::ggplot(data=qed, aes(x = mid)) +
  geom_histogram(bins = 50)

# Segment chromosome into 1Mb bins, counting nGenes in each
binned <- 
  table(qed$mid %>% cut(breaks = seq(0, max(.), by =1000000))) %>% 
  as_tibble(.name_repair = "unique") %>% 
  rownames_to_column("bin") %>% 
  mutate(bin = as.integer(bin)) %>%
  separate("...1", into = c("start", "end"), sep = ",") %>%
  mutate(start = as.numeric(gsub("\\(","",start))) %>% 
  mutate(end = as.numeric(gsub("]","",end))) %>% 
  mutate(meanlen = round(1000000/n))
binned
  
# Summaries used for building distribution to sample from in SLiM
length(binned$n)
mean(binned$n)
median(binned$n)
sd(binned$n)


length(binned)

# Plot similarly to above, but using bins rather than histogram
binplot <- ggplot2::ggplot(data=binned, aes(x = bin, y = n)) +
  geom_bar(stat = "identity")
binplot
