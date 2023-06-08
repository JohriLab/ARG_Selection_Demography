#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

nSamples <- as.integer(args[1]) #e.g. 10
directory <- as.character(args[2]) #e.g. 'foo/bar'
prefix <- as.character(args[3]) #e.g. 'N250'

mywd <- getwd()
if(gsub('.*/','', mywd) != directory && '.' != directory){
setwd(directory)
}

library(tidyverse)

#for i in {1..10}
#do
#cat N250_R${i}_popsize.coal | tail -n +2 | sed 's/^0 0 //' | sed 's/ $//g' | sed 's/ /,/g' > N250_R${i}.recoal
#done
dat_Ne <- tibble::tibble()
for (i in 1:nSamples){
dat_Ne  <- readr::read_csv(file = base::paste0(prefix,'_R',i,'.recoal'),
           col_names = F, show_col_types = F) %>%
    t() %>% tibble::as_tibble() %>%
    dplyr::rename("nGen" = "V1", "Coal" = "V2") %>%
    dplyr::mutate(!!paste0("R",i) := 0.5 / Coal) %>%
    dplyr::select(-Coal)  %>%
    tidyr::gather("Replicate", "Ne", 2) %>%
    base::rbind(dat_Ne, (.))
}

allNe_plot <- ggplot2::ggplot(dat_Ne) +
  ggplot2::geom_line(ggplot2::aes(x = nGen, y = Ne, colour = Replicate)) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10()

ggplot2::ggsave(paste0(prefix,'_all.pdf'),
                allNe_plot,
                device = 'pdf',
                dpi = 300,
                height = 9,
                width = 16,
                units = "cm")
