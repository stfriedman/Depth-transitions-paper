args <- as.integer(commandArgs(trailingOnly = TRUE))
id <- args[1]

print(sprintf("SLURM ARRAY ID: %d", id))

library(phytools)
library(tidyverse)

select <- dplyr:::select
map <- purrr:::map


dir <- "/"
path <- "results/"

# #to run on my machine
# dir <- "~/Documents/Rabsoky_proj/"
# path <- "~/Documents/Rabsoky_proj/slurm/results/"
# id <- 1

source(paste0(dir, "slurm/pcd_func.R"))
source(paste0(dir, "optim_pcd.R"))
d <- readRDS(paste0(dir, "slurm/pcd_input_files_updated.Rdata"))

files <- list.files(path, pattern='pcd_cell')
done_cells <- str_extract(files, "V[0-9]+")

chunk <- 1
pcd_df <- d$pcd_df %>%
  filter(!cell %in% done_cells) %>%
  mutate(nsp = map_dbl(data, nrow)) %>%
  filter(nsp > 2) %>%
  dplyr::select(-nsp) %>%
  mutate(chunk = rep(1:ceiling(nrow(.)/chunk),each=chunk)[1:nrow(.)]) %>%
  split(.$chunk)

out <- pcd_df[[id]] %>%
  mutate(res = map2(cell, data, possibly(~pcd_func(.x, .y, d$syntree, path = path), 
                                         otherwise = NA)))
