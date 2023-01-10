library(phytools)
library(tidyverse)
library(furrr)
library(progressr)

select <- dplyr:::select
map <- purrr:::map

#to run on my machine
dir <- "~/Documents/Rabsoky_proj/"
path <- "~/Documents/Rabsoky_proj/slurm/results/"
id <- 1

source(paste0(dir, "slurm/pcd_func.R"))
source(paste0(dir, "optim_pcd.R"))
d <- readRDS(paste0(dir, "slurm/pcd_input_files_updated.Rdata"))

files <- list.files(path, pattern='pcd_cell')
done_cells <- str_extract(files, "V[0-9]+")

pcd_df <- d$pcd_df %>%
  filter(!cell %in% done_cells) %>%
  mutate(nsp = map_dbl(data, nrow)) %>%
  filter(nsp > 2) %>%
  ungroup() %>%
  dplyr::select(-nsp)


plan(multisession, workers = 4) # four rsessions running simultaneously
with_progress({   #progress bar
  p <- progressor(steps = length(pcd_df))
  
  out <- pcd_df %>%
    mutate(res = future_map2(cell, data, possibly(~pcd_func(.x, .y, d$syntree, path = path), 
                                                  otherwise = NA), p = p))
})


# email me once script completed
source("~/Documents/Rfunctions/EmailfromR.R")

