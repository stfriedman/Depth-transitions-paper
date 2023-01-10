# Code run on cluster not on my computer. Took months to run and continual restarting of the job after it would time out and crash (very frustrating). Therefore, this is a standalone script. 

library(tidyverse)
library(rfishbase)
library(picante)
library(furrr)
library(here)

select <- dplyr:::select
map <- purrr:::map

source(here("Rcode/functions.R")) # loading helper functions 

dir <- here("Rabosky2019_suppmat/dataFiles")
rawtree <- read.tree(paste0(dir, "/trees/ultrametric_12K/actinopt_12k_treePL.tre"))


# replacing species in tree and data with accepted synonyms
syns <- read_csv(paste0(dir, "taxonConversion.csv"))
colnames(syns) <- c("treeTaxon", "Species")
  
tr_syns <- synonyms(gsub("_", " ", rawtree$tip.label)) %>%
  as_tibble() %>%
  select(synonym, Status, Species) %>%
  mutate(synonym = gsub(" ", "_", synonym),
         Species = gsub(" ", "_", Species)) %>%
  filter(Status != "misapplied name") %>%
  
  # adding in syns list from Rabosky documentation
  left_join(syns) %>%
  mutate(Species = ifelse(!is.na(treeTaxon), treeTaxon, Species)) %>%
  arrange(match(synonym, rawtree$tip.label)) %>%

  # dealing with duplicated synonyms 
  distinct(synonym, Species, .keep_all = TRUE) %>%
  group_by(synonym) %>%
  add_count() %>%
  mutate(filt = ifelse(n > 1 & synonym != Species & Status != "accepted name", T, F),
         filt = replace(filt, 
                          synonym %in% c("Apogon_striatus", "Chrysophrys_auratus", "Barbus_miolepis") &
                            Status == "synonym", F)) %>%
  filter(filt == F) %>%
  select(-treeTaxon, -n, -filt) 

syntree <- rawtree    
for(i in 1:nrow(tr_syns)){
  if(syntree$tip.label %in% tr_syns$synonym){
    syntree$tip.label[which(syntree$tip.label == tr_syns$synonym[[i]])] <- tr_syns$Species[[i]]
  }
}

# species by cell occurrence
sp_occur <- read_csv(paste0(dir, "/species_by_cell_PAmat_merged.csv")) %>%
  left_join(tr_syns, by = c("X1" = "synonym")) %>%
  select(X1, tree_name = Species, everything(), -Status) %>%
  filter(tree_name %in% syntree$tip.label) %>%
  mutate(tree_name = gsub("_", " ", tree_name)) 

# getting species depth and taxonomic info
depthdat <- species(sp_occur$tree_name, fields = c("Species", "DepthRangeShallow", "DepthRangeDeep")) %>%
  mutate(depth = case_when(
    DepthRangeDeep <= 200 ~ "photic",
    DepthRangeShallow >= 200 ~ "aphotic",
    DepthRangeDeep >= 200 & DepthRangeShallow <= 200 ~ "both")) 


# running phylogenetic dissimilarity metric
pcd_func <- function(cell, data, tree, path, run_dups = FALSE){
  if(!hasArg(path)) path <- getwd()
  if(run_dups == FALSE){
    input_files <- list.files(path, pattern='pcd_cell')
    done_cells <- str_extract(input_files, "V[0-9]+")
  }
  
  # trying to save time and not run duplicate cells
  if(cell %in% done_cells & run_dups == FALSE){
    cat(paste0(cell, ": getting data from existing file\n"))
    out <- readRDS(paste0(path, "/pcd_cell_", cell, ".Rdata")) %>%
      select(res) %>%
      unnest()
  } else {
    df <- data %>% 
      filter(count != 0) %>%
      mutate(
        photic = ifelse(depth %in% c("photic", "both"), 1, 0),
        aphotic = ifelse(depth %in% c("aphotic", "both"), 1, 0)
      ) %>%
      select(tree_name, photic, aphotic)
    df_tree <- drop.tip(tree, setdiff(tree$tip.label, df$tree_name))
    
    if(length(df$tree_name) > 0){
    x <- df %>%
      arrange(match(tree_name, df_tree$tip.label)) %>%
      pivot_longer(photic:aphotic, names_to = "area", values_to = "count") %>%
      pivot_wider(names_from = tree_name, values_from = count) %>%
      column_to_rownames("area") %>%
      as.matrix()
    
    cat(paste0(cell, ": calculating pcd...\n"))
    res <- optim_pcd(x, df_tree)
    out <- tibble(metric = c("PCD", "PCDc", "PCDp"),
                  values = c(res$PCD[[1]], res$PCDc[[1]], res$PCDp[[1]]))
    } else { 
    #if cell has no species
     out <- tibble(metric = c("PCD", "PCDc", "PCDp"),
                   values = rep(NA, 3))
    }
    save <- tibble(cell = cell, nest(out, .key = "res"))
    saveRDS(save, paste0(path, "/pcd_cell_", cell, ".Rdata"))
  }
  out
}


# only grabbing every 100th cell for now to speed up results
pcd_df <- sp_occur %>%
  left_join(depthdat, by = c("tree_name" = "Species")) %>%
  filter(!duplicated(tree_name)) %>% # Gadus macrocephalus is duplciated
  select(-c(X1, DepthRangeShallow, DepthRangeDeep)) %>%
  pivot_longer(names_to = "cell", values_to = "count", -c(tree_name, depth)) %>%
  mutate(tree_name = gsub(" ", "_", tree_name)) %>%
  filter(!is.na(depth)) %>%
  group_by(cell) %>%
  nest() %>%
  ungroup() %>%
  mutate(cell_num = 1:nrow(.)) 

path <- here("data/PCD_output_bycell")
pcd_df %>%
  filter(cell_num %% 100 == 0) %>%
  select(-cell_num) %>%
  mutate(res = future_map2(cell, data, possibly(~pcd_func(.x, .y, syntree), otherwise = NA)))


# # to run all cells
# pcd_df <- sp_occur %>%
#   rename(Species = X1) %>%
#   left_join(depthdat) %>%
#   pivot_longer(names_to = "cell", values_to = "count", -c(Species, depth)) %>%
#   mutate(Species = gsub(" ", "_", Species)) %>%
#   filter(!is.na(depth)) %>%
#   group_by(cell) %>%
#   nest() %>%
#   mutate(res = future_map2(cell, data, possibly(~pcd_func(.x, .y, tree), otherwise = NA)))

#saveRDS(pcd_df, here("data/pcd_df.Rdata"))

pcd_df %>%
  select(cell, res) %>%
  unnest() %>%
  filter(metric == "PCD")


# input_files <- list.files(here("data/PCD_output_bycell"),
#                           pattern='pcd_cell', full.names=TRUE)
# pcd_df <- tibble(file = input_files) %>%
#   mutate(data = purrr:::map(file, readRDS)) %>%
#   select(-file) %>%
#   unnest() 
# saveRDS(pcd_df, here("data/pcd_df.Rdata"))
