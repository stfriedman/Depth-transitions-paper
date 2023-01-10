library(tidyverse)
library(rfishbase)
library(picante)

select <- dplyr:::select
map <- purrr:::map

dir <- "~/Documents/Rabsoky_proj/Rabosky2019_suppmat/dataFiles"
tree <- read.tree(paste0(dir, "/trees/ultrametric_12K/actinopt_12k_treePL.tre"))

# replacing species in tree and data with accepted synonyms
syns <- read_csv("~/Documents/Rabsoky_proj/Rabosky2019_suppmat/dataFiles/taxonConversion.csv") 
colnames(syns) <- c("treeTaxon", "Species")
  
tr_syns <- synonyms(gsub("_", " ", tree$tip.label)) %>%
  as_tibble() %>%
  select(synonym, Status, Species) %>%
  mutate(synonym = gsub(" ", "_", synonym),
         Species = gsub(" ", "_", Species)) %>%
  filter(Status != "misapplied name") %>%
  
  # adding in syns list from Rabosky documentation
  left_join(syns) %>%
  mutate(Species = ifelse(!is.na(treeTaxon), treeTaxon, Species)) %>%
  arrange(match(synonym, tree$tip.label)) %>%

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

syntree <- tree    
for(i in 1:nrow(tr_syns)){
  if(syntree$tip.label %in% tr_syns$synonym){
    syntree$tip.label[which(syntree$tip.label == tr_syns$synonym[[i]])] <- tr_syns$Species[[i]]
  }
}

# species by cell occurance
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


# only grabbing every 100th cell for now to speed up results
pcd_df <- sp_occur %>%
  left_join(depthdat, by = c("tree_name" = "Species")) %>%
  filter(!duplicated(tree_name)) %>% # Gadus macrocephalus is duplicated
  select(-c(X1, DepthRangeShallow, DepthRangeDeep)) %>%
  pivot_longer(names_to = "cell", values_to = "count", -c(tree_name, depth)) %>%
  mutate(tree_name = gsub(" ", "_", tree_name)) %>%
  filter(!is.na(depth) & count > 0) %>%
  group_by(cell) %>%
  nest() %>%
  ungroup() 

saveRDS(list(pcd_df = pcd_df, syntree = syntree), 
        "~/Documents/Rabsoky_proj/slurm/pcd_input_files_updated.Rdata")

