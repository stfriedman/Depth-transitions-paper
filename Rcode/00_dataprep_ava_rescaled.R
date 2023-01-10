rerun_scale <- FALSE

# defining basepath for Rabosky supplemental materials 
basepath <- here("data/Rabosky2019_suppmat/dataFiles")



# RESCALING PHYLOGENY --------------------------------------------------

if(rerun_scale){
# loading Rabosky tree
rabtree <- read.tree(here(basepath, "trees/ultrametric_12K/actinopt_12k_treePL.tre")) %>%
  force.ultrametric()


# loading Ava tree
rawavatree <- read.nexus(here("data/Avatree.tre")) 
rawavatree$tip.label <- str_extract(str_to_title(rawavatree$tip.label), "[A-Za-z]+_[a-z]+")
dups <- avatree$tip.label[which(duplicated(rawavatree$tip.label))]
avatree <- drop.tip(rawavatree, c(setdiff(rabtree$tip.label, rawavatree$tip.label), dups)) %>%
  ladderize()


# ref tree to genus-level
# finding synonyms to collapse taxa correctly
avasyns <- synonyms(gsub("_", " ", avatree$tip.label)) %>%
  clean_names() %>%
  filter(!grepl("misapplied|ambiguous", status)) %>%
  select(fb_species = species, tree_name = synonym) %>%
  mutate(species = ifelse(is.na(fb_species), tree_name, fb_species)) %>%
  filter(!duplicated(species)) %>%
  select(tree_name, species)

tax <- load_taxa() %>% 
  as_tibble() %>%
  clean_names() %>%
  select(genus, family, order) %>%
  unique()

avatax <- avasyns %>%
  mutate(species = gsub(" ", "_", species),
         genus = str_extract(species, "^[A-Za-z]+")) %>%
  left_join(tax, by = "genus") %>%
  filter(!duplicated(species)) %>%
  select(-species) %>%
  mutate(tree_name = gsub(" ", "_", tree_name),
         family = ifelse(family %in% c("Gobiidae", "Serranidae"), #removing non-mono grps
                         paste0("DELETE", 1:20),
                         family)
         ) %>%
  mutate(family = case_when(
    genus == "Emblemariopsis" ~ "Chaenopsidae",
    family == "Scaridae" ~ "Labridae",
    TRUE ~ family
  )) %>%
  column_to_rownames("tree_name") %>%
  as.matrix


# modified code to drop taxa from non-monophyletic groups and taxa without taxonomic info
ref <- geiger:::subset.phylo(avatree, tax = avatax, rank = "family", drop.tax = TRUE)
ref <- drop.tip(ref, ref$tip.label[grep("DELETE", ref$tip.label)])


# installed treePL via homebrew, but wasn't able to add path to system environment like the help file for congruify.phylo suggests, so I moved the executable to an existing path that I found here: Sys.getenv("PATH") and it worked!
res_tax <- tibble(species = rabtree$tip.label) %>%
  mutate(species = gsub("_", " ", species)) %>%
  left_join(load_taxa() %>% 
              as_tibble() %>%
              clean_names(), by = "species") %>%
  select(species, genus, family, order) %>%
  mutate(species = gsub(" ", "_", species)) %>%
  column_to_rownames("species") %>%
  as.matrix()


  restree <- congruify.phylo(ref, rabtree, taxonomy = res_tax, scale = "treePL")
  saveRDS(restree, here("data/prog_data/rescaledAva_rabosky.Rdata"))
} else {
  restree <- readRDS(here("data/prog_data/rescaledAva_rabosky.Rdata"))
}


wholetree <- restree$phy





# LOAD DATA AND FILTER --------------------------------------
rab_fams <- c("Notothenioid", "Sebastidae", "Liparidae", "Zoarcidae", "Pleuronectidae")

rawdata <- read_csv(here(basepath, "rate_lat_stats_by_sp_fixed0.5.csv")) %>%
  clean_names() %>%
  # only marine species in the phylogeny
  filter(freshwater == 0 & marine == 1 |
           freshwater == 1 & marine == 1) %>%
  filter(!is.na(depth_range_shallow) & !is.na(depth_range_deep)) %>%
  filter(sp %in% wholetree$tip.label) %>%
  arrange(match(sp, wholetree$tip.label)) %>%
  
  # default depth data wrong for these species, updated with info in species notes on fishbase
  mutate(depth_range_deep = replace(depth_range_deep, sp == "Gymnoscopelus_hintonoides", 800),
         depth_range_shallow = replace(depth_range_shallow, sp == "Gymnoscopelus_hintonoides", 328)) %>%
  mutate(depth_range_deep = replace(depth_range_deep, sp == "Lampanyctus_jordani", 1000),
         depth_range_shallow = replace(depth_range_shallow, sp == "Lampanyctus_jordani", 0)) %>%
  
  # creating depth categories
  mutate(tridepth = case_when(
    depth_range_deep > 0 & depth_range_deep < 200 ~ "shallow",
    depth_range_deep >= 200 & depth_range_deep < 1000 ~ "intermediate",
    depth_range_deep >= 1000 ~ "deep"
  )) %>%
  rowwise() %>%
  dplyr::mutate(med_depth = median(c(depth_range_shallow, depth_range_deep), na.rm = TRUE)) %>%

  #creating latitude categories
  mutate(lat_cat = case_when(
    lat_centroid >-30 & lat_centroid < 30 ~ "tropical",
    abs(lat_centroid) > 30 & abs(lat_centroid) < 60 ~ "temperate",
    abs(lat_centroid) > 60 ~ "polar"
  )) 



# families in notothenioids
notothen <- c("Nototheniidae", "Bathydraconidae", "Channichthyidae", "Artedidraconidae", "Harpagiferidae")


# removing problematic species and fixing some taxonomic issues
alldata <- rawdata %>%
  
  # removing species with known placement issues in tree
  filter(!sp %in% c("Niphon_spinosus", "Trachyrincus_scabrus", "Trachyrincus_murrayi",
                     "Squalogadus_modificatus", "Triplophos_hemingi")) %>% 

  # making notothen clade same as Rabosky et al. 2018
  mutate(family = replace(family, family %in% notothen, "Notothenioid")) %>%
  
  # fixing incorrect family
  mutate(family = replace(family, sp == "Nautichthys_pribilovius", "Hemptripteridae"),
         family = replace(family, family == "Engraulidae", "Clupeidae")) %>%
  
  # filtering out these genera to make pleuro clade same as Rabosky et al. 2018,
  # last two genera removed to make Paralichthyidae monophyletic
  filter(!grepl("Pleuronichthys|Hypsopsetta|Atheresthes|Cyclopsetta|Syacium", sp)) 



# MATCHING TREE TO UPDATED DATASET --------------------------------------
tree <- drop.tip(wholetree, setdiff(wholetree$tip.label, alldata$sp))
stopifnot(all(tree$tip.label == alldata$sp))



# trimmed dataset for family-level analyses --------------------------------------
trim_data <- alldata %>%
  group_by(family) %>%
  add_count(name = "nsp") %>%
  
  # removing families with less than 20 species
  filter(nsp > 20) %>%
  
  #removing paraphyletic families
  filter(!family %in% c("Chaenopsidae", "Congridae", "Cynoglossidae",
                        "Stichaeidae", "Scorpaenidae"))

trim_tree <- drop.tip(tree, setdiff(tree$tip.label, trim_data$sp))
stopifnot(all(trim_tree$tip.label == trim_data$sp))



# removing unnecessary objects
rm(rawdata, wholetree, basepath)

