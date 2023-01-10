##################################################
# creating unique identifier for files saved from this script
suffix <- "avatree"

# toggle to rerun all or just load cached results
rerun <- FALSE
##################################################


# LOAD LIBRARIES --------------------------------------

PKG <- c(
        # phylogenetics packages
        "phytools",
         "geiger",
         "geomorph",
         "rfishbase",
         "diversitree",
         
         #tidy/data formatting packages
         "tidyverse",
         "here",
         "janitor",
         "kableExtra",
         
         #plotting packages
         "wesanderson",
         "ggrepel",
         "cowplot",
         "viridis")

PKG <- unique(PKG)
for (p in PKG) {
  if(!require(p,character.only = TRUE)) {
    install.packages(p)
    require(p,character.only = TRUE)}
}
rm(PKG, p)


# specifying packages for common functions
select <- dplyr::select
map <- purrr::map
filter <- dplyr::filter



# LOAD DATA AND FUNCTIONS --------------------------------------

source(here("Rcode/functions.R")) # loading helper functions 
source(here("Rcode/00_dataprep_avatree.R")) # loading prepped data


# double checking that data and phylogeny match 
stopifnot(all(tree$tip.label == alldata$sp))



# PLOT FIGURE 2 ---------------------------------------------------------
# code to run Mk models + plot and export figure 2/chord plot from the manuscript 
source(here("Rcode/plot_fig2_avatree.R"))


