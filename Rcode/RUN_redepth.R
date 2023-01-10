##################################################
# creating unique identifier for files saves from this script
suffix <- "redepth"

# toggle to rerun all or just load cached results
rerun <- FALSE
##################################################


# LOAD LIBRARIES --------------------------------------

PKG <- c("phytools",
         "geiger",
         "geomorph",
         "rfishbase",
         
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
source(here("Rcode/00_dataprep_redepth.R")) # loading prepped data


# double checking that data and phylogeny match 
stopifnot(all(tree$tip.label == alldata$sp))



# TRANSITIONS BY CLADE -------------------------------------------
# script runs the clade-level transition analyses and null simulations for cladewise comparison
source(here("Rcode/02_transition_analysis_redepth.R"))



# PLOT FIGURE 1 ---------------------------------------------------------
# code to plot and export figure 1 from the manuscript
source(here("Rcode/plot_fig1.R"))

