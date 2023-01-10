##################################################
# creating unique identifier for files saved from this script
suffix <- "resamp_tips"

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
source(here("Rcode/00_dataprep_resamp_tips.R")) # loading prepped data




# TRANSITIONS BY CLADE -------------------------------------------
# script runs the clade-level tranisiton analyses and null simulations for cladewise comparison
source(here("Rcode/02_transition_analysis_resamp_tips.R"))




# PGLS ANALYSES ---------------------------------------------------------
# script runs both family-level and species-level pgls of depth x latitude
source(here("Rcode/03_pgls.R"))



# PLOT FIGURE 1 ---------------------------------------------------------
# code to plot and export figure 1 from the manuscript
source(here("Rcode/plot_fig1_resamp_tips.R"))

