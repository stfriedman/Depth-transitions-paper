##################################################
# creating unique identifier for files saved from this script
suffix <- "revisions"

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
source(here("Rcode/00_dataprep.R")) # loading prepped data


# double checking that data and phylogeny match 
stopifnot(all(tree$tip.label == alldata$sp))




# SUMMARIZE RAW DATA --------------------------------------
# script creates Table S3, raw data table for supplement, and supplemental figure summarizing depth distribution by family
source(here("Rcode/01_summ_data.R"))




# TRANSITIONS BY CLADE -------------------------------------------
# script runs the clade-level tranisiton analyses and null simulations for cladewise comparison
source(here("Rcode/02_transition_analysis.R"))




# PGLS ANALYSES ---------------------------------------------------------
# script runs both family-level and species-level pgls of depth x latitude
source(here("Rcode/03_pgls.R"))



# RATE ANALYSES ---------------------------------------------------------
# script to calculate rates of depth evolution by clade and generate Fig. S??; also contains code to compare different speciation rate metrics and generate Fig. S??
source(here("Rcode/04_rates.R"))



# PLOT FIGURE 1 ---------------------------------------------------------
# code to plot and export figure 1 from the manuscript
source(here("Rcode/plot_fig1.R"))



# PLOT FIGURE 2 ---------------------------------------------------------
# code to run Mk models + plot and export figure 2/chord plot from the manuscript and Q matrix for supplemental materials
source(here("Rcode/plot_fig2.R"))



# PCD METRIC ACROSS GRID --------------------------------------------------
# IMPORTANT: Code run on cluster not on my computer. Took many months to run and continual restarting of the job after it would time out and crash (very frustrating). Therefore, this is structured as a standalone script and is not meant to be actually run from this "run" script. Putting it here for posterity and for clarity of where this script fits into study's analysis.  
# source(here("Rcode/pcd_bycell.R"))
# source(here("Rcode/plot_fig4.R"))


# # BISSINESS ANALYSIS  --------------------------------------------------
# source(here("Rcode/05_bisseness.R"))



# COMPILES SUPPLEMENTAL MATERIALS  --------------------------------------------
rmarkdown::render(here("Rcode/supplemental_materials.Rmd"),
                    output_dir = here("results/"),
                    output_file = "supplemental_materials")

