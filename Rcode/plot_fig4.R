library(tidyverse)
library(wesanderson)
library(rnaturalearth)
library(here)
require(raster)
require(maptools)
require(rgeos)
require(cleangeo)

pal <- wes_palette("Zissou1", 100, type = "continuous")
basepath <- here("data/Rabosky2019_suppmat/")
cell_info <- read_csv(paste0(basepath, "cellinfo_0.5.csv")) %>%
  select(cell, DepthMax)
cell_file <- read_csv(paste0(basepath, 'dataFiles/cellValues_fixed0.5.csv')) %>%
  left_join(cell_info, by = "cell")

input_files <- list.files(here("analyses/slurm/results"),
                          pattern='pcd_cell', full.names=TRUE)
pcd_df <- tibble(file = input_files) %>%
  mutate(data = purrr:::map(file, readRDS)) %>%
  dplyr:::select(-file) %>%
  unnest()


pcd_bycell <- pcd_df %>%
  unnest() %>%
  pivot_wider(names_from = metric, values_from = values) %>%
  
  # making PCD 0 when all species shared between photic/aphotic
  mutate(PCD = ifelse(is.na(PCDc) & !is.na(PCD), 0, PCD)) %>%
  
  mutate(cell = gsub("V", "cell", cell)) %>%
  right_join(cell_file, by = "cell") %>%

  #mutate(values = ifelse(marineCells == 0, NA, values)) %>%
  arrange(match(cell, cell_file$cell))


# ## species in each cell
# species_occ <- read_csv(paste0(basepath, 'dataFiles/species_by_cell_PAmat_merged.csv')) 
# 
# cell_sp <- species_occ %>%
#   rename(species = `...1`) %>%
#   pivot_longer(-species, names_to = "cell", values_to = "presence") %>%
#   filter(presence != 0) %>%
#   select(-presence) %>%
#   group_by(cell) %>%
#   mutate(species = list(as.character(species)),
#          cell = gsub("V", "cell", cell)) %>%
#   unique() %>%
#   arrange(cell)


# # does species richness drive PCD? NOPE
# d <- readRDS(here("analyses/slurm/pcd_input_files.Rdata"))
# pcd_df <- d$pcd_df %>%
#   unnest() %>%
#   filter(count > 0) %>%
#   mutate(cell = gsub("V", "cell", cell)) %>%
#   group_by(cell) %>%
#   add_count() %>%
#   select(cell, nsp = n) %>%
#   unique()
# 
# pcd_bycell %>%
#   left_join(pcd_df, by = "cell") %>%
#   ggplot(aes(x = nsp, y = PCD)) +
#   geom_point() +
#   theme_classic()


# plotting PCD on map projection
template <- raster(paste0(basepath, 'dataFiles/rasterTemplate_150km.tif'))
source(paste0(basepath, 'scripts/supporting_fxns/sourceForMaps.R'))

values(template) <- pcd_bycell$PCD
values(template)[values(template) > 4] <- NA
#values(template)[pcd_bycell$DepthMax < 200] <- NA

template <- disaggregate(template, 8)
template <- unprojectMap(template, method='ngb')
ext <- projectExtent(template, crs=mollCRS)
pcd_moll <- projectRaster(template, to=ext, over=TRUE)


# interpolation of cells that were computationally limited; estimating as mean of surrounding cells
temp <- raster(paste0(basepath, 'dataFiles/rasterTemplate_150km.tif'))

values(temp) <- pcd_bycell$PCD
values(temp)[values(temp) > 4] <- NA

f <- focal(temp, w=matrix(1,nrow=3, ncol=3), fun=mean, NAonly=TRUE, na.rm=TRUE) 

intertemp <- disaggregate(f, 8)
intertemp <- unprojectMap(intertemp, method='ngb')
interext <- projectExtent(intertemp, crs=mollCRS)
inter_pcd_moll <- projectRaster(intertemp, to=interext, over=TRUE)



# PLOTTING RESULTS -----------------------------------
png("results/figures/Fig4_globe.png", 
    width     = 3.25,
    height    = 3.25,
    units     = "in",
    res       = 1200,
    pointsize = 4)
par(mar = rep(0.5, 4), oma = rep(0.5, 4))
image(inter_pcd_moll, col=pal, axes=FALSE, xlab='', ylab='', asp=1.2)
plot(wrldMoll, col='grey40', border=NA, add=TRUE) #draws countries
lines(outerPtsMoll, xpd=NA, col = "grey40") # line around globe
dev.off()


# reading globe in as a gg obj
globe <- ggdraw() +
  draw_image(
    "results/figures/Fig4_globe.png", scale = 1
  ) 


# function for rescaling values so color is scaled appropriately on the facetted plot
scale_func <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


p1 <- pcd_bycell %>%
  filter(marineCells == 1) %>%
  filter(PCD < 4) %>%  #removing anomalous cell22349, which only has 4 spp. 
  ggplot(aes(x = midLat, y = PCD)) +
  geom_point(aes(col = PCD), alpha = 0.6, cex = 1.5) +
  theme_classic() +
  theme(aspect.ratio = 1.3,
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.1),
        axis.title = element_text(size = 12)) +
  ylab("Phylogenetic community dissimarity (photic vs. aphotic)") +
  xlab("Latitude (cell midpoint)") +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pal) 


pdf("results/figures/Fig4.pdf")
plot_grid(p1, globe, ncol = 2, align = "h", rel_widths = c(1, 0.8))
dev.off()
