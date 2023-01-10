library(tidyverse)
library(wesanderson)
library(rnaturalearth)

select <- dplyr::select

pal <- wes_palette("Zissou1", 100, type = "continuous")
basepath <- '~/Documents/Rabsoky_proj/Rabosky2019_suppmat/'
cell_info <- read_csv(paste0(basepath, "cellinfo_0.5.csv")) %>%
  select(cell, DepthMax)
cell_file <- read_csv(paste0(basepath, 'dataFiles/cellValues_fixed0.5.csv')) %>%
  left_join(cell_info, by = "cell")

input_files <- list.files("~/Documents/Rabsoky_proj/slurm/results",
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


pcd_bycell %>%
  filter(marineCells == 1) %>%
  filter(PCD < 4) %>%
  #filter(DepthMax > 200) %>%
  ggplot(aes(x = midLat, y = PCD)) +
  geom_point(aes(col = PCD), alpha = 0.8, cex = 2) +
  theme_classic() +
  ylab("Phylogenetic community dissimarity (photic vs. aphotic)") +
  xlab("Latitude (cell midpoint)") +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pal) 


#### plotting PCD on map projection
require(raster)
require(maptools)
require(rgeos)
require(cleangeo)
library(rnaturalearth)

template <- raster(paste0(basepath, 'dataFiles/rasterTemplate_150km.tif'))
source(paste0(basepath, 'scripts/supporting_fxns/sourceForMaps.R'))

values(template) <- pcd_bycell$PCD
values(template)[values(template) > 4] <- NA
#values(template)[pcd_bycell$DepthMax < 200] <- NA


template <- disaggregate(template, 8)
template <- unprojectMap(template, method='ngb')
ext <- projectExtent(template, crs=mollCRS)
pcd_moll <- projectRaster(template, to=ext, over=TRUE)





###########################################################################333

### interpolation of NA values
temp <- raster(paste0(basepath, 'dataFiles/rasterTemplate_150km.tif'))

values(temp) <- pcd_bycell$PCD
values(temp)[values(temp) > 4] <- NA

f <- focal(temp, w=matrix(1,nrow=3, ncol=3), fun=mean, NAonly=TRUE, na.rm=TRUE) 

intertemp <- disaggregate(f, 8)
intertemp <- unprojectMap(intertemp, method='ngb')
interext <- projectExtent(intertemp, crs=mollCRS)
inter_pcd_moll <- projectRaster(intertemp, to=interext, over=TRUE)





## plot both maps together ##
par(mfrow = c(2,1), oma=c(0.5,0.5,0.5,0.5))

## og map ##
image(pcd_moll, col=pal, axes=FALSE, xlab='', ylab='', asp=1.2)
plot(wrldMoll, col='grey40', border=NA, add=TRUE) #draws countries
lines(outerPtsMoll, xpd=NA, col = "grey40") # line around globe

# interpolated map
image(inter_pcd_moll, col=pal, axes=FALSE, xlab='', ylab='', asp=1.2)
plot(wrldMoll, col='grey40', border=NA, add=TRUE) #draws countries
lines(outerPtsMoll, xpd=NA, col = "grey40") # line around globe







##### ggplot of interpolated data 
# # reprojecting data into lat/long
# r <- unprojectMap(f) #project equal area raster to unprojected long lat
# plot(r)
# 
# #removing land from interpolated points
# land <- rasterize(wrld, template)
# # values over land need to be NA for mask to work
# values(land)[!is.na(values(land))] <- 0
# values(land)[is.na(values(land))] <- 1
# values(land)[values(land) == 0] <- NA
# plot(land)
# 
# 
# inter_res <- mask(f, land)
# plot(inter_res)
# plot(temp)
# 
# as_tibble(rasterToPoints(inter_res)) %>%
#   ggplot(aes(x = y, y = layer)) +
#   geom_point(aes(col = layer), alpha = 0.8, cex = 2) +
#   theme_classic() +
#   ylab("Phylogenetic community dissimarity (photic vs. aphotic)") +
#   xlab("Latitude (cell midpoint)") +
#   theme(legend.position = "none") +
#   scale_color_gradientn(colors = pal) 


pcd_bycell %>%
  mutate(new_PCD = values(f)) %>%
  filter(marineCells == 1) %>%
  filter(new_PCD < 4) %>%
  ggplot(aes(x = midLat, y = new_PCD)) +
  geom_point(aes(col = new_PCD), alpha = 0.8, cex = 2) +
  theme_classic() +
  ylab("Phylogenetic community dissimarity (photic vs. aphotic)") +
  xlab("Latitude (cell midpoint)") +
  theme(legend.position = "none") +
  scale_color_gradientn(colors = pal) 


