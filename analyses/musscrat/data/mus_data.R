dir <- "~/Documents/Rabsoky_proj/"
source(paste0(dir, "dataprep.R")) # loading prepped data

stopifnot(all(data$sp == tree$tip.label))

lat_cat <- matrix(data$lat_cat)
rownames(lat_cat) <- data$sp
lat_cat[,1][lat_cat[,1] == "tropical"] <- 0
lat_cat[,1][lat_cat[,1] == "temperate"] <- 1
lat_cat[,1][lat_cat[,1] == "polar"] <- 2

med_depth <- matrix(data$med_depth)
rownames(med_depth) <- data$sp


write.tree(tree, file = paste0(dir, "musscrat/data/tree.tre"))
write.nexus.data(lat_cat, file = paste0(dir, "musscrat/data/lat_cat.nex"))
write.nexus.data(med_depth, file = paste0(dir, "musscrat/data/med_depth.nex"))
