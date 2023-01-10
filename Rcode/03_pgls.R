# SPECIES-LEVEL PGLS -------------------------------------------

gdf <- geomorph.data.frame(med_depth = setNames(alldata$med_depth, alldata$sp),
                           lat_centroid = setNames(alldata$lat_centroid, alldata$sp),
                           phy = tree)
depth_lat_pgls <- procD.pgls(med_depth ~ lat_centroid, phy = phy, 
                             data = gdf, iter = 999)

cat("species-level depth x latitude results: ")
summary(depth_lat_pgls)




# FAMILY-LEVEL PGLS -------------------------------------------
if(!exists("ta_results")){
  ta_results <- readRDS(here(paste0("data/prog_data/ta_results_", suffix, ".Rdata")))
}

stopifnot(all(ta_results$family == fam_tree$tip.label))

rate_gdf <- geomorph.data.frame(trans_rate = setNames(ta_results$trans_rate, ta_results$family),
                                med_lat =  setNames(ta_results$med_lat, ta_results$family),
                                med_lambda = setNames(ta_results$med_lambda, ta_results$family),
                                phy = fam_tree)

rate_lat_lm <- procD.lm(trans_rate ~ med_lat, data = rate_gdf, iter = 999)
rate_lat_pgls <- procD.pgls(trans_rate ~ med_lat, data = rate_gdf, 
                        phy = phy, iter = 999)
cat("family-level transition rate x latitude results: ")
summary(rate_lat_lm)
summary(rate_lat_pgls)



rate_lm <- procD.lm(trans_rate ~ med_lambda, data = rate_gdf, iter = 999)
rate_pgls <- procD.pgls(trans_rate ~ med_lambda, data = rate_gdf, 
                        phy = phy, iter = 999)

cat("family-level transition rate x speciation rate results: ")
summary(rate_lm)
summary(rate_pgls)




# FAMILY-LEVEL PGLS WITHOUT ZOARCIDS -------------------------------------------
# ta_res <- ta_results %>% 
#   filter(family != "Zoarcidae")
# tr <- drop.tip(fam_tree, "Zoarcidae")
# 
# rate_gdf_noZ <- geomorph.data.frame(trans_rate = setNames(ta_res$trans_rate, ta_res$family),
#                                 med_lambda = setNames(ta_res$med_lambda, ta_res$family),
#                                 phy = tr)
# rate_lm_noZ <- procD.lm(trans_rate ~ med_lambda, data = rate_gdf_noZ, iter = 999)
# rate_pgls_noZ <- procD.pgls(trans_rate ~ med_lambda, data = rate_gdf_noZ, 
#                         phy = phy, iter = 999)
# 
# summary(rate_lm_noZ)
# summary(rate_pgls_noZ)
