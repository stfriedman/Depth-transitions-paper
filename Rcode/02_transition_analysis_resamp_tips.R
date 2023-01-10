# FAMILY-LEVEL DATASET -------------------------------------------

fam_data <- trim_data %>%
  select(perc, iter, species = sp, family = family, lat_centroid, 
         tv_lambda, tridepth, nsp, tree) %>%
  group_by(perc, iter, family) %>%
  mutate(species = list(as.character(species))) %>%
  add_count(tridepth, name = "ncat") %>%
  mutate(med_lat = median(lat_centroid),
         med_lambda = median(tv_lambda)) %>%
  select(-lat_centroid, -tv_lambda) %>%
  unique() %>%
  nest(depth = c(tridepth, ncat)) 




# GETTING MRCAS FOR EACH CLADE -------------------------------------------

mrcas <- fam_data %>%
  mutate(mrca = map2_dbl(tree, species, getMRCA),
         age = map2_dbl(tree, mrca, ~branching.times(.x)[[as.character(.y)]]))



# MAKING SIMMAPS AND COUNTING TRANSITIONS -------------------------------------------

simtree_foo <- function(data, tree){
  tmp <- data %>%
    arrange(match(sp, tree$tip.label)) 
  
  stopifnot(all(tmp$sp == tree$tip.label))
  
  cat_depths <- setNames(tmp$tridepth, tmp$sp)
  simtree <- make.simmap(tree, cat_depths, model = "ARD", nsim = 1)
  simtree
}


#if(rerun){
  alldata_simtree <- alldata %>%
    mutate(simtree = map2(samp, tree, simtree_foo))
  
  saveRDS(alldata_simtree, here(paste0("data/prog_data/depthsimtree_resamp_tips.Rdata")))
#}



# counting transitions by clade
clade_transitions <- alldata_simtree %>%
  left_join(mrcas, by = c("perc", "iter")) %>%
  group_by(perc, iter) %>%
  nest(mrcas = c(species:age)) %>%
  rename(tree = tree.x) %>%
  mutate(trans = map2(mrcas, simtree, ~num_trans(.x$mrca, .y)))


summ_clade_transitions <- clade_transitions %>%
  select(perc, iter, trans) %>%
  unnest(cols = c(trans)) %>%
  left_join(mrcas, by = c("perc", "iter", "mrca")) %>%
  select(-desc, -tree) %>%
  group_by(perc, family) %>%
  summarize(trans_num = median(trans_num, na.rm = TRUE)) 




# RUNNING SIMULATIONS -----------------------------------------
nsims <- 1 #simmaps for each tree
sim_x <- clade_transitions %>%
  select(perc, iter, tree, simtree, mrcas) %>%
 # filter(iter == 1)  %>% # just running under first Qmat for each percentage sampling
  mutate(q_emp = map(simtree, ~.$Q),
         simchar = map2(tree, q_emp, ~sim.char(.x, .y, model = "discrete", nsim = nsims))) %>%
  ungroup()


# running simulations of transitions under empirical Q matrix; takes a long time
sim_output <- sim_x %>%
  select(-simtree) %>%
  select(perc, tree, iter, simchar) %>%
  expand_grid(sim_id = 1:nsims) %>%
  mutate(sim_x = map2(simchar, sim_id, ~.x[,,.y]),
         simmap = map2(tree, sim_x, make.simmap))
  
saveRDS(sim_output, here(paste0("data/prog_data/sim_res_resamp_tips.Rdata")))



# summarizing results of simulations
sim_result <- sim_output %>%
  ungroup() %>%
  left_join(mrcas, by = c("perc", "iter")) %>%
  group_by(perc, iter) %>%
  nest(mrcas = c(species:age)) %>%
  rename(tree = tree.x) %>%
  mutate(trans = map2(mrcas, simmap, ~num_trans(.x$mrca, .y)))


sim_result_summ <- sim_result %>%
  select(perc, iter, trans) %>%
  unnest(cols = c(trans)) %>%
  left_join(mrcas, by = c("perc", "iter", "mrca")) %>%
  ungroup() %>%
  select(perc, family, trans_num, age, med_lat, med_lambda) %>%
  group_by(perc, family) %>%
  mutate(sim_0.05 = quantile(trans_num, probs = 0.05, na.rm = TRUE),
         sim_median = median(trans_num, na.rm = TRUE),
         sim_0.95 = quantile(trans_num, probs = 0.95, na.rm = TRUE),
         age = mean(age, na.rm = TRUE),
         med_lat = mean(med_lat, na.rm = TRUE),
         med_lambda = mean(med_lambda, na.rm = TRUE)) %>%
  select(-trans_num) %>%
  unique() 




# RESULTS OF CLADEWIDE TRANSITION ANALYSIS  -----------------------------------------
if(rerun){
ta_results <- left_join(summ_clade_transitions, sim_result_summ, by = c("perc", "family")) %>%
    mutate(n_diff = trans_num - sim_median,
         trans_rate = trans_num/age,
         col = case_when(
           trans_num >= sim_0.05 & trans_num <= sim_0.95 ~ "mid",
           trans_num > sim_0.95 ~ "big",
           trans_num < sim_0.05 ~ "small"),
         famcol = ifelse(family %in% rab_fams, "#F8766D", "grey40"))
saveRDS(ta_results, here(paste0("data/prog_data/ta_results_resamp_tips.Rdata")))
} else {
  ta_results <- readRDS(here(paste0("data/prog_data/ta_results_resamp_tips.Rdata")))
}
