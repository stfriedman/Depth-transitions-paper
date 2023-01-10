# FAMILY-LEVEL DATASET -------------------------------------------

trim_taxonomy <- trim_data %>% 
  select(sp, family) %>%
  column_to_rownames("sp") %>%
  data.frame()

stopifnot(all(rownames(trim_taxonomy) == trim_tree$tip.label))


# making family-level tree 
fam_tree <- geiger:::subset.phylo(trim_tree, tax = trim_taxonomy, rank = "family")

fam_data <- trim_data %>%
  select(species = sp, family = family, lat_centroid, 
         tv_lambda, tridepth, nsp) %>%
  group_by(family) %>%
  mutate(species = list(as.character(species))) %>%
  add_count(tridepth, name = "ncat") %>%
  mutate(med_lat = median(lat_centroid),
         med_lambda = median(tv_lambda)) %>%
  select(-lat_centroid, -tv_lambda) %>%
  unique() %>%
  nest(depth = c(tridepth, ncat)) %>%
  arrange(match(family, fam_tree$tip.label)) 

stopifnot(all(fam_data$family == fam_tree$tip.label))



# GETTING MRCAS FOR EACH CLADE -------------------------------------------
mrcas <- fam_data %>%
  mutate(mrca = map_dbl(species, ~getMRCA(tree, .)),
         age = map_dbl(mrca, ~branching.times(tree)[[as.character(.)]]))



# MAKING SIMMAPS  -------------------------------------------

# making simmaps 
stopifnot(all(alldata$sp == tree$tip.label))
cat_depths <- setNames(alldata$tridepth, alldata$sp)

if(rerun){
  simtree <- make.simmap(tree, cat_depths, model = "ARD", nsim = 100)
  saveRDS(simtree, here(paste0("data/prog_data/depthsimtree_", suffix, ".Rdata")))

# need to reformat to put into tibble format
simtree <- structure(simtree, class = c("list","multiSimmap","multiPhylo"))


# counting transitions by clade
clade_transitions <- tibble(tree_id = 1:length(simtree)) %>%
  mutate(simmap = simtree,
         trans = map(simmap, ~num_trans(mrcas$mrca, .x))) %>%
  unnest(trans) %>%
  select(-simmap, -desc) %>%
  group_by(mrca) %>%
  summarize(trans_num = median(trans_num)) 
}


# RUNNING SIMULATIONS -----------------------------------------

# running simulations of transitions under empirical Q matrix; takes a long time
if(rerun){
  nsims <- 100
  q_emp <- simtree[[1]]$Q
  sim_x <- sim.char(tree, q_emp, model = "discrete", nsim = nsims)
  sim_output <- tibble(sim_id = 1:nsims) %>%
    mutate(x = map(sim_id, ~sim_x[,,.]),
           simmap = map(x, ~make.simmap(tree, .)))
  saveRDS(sim_output, here(paste0("data/prog_data/sim_res_", suffix, ".Rdata")))


# summarizing results of simulations
sim_result <- sim_output %>% 
  mutate(trans = map(simmap, ~num_trans(mrcas$mrca, .x))) %>%
  select(-simmap, -x) %>% 
  unnest() 

sim_result_summ <- sim_result %>%
  select(mrca, trans_num) %>%
  group_by(mrca) %>%
  mutate(sim_0.05 = quantile(trans_num, probs = 0.05),
         sim_median = median(trans_num),
         sim_0.95 = quantile(trans_num, probs = 0.95)) %>%
  select(-trans_num) %>%
  unique() 
}



# RESULTS OF CLADEWIDE TRANSITION ANALYSIS  -----------------------------------------
if(rerun){
ta_results <- left_join(clade_transitions, sim_result_summ, by = "mrca") %>%
  left_join(mrcas, by = "mrca") %>%
  unique() %>%
  mutate(n_diff = trans_num - sim_median,
         trans_rate = trans_num/age,
         col = case_when(
           trans_num >= sim_0.05 & trans_num <= sim_0.95 ~ "mid",
           trans_num > sim_0.95 ~ "big",
           trans_num < sim_0.05 ~ "small"),
         famcol = ifelse(family %in% rab_fams, "#F8766D", "grey40"),
         Family = fct_reorder(family, n_diff))
saveRDS(ta_results, here(paste0("data/prog_data/ta_results_", suffix, ".Rdata")))
} else {
  ta_results <- readRDS(here(paste0("data/prog_data/ta_results_", suffix, ".Rdata")))
}
