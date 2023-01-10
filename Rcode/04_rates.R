# RATE OF DEPTH EVOLUTION BY CLADE ----------------------
med_depth <- trim_data %>% 
  select(sp, family, med_depth) %>%
  mutate(dumm = 1) %>%
  column_to_rownames("sp") 
trim_fam_vect <- setNames(med_depth$family, rownames(med_depth))
depth_rates <- compare.evol.rates(med_depth[,-1], gp = trim_fam_vect, 
                                  phy = trim_tree, iter = 999)

rate_results <- tibble(family = names(depth_rates$sigma.d.gp), 
                  rate = depth_rates$sigma.d.gp)



# DIFFERENT SPECIATION RATE METRICS ----------------------
species_rate <- ta_results %>% 
  mutate(famtree = map(species, ~drop.tip(tree, setdiff(tree$tip.label, .)))) %>%
  mutate(tip_nd = map(famtree, nodeDensity),
         nd = map_dbl(tip_nd, mean),
         tip_dr = map(famtree, DRstat),
         dr = map_dbl(tip_dr, mean)) %>%
 # mutate(focal_clade = family %in% rab_fams)
  mutate(focal_clade = case_when(
    family == "Notothenioid" ~ "Notothenioid",
    family == "Sebastidae" ~ "Sebastidae",
    family == "Liparidae" ~ "Liparidae",
    family == "Zoarcidae" ~ "Zoarcidae",
    family == "Pleuronectidae" ~ "Pleuronectidae",
    TRUE ~ "none"
  ))



# QUANTIFYING DIFFERENT RATE METHODS ----------------------

summary(lm(dr ~ med_lambda, data = species_rate))
summary(lm(nd ~ med_lambda, data = species_rate))


summary(lm(dr ~ med_lat, data = species_rate))
summary(lm(nd ~ med_lat, data = species_rate))
summary(lm(med_lambda ~ med_lat, data = species_rate))

# all are very significant!



# PLOTTING -------------------------- 
labs <- c(
  dr = "DR statistic",
  med_lambda = "BAMM",
  nd = "Node Density"
)


p3 <- species_rate %>%
  pivot_longer(c(med_lambda, nd, dr), names_to = "method", values_to = "rate") %>%
  ggplot(aes(y = rate, x = trans_rate)) +
  facet_wrap(~method, scales = "free_y", strip.position = "left", 
             labeller = labeller(method = labs)) +
  theme_classic() +
  geom_point(aes(col = focal_clade), alpha = 0.8, cex = 2.2) +
  xlab("Rate of Depth Transitions") + ylab("") +
  theme(aspect.ratio = 1,
        legend.position = "none",
        strip.placement = "outside",
        strip.background = element_blank(),
        axis.title = element_text(size = 12),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 8),
        strip.text = element_text(size = 12)) +
  scale_color_manual(values = c("#4B8460", "grey60", "#804D8A", 
                                "#203F60", "#E6AF4B", "#AB4932")) 
#stat_smooth(method = "lm", se = FALSE, col = "grey40", cex = 0.5)

