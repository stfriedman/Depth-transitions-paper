# plot rates for figure 3
if(!exists("ta_results")){
  ta_results <- readRDS(here("data/prog_data/ta_results_revisions.Rdata"))
}

param_run1 <- read_csv(here("analyses/musscrat/output/params_run1_prior500.csv")) %>%
  rename(tropical = `zeta[1]`, temperate = `zeta[2]`, polar = `zeta[3]`)

out <- readRDS(here("analyses/musscrat/output/Edrun1_rates.Rdata"))
clade_rates <- out$clade_rates
mean_rates <- lapply(clade_rates, function(x) mean(unlist(x)))
rates_df <- enframe(unlist(mean_rates), name = "family", value = "depth_rate")

res <- ta_results %>%
  left_join(rates_df, by = "family") %>%
  mutate(lat_cat = case_when(
    med_lat >-30 & med_lat < 30 ~ "tropical",
    abs(med_lat) > 30 & abs(med_lat) < 60 ~ "temperate",
    abs(med_lat) > 60 ~ "polar"
  )) %>%
  mutate(focal_clade = case_when(
    family == "Notothenioid" ~ "Notothenioid",
    family == "Sebastidae" ~ "Sebastidae",
    family == "Liparidae" ~ "Liparidae",
    family == "Zoarcidae" ~ "Zoarcidae",
    family == "Pleuronectidae" ~ "Pleuronectidae",
    TRUE ~ "none"
  ))



p1 <- param_run1 %>%
  pivot_longer(tropical:polar, names_to = "zeta", values_to = "rate") %>%
  ggplot(aes(x = rate, y = ..scaled..)) +
  geom_density(aes(fill = zeta), alpha = 0.8, col = NA) +
  theme_classic() +
  theme(legend.position = "none",
        aspect.ratio=1,
        legend.title = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12)) +
  scale_fill_manual(values = c("#472D7B", "#20A486", "#E6AE48FF")) +
  xlab("Rate of Depth Evolution") + ylab("Density")

p2 <- ggplot(res, aes(x = depth_rate, y = med_lambda)) +
  geom_point(aes(col = focal_clade), alpha = 0.8, cex = 2.4) +
  theme_classic() +
  xlab("") +
  ylab("Speciation Rate") +
  theme(legend.position = "none", 
        aspect.ratio=1,
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12)) +
  scale_color_manual(values = c("#4B8460", "grey60", "#804D8A", 
                                "#203F60", "#E6AF4B", "#AB4932")) 
# geom_label_repel(data = filter(res, med_lambda > 0.2 | depth_rate > 0.3),
#                  aes(label = Family),
#                  cex = 3, label.size = NA, fill = NA) 


p4 <- ggplot(res, aes(x = depth_rate, y = abs(med_lat))) +
  geom_point(aes(col = focal_clade), alpha = 0.8, cex = 2.4) +
  theme_classic() +
  xlab("Rate of Depth Evolution") +
  ylab("Latitude") +
  theme(legend.position = "none", 
        aspect.ratio=1,
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.1),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 12)) +
  scale_color_manual(values = c("#4B8460", "grey60", "#804D8A", 
                                "#203F60", "#E6AF4B", "#AB4932")) 


if(!exists(p3)){
  source(here("Rcode/04_rates.R"))
}


a <- plot_grid(p1, plot_grid(p2, p4, nrow = 2),
               ncol = 2, rel_widths = c(1, 0.8), align = "v")

pdf(here("results/figures/fig3_v2.pdf"))
plot_grid(
  a, 
  p3, 
  nrow = 2, align = "v", axis = "r")
dev.off()



# param_run1 %>% 
#   select(tropical:polar) %>%
#   pivot_longer(tropical:polar, names_to = "zeta", values_to = "rate") %>%
#   group_by(zeta) %>%
#   summarize(mean = mean(rate))
# 
# mean(param_run1$is_state_dependent)
# mean(param_run1$total_num_changes)


# # number of transition
# param_run1 %>%
#   ggplot(aes(x = total_num_changes, y = ..scaled..)) +
#   geom_density(alpha = 0.8, fill = "grey80", col = "grey70", lwd = 1) +
#   theme_classic() +
#   theme(legend.position = "none",
#         legend.title = element_blank(),
#         axis.line = element_line(size = 0.2),
#         axis.ticks = element_line(size = 0.1),
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 14)) +
#   xlab("Number of Transitions") + ylab("Density")

