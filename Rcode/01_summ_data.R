# TABLE S3 --------------------------------
# getting taxonomic information for every species from fishbase
all_fishbase <- load_taxa() %>%
  as_tibble() %>%
  clean_names() %>%
  group_by(family) %>%
  add_count(name = "total_sp") %>%
  select(family, total_sp) %>%
  unique() 



tableS1 <- trim_data %>%
  select(Species = sp, Family = family, lat_centroid, tridepth, nsp) %>%
  group_by(Family) %>%
  add_count(tridepth, name = "ncat") %>%
  mutate(med_lat = median(lat_centroid),
         ncat = (ncat/nsp)*100) %>%
  select(-Species, -lat_centroid) %>%
  unique() %>%
  pivot_wider(names_from = tridepth, values_from = ncat, values_fill = 0) %>%
  arrange(Family) %>%
  left_join(all_fishbase, by = c("Family" = "family")) %>%
  mutate(perc_cover = (nsp/total_sp) *100) %>%
  select(Family, nsp, perc_cover, med_lat, shallow, intermediate, deep) %>%
  mutate(across(perc_cover:deep, round, 1)) %>%
  rename(n_species = nsp, "% cover" = perc_cover, "median lat" = med_lat) 

kable(tableS1, booktab = TRUE, format = "latex", #longtable = TRUE, #escape = FALSE, 
      align = c("l", rep("c", 7)), linesep = "", digits = 1) %>%
 # kable_styling(font_size = 10, latex_options = "repeat_header") %>%
  save_kable(here("results/tables/TableS3.pdf"))



# FAMILY DEPTH SUMMARY SUPP FIG --------------------------------

# Figure summarizing depths for each family 
trim_data %>%
  select(species = sp, family, depth_range_shallow, depth_range_deep,
         tridepth, med_depth, nsp) %>%
  ggplot() +
  geom_boxplot(aes(x = family, y = med_depth), col = "grey40") +
  theme_classic() +
  coord_flip() +
  xlab("") + ylab("Median Depth by Species (m)")
ggsave(here("results/figures/", "depth_summ_fig.pdf"))



# RAW DATA SUPP. TABLE --------------------------------

# raw data supplement table
fb_supp_table <- alldata %>%
  group_by(family) %>%
  add_count(name = "nsp") %>%
  select(species = sp, family, freshwater, marine,
         depth_min = depth_range_shallow,
         depth_max = depth_range_deep, depth_cat = tridepth) %>%
  arrange(family, species)
write_csv(fb_supp_table, here("results/tables/", "fishbase_data.csv"))



# random sample to validate fishbase depth estimates
# fb_supp_table %>%
#   ungroup() %>%
#   slice_sample(n = 100, replace = FALSE) %>%
#   write.csv("results/tables/validate_depths.csv")
