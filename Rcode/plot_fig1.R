if (suffix == "redepth") {
  # reading in results from normal depth categorization
  ta_results_revisions <- readRDS(here("data/prog_data/ta_results_revisions.Rdata")) %>%
    mutate(cat = "norm") %>%
    select(mrca, cat, everything()) %>%
    mutate(family = as.factor(family),
           family = reorder(family, -n_diff))
  
  # combining with results from depth recategorizations
  fig1_data <- bind_rows(ta_results, ta_results_revisions) %>%
    ungroup() %>%
    mutate(
      family = as.factor(family),
      family = fct_relevel(family, rev(levels(ta_results_revisions$family)))
    ) %>%
    mutate(cat = case_when(
      cat == "small" ~ "shallow: 0-100m",
      cat == "norm" ~ "shallow: 0-200m",
      cat == "big" ~ "shallow: 0-300m"
    ))
  
  # setting group colors/aesthetics
  famcol <- ifelse(levels(fig1_data$family) %in% rab_fams, "#F8766D", "grey40")
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  depth_cols <- setNames(
    c("powderblue", "#2C8EB5", "#16465B"),
    c("shallow", "intermediate", "deep")
  )
  

  ggplot(fig1_data, aes(x = family, y = trans_num, group = cat)) +
    geom_point(aes(fill = col),
      pch = 21, cex = 2, stroke = 0,
      alpha = 0.8, position = position_dodge(0.7)
    ) +
    geom_linerange(aes(
      x = family, ymin = sim_0.05,
      ymax = sim_0.95, col = cat
    ),
    lwd = 1, alpha = 0.8,
    position = position_dodge(0.7)
    ) +
    theme_classic() +
    coord_flip() +
    theme(
      axis.text.y = element_text(colour = famcol),
      plot.margin = unit(c(0, 0.4, 0, 0), "cm"),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.line = element_line(size = 0.2),
      axis.ticks.x = element_line(size = 0.1),
      axis.ticks.y = element_blank(),
      axis.title = element_text(size = 12),
      axis.title.y = element_blank()
    ) +
    labs(colour = "Depth categories") +
    ylab("Number of Transitions") +
    scale_color_manual(values = c("black", "grey40", "grey70")) +
    scale_fill_manual(values = c(
      wes_palette("Darjeeling1")[2], "grey40",
      wes_palette("Darjeeling1")[3]
    ), guide = "none")

  ggsave(here("results/figures/Fig1_redepth.pdf"))
  
} else {

  # setting group colors/aesthetics
  fig1_data <- ta_results %>%
    mutate(
      family = as.factor(family),
      family = reorder(family, n_diff)
    )
  
  famcol <- ifelse(levels(fig1_data$family) %in% rab_fams, "#F8766D", "grey40")
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  depth_cols <- setNames(
    c("powderblue", "#2C8EB5", "#16465B"),
    c("shallow", "intermediate", "deep")
  )
  
  
  # barplot of depth occupation for each clade
  p1 <- fig1_data %>%
    select(family, depth) %>%
    unnest(cols = c(depth)) %>%
    ggplot(aes(fill = tridepth, y = ncat, x = family)) +
    geom_bar(position = "fill", stat = "identity", width = 0.7) +
    coord_flip() +
    theme_classic() +
    theme(
      axis.text.y = element_text(colour = famcol),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.position = "none"
    ) +
    xlab("") +
    ylab("") +
    scale_fill_manual(values = depth_cols)


  # figure of tranisiton by clade compared to null
  p2 <- ggplot(fig1_data, aes(x = family, y = trans_num)) +
    geom_segment(aes(x = family, xend = family, y = sim_0.05, yend = sim_0.95),
      col = "grey90", lwd = 3, lineend = "round"
    ) +
    geom_segment(aes(x = family, xend = family, y = sim_median, yend = trans_num),
      col = "grey50", lwd = 0.6
    ) +
    geom_point(aes(x = family, y = sim_median),
      pch = 21,
      fill = "white", col = "grey50", size = 3, stroke = 1
    ) +
    geom_point(aes(col = col), size = 3.3) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      plot.margin = unit(c(0, 0.4, 0, 0), "cm"),
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      panel.grid.major.y = element_line(size = 0.4),
      axis.line.x = element_line(size = 0.2),
      axis.line.y = element_blank(),
      axis.ticks.x = element_line(size = 0.1),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      axis.title = element_text(size = 12),
      axis.title.y = element_blank()
    ) +
    coord_flip() +
    ylab("Number of Transitions") +
    labs(col = "Speciation Rate") +
    scale_color_manual(values = c(
      wes_palette("Darjeeling1")[2], "grey40",
      wes_palette("Darjeeling1")[3]
    ))


  # figure of transitions vs. expectation and latitude
  p3 <- ggplot(fig1_data, aes(y = n_diff, x = abs(med_lat))) +
    theme_classic() +
    geom_hline(yintercept = 0, lty = "dashed", col = "grey70", cex = 0.3) +
    geom_point(aes(fill = med_lambda), pch = 21, size = 4, stroke = 0) +
    scale_fill_gradientn(colors = pal) +
    scale_color_manual(values = c("#F8766D", "grey40"), guide = "none") +
    geom_label_repel(
      data = filter(ta_results, n_diff > 10 | n_diff < -20 | med_lat < -20 |
        med_lat > 30 | med_lambda > 0.3),
      aes(label = family, color = famcol),
      cex = 3.2, label.size = NA, fill = NA
    ) +
    labs(fill = "Speciation\nRate") +
    theme(
      legend.position = c(0.82, 0.15),
      legend.background = element_rect(size = 0.2, colour = "grey70"),
      axis.line = element_line(size = 0.2),
      axis.ticks = element_line(size = 0.1),
      axis.title = element_text(size = 12)
    ) +
    xlab("Median Latitude") +
    ylab("Observed transitions - expected")

  
  pdf(here(paste0("results/figures/Fig1_", suffix, ".pdf")), width = 10, height = 8)
  plot_grid(p1, p2, p3, align = "h", nrow = 1, rel_widths = c(0.4, 0.8, 1.2))
  dev.off()
}
