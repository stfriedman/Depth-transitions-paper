# setting group colors/aesthetics
  fam_ord <- ta_results %>%
    filter(perc == 0.9) %>%
    ungroup() %>%
    mutate(
      family = as.factor(family),
      family = reorder(family, n_diff)
    )

  fig1_data <- ta_results %>%
    mutate(perc_title = paste0(perc*100, "% sampled")) 
  fig1_data$family <- factor(fig1_data$family, levels=levels(fam_ord$family))
  
  famcol <- ifelse(levels(fig1_data$family) %in% rab_fams, "#F8766D", "grey40")
  pal <- wes_palette("Zissou1", 100, type = "continuous")


#pdf(here(paste0("results/figures/Fig1_", suffix, ".pdf")), width = 12, height = 10)

  # figure of transition by clade compared to null
  ggplot(fig1_data, aes(x = family, y = trans_num)) +
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
    facet_wrap(~perc_title) + 
    theme_classic() +
    theme(
      #axis.text.y = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      axis.text.y = element_text(colour = famcol, size = 8),
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

# 
#   # figure of transitions vs. expectation and latitude
#   p3 <- ggplot(fig1_data, aes(y = n_diff, x = abs(med_lat))) +
#     theme_classic() +
#     geom_hline(yintercept = 0, lty = "dashed", col = "grey70", cex = 0.3) +
#     geom_point(aes(fill = med_lambda), pch = 21, size = 4, stroke = 0) +
#     scale_fill_gradientn(colors = pal) +
#     scale_color_manual(values = c("#F8766D", "grey40"), guide = "none") +
#     geom_label_repel(
#       data = filter(ta_results, n_diff > 10 | n_diff < -20 | med_lat < -20 |
#         med_lat > 30 | med_lambda > 0.3),
#       aes(label = family, color = famcol),
#       cex = 3.2, label.size = NA, fill = NA
#     ) +
#     labs(fill = "Speciation\nRate") +
#     theme(
#       legend.position = c(0.82, 0.15),
#       legend.background = element_rect(size = 0.2, colour = "grey70"),
#       axis.line = element_line(size = 0.2),
#       axis.ticks = element_line(size = 0.1),
#       axis.title = element_text(size = 12)
#     ) +
#     xlab("Median Latitude") +
#     ylab("Observed transitions - expected")

    dev.off()
