# LOAD LIBRARIES ---------------------------------------------------------
library(circlize)
library(chorddiag)  #devtools::install_github("mattflor/chorddiag")


# EVALUATING MK MODELS -----------------------------------

# setting up lat/depth categories
grps <- as.factor(setNames(paste0(alldata$lat_cat, "_", alldata$tridepth), alldata$sp))

stopifnot(all(tree$tip.label == names(grps)))



if(rerun){
  fit_tib <- tibble(model_type = c("ER", "SYM", "ARD"),
                    model = list("ER", "SYM", "ARD")) %>%
    mutate(fit = map(model, ~fitMk(tree, grps, model = .x)),
           AIC = map_dbl(fit, AIC))
  saveRDS(fit_tib, here(paste0("data/prog_data/fit_tib_", suffix, ".Rdata")))
} else {
  fit_tib <- readRDS(here(paste0("data/prog_data/fit_tib_", suffix, ".Rdata")))
}

bestfit <- fit_tib$fit[which.min(fit_tib$AIC)][[1]]


# extracting and reformatting Q matrix for long format
Q <- data.frame(Q.mat(bestfit)) %>%
  rownames_to_column("to") %>%
  gather(key = "from", value = "value", -to) %>%
  filter(from != to) %>%
  arrange(from) %>%
  filter(value != 0)


# Q matrix formatting for supplemental materials
supp_q <- matrix(NA,length(bestfit$states),length(bestfit$states))
supp_q[] <- c(0,bestfit$rates)[bestfit$index.matrix+1]
diag(supp_q)<-0
diag(supp_q)<--rowSums(supp_q)
colnames(supp_q)<-rownames(supp_q)<-gsub("_", " ", bestfit$states)

kable(supp_q, booktab = TRUE, format = "latex", escape = FALSE,
      align = rep("c", 9), linesep = "", digits = 3) %>%
  kable_styling(font_size = 10, latex_options="scale_down") %>%
  save_kable("results/tables/TableS4_Qmat.pdf")



# AESTHETICS FOR PLOT -----------------------------------------
# grouping by latitude
nm <- unique(Q$from)
group <- structure(gsub("_[a-z]+", "", nm), names = nm)


# cols for grids
depthcol <- structure(gsub("[a-z]+_", "", nm), names = nm)
depthcol[grepl("shallow", depthcol)] <- "#B0E0E6"
depthcol[grepl("intermediate", depthcol)] <- "#2C8EB5"
depthcol[grepl("deep", depthcol)] <- "#16465B"


## setting colors for arrows
depthcol2 <- setNames(Q$to, Q$to)
pol <- which(grepl("polar", names(depthcol2)))
temp <- which(grepl("temperate", names(depthcol2)))
trop <- which(grepl("tropical", names(depthcol2)))

depthcol2[grepl("shallow", depthcol2)] <- "#b0e0e6"
depthcol2[grepl("intermediate", depthcol2)] <- "#2C8EB5"
depthcol2[grepl("deep", depthcol2)] <- "#16465B"


# loop to color only one latitude at a time
ar_cols <- rep(list(depthcol2), 3)
n <- length(depthcol2)
id <- list(pol, temp, trop)
x <- combn(1:3, 2)
for(i in 1:3){
  int <- unlist(id[x[,i]])
  fo <- seq_len(n)[!seq_len(n) %in% int] 
  ar_cols[[i]][int] <- sapply(ar_cols[[i]][int], function(x) paste0(x, "0d")) 
  ar_cols[[i]][fo] <- sapply(ar_cols[[i]][fo], function(x) paste0(x, "d9")) 
}




# PLOTTING FIGURE 2 -----------------------------------------
pdf(here(paste0("results/figures/fig2_", suffix, ".pdf")), width=12, height=8.5)

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, points.overflow.warning = FALSE)
par(mar = rep(0, 4), mfrow = c(1,3))

for(i in 1:3){
# Base plot
chordDiagram(
  x = Q, 
  grid.col = depthcol,
  col = ar_cols[[i]],
  # grid.col = latcol,
  # col = depthcol,
  group = group,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.08, 0.1),
  link.arr.type = "big.arrow", 
  # link.sort = TRUE, 
  # link.decreasing = TRUE,
  link.largest.ontop = TRUE,
  preAllocateTracks = list(
    track.height = 0.1,
    track.margin = c(0.01, 0)
  ))

# Add text and axis
circos.trackPlotRegion(
  track.index = 2, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    #sector.index = get.cell.meta.data("sector.index")
    sector.index = gsub("[a-z]+_", "", get.cell.meta.data("sector.index"))
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 0.5, 
      col = "white",
      labels = sector.index, 
      facing = "bending", 
      cex = 1,
      niceFacing = TRUE
    )
  }
)


highlight.sector(names(group[7:9]), track.index = 1, facing = "bending", font = 2,
                 col = ifelse(i == 1, "#E6AE48FF", "#E6AE483a"),
                 border = ifelse(i == 1, TRUE, FALSE), 
                 lwd = ifelse(i == 1, 2, 0.01),
                 text = "tropical", cex = 1.5, text.col = "white", niceFacing = TRUE)
highlight.sector(names(group[4:6]), track.index = 1, facing = "bending", font = 2,
                 border = ifelse(i == 2, TRUE, FALSE),
                 lwd = ifelse(i == 2, 2, 0.01),
                 col = ifelse(i == 2, "#20A486", "#20A4863a"),
                 text = "temperate", cex = 1.5, text.col = "white", niceFacing = TRUE)
highlight.sector(names(group[1:3]), track.index = 1, facing = "bending", font = 2,
                 col = ifelse(i == 3, "#472D7B", "#472D7B3a"),
                 border = ifelse(i == 3, TRUE, FALSE),
                 lwd = ifelse(i == 3, 2, 0.01),
                 text = "polar", cex = 1.5, text.col = "white", niceFacing = TRUE)
}
dev.off()
