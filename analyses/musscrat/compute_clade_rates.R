library(RevGadgets)
library(phytools)
library(tidyverse)

computeCladeRates = function(tree, tips, samples, parameter, inverse=FALSE) {

    # match R and Rev indexes
    map = RevGadgets:::matchNodes(tree)
    map = map[-(tree$Nnode + 2),]

    # get the MRCA and all the descendant branches of the MRCA
    mrca  = ape::getMRCA(tree, tips)
    descs = phytools::getDescendants(tree, mrca)

    # retrieve the samples
    if ( inverse == FALSE ) {
        # get the edges for this clade
        rev_nodes = map$Rev[map$R %in% descs]
    } else {
        # get the edges outside this clade
        rev_nodes = map$Rev[map$R %in% descs == FALSE]
    }

    param_names = paste0(parameter,"[",rev_nodes,"]")
    return(samples[,param_names])

}


# define the clade (with two species that span the MRCA)
dir <- "~/Documents/Rabsoky_proj/"
source(paste0(dir, "dataprep.R")) # loading prepped data

stopifnot(all(data$sp == tree$tip.label))

clades <- data %>%
    select(sp, family) %>%
    group_by(family) %>%
    add_count(name = "nsp") %>%
    filter(nsp > 20) %>%
    filter(!family %in% c("Cynoglossidae", "Congridae", "Engraulidae")) %>%  #removing paraphyletic families   
    nest() 


# read the samples
samples = read.table(paste0(dir,"musscrat/output/Edrun1_output_UCLN_state_dependent_1params.log"),
                     header=TRUE, stringsAsFactors=FALSE, sep="\t", check.names=FALSE)

# compute the rates (options: state_branch_rate, background_branch_rate, etc)
clade_rates <- list()
background_rates <- list()
for(i in 1:nrow(clades)){
    cat(paste0("Running family: ", clades$family[[i]], "\n"))
    clade_rates[[i]] = computeCladeRates(tree, clades$data[[i]]$sp, 
                                           samples, "state_branch_rate", inverse = FALSE)
    background_rates[[i]] = computeCladeRates(tree, clades$data[[i]]$sp, 
                                           samples, "background_branch_rate", inverse = FALSE)
    cat("Finished\n")
}
names(clade_rates) <- names(background_rates) <- clades$family

# out <- list(clade_rates = clade_rates, background_rates = background_rates)
# saveRDS(out, paste0(dir, "musscrat/output/Edrun1_rates.Rdata"))


mean_rates <- lapply(clade_rates, function(x) mean(unlist(x)))
rates_df <- enframe(unlist(mean_rates), name = "Family", value = "depth_rate")

all_res <- readRDS(paste0(dir, "prog_data/data_fig1.Rdata")) %>%
    left_join(rates_df, by = "Family") 

source("~/Documents/Rfunctions/collapse_tree.R")
tax_df <- all_res %>%
    select(Genus_species = Species, Family) %>%
    unnest()
fam_tree <- collapse_tree(
    drop.tip(tree, setdiff(tree$tip.label, tax_df$Genus_species)),
    tax_df, clade_level = "Family")

all_res <- all_res %>%
    arrange(match(Family, fam_tree$tip.label)) 

stopifnot(all(all_res$Family == fam_tree$tip.label))

library(ggrepel)
ggplot(all_res, aes(x = med_lambda, y = depth_rate)) +
    geom_point(aes(col = famcol), alpha = 0.8, cex = 2.5) +
    theme_classic() +
    ylab("Rate of Depth Evolution") +
    xlab("Speciation Rate") +
    theme(legend.position = "none",
          axis.line = element_line(size = 0.2),
          axis.ticks = element_line(size = 0.1),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) +
    geom_label_repel(data = filter(all_res, med_lambda > 0.2 | depth_rate > 0.3), 
                     aes(label = Family, color = famcol), 
                     cex = 5, label.size = NA, fill = NA) +
    stat_smooth(method = "lm", se = FALSE) +
    scale_color_manual(values = c("#F8766D", "grey40"))


ggplot(all_res, aes(x = abs(med_lat), y = depth_rate)) +
    geom_point(aes(col = famcol), alpha = 0.8, cex = 3) +
    theme_classic() +
    ylab("Rate of Depth Evolution") +
    xlab("Latitude") +
    theme(legend.position = "none",
          axis.line = element_line(size = 0.2),
          axis.ticks = element_line(size = 0.1),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12)) +
    geom_label_repel(data = filter(all_res, med_lat > 40 | depth_rate > 0.3), 
                     aes(label = Family, color = famcol), 
                     cex = 5, label.size = NA, fill = NA) +
    scale_color_manual(values = c("#F8766D", "grey40"))


library(geomorph)
lm <- procD.lm(depth_rate ~ med_lambda, data = all_res, iter = 999)
summary(lm)
pgls <- procD.pgls(depth_rate ~ med_lambda, data = all_res, phy = fam_tree, iter = 999)
summary(pgls)



lm2 <- procD.lm(depth_rate ~ abs(med_lat), data = all_res, iter = 999)
summary(lm2)
pgls2 <- procD.pgls(depth_rate ~ abs(med_lat), data = all_res, phy = fam_tree, iter = 999)
summary(pgls2)





# compute mean
i <- 10
mean( unlist(clade_rates[[i]]) )
mean( unlist(background_rates[[i]]) )

# make a histogram of the mean rates per lineage?
hist( colMeans( clade_rates[[i]] ), breaks=20 )

# compare inside and outside rates
plot(density(unlist(clade_rates[[i]])), col="red")
lines(density(unlist(background_rates[[i]])), col="blue")
