# All helper functions needed to run analyses/recreate plots for depth paper
PKG <- c(
  "picante",
  "RandomFieldsUtils",
  "rfishbase",
  "tidyverse",
  "ape",
  "viridis",
  'plotrix'
)

PKG <- unique(PKG)
for (p in PKG) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    require(p, character.only = TRUE)
  }
}
rm(PKG, p)


optim_pcd <- function(comm, tree, PSVmncd = NULL, PSVpool = NULL, reps = 10^4) {
  SSii <- PSVmncd
  SCii <- PSVpool
  comm[comm > 0] <- 1
  if (is(tree)[1] == "phylo") {
    if (is.null(tree$edge.length)) {
      tree <- compute.brlen(tree, 1)
    }
    tree <- prune.sample(comm, tree)
    V <- vcv.phylo(tree, corr = TRUE)
    comm <- comm[, tree$tip.label]
  } else {
    V <- tree
    species <- colnames(comm)
    preval <- colSums(comm) / sum(comm)
    species <- species[preval > 0]
    V <- V[species, species]
    comm <- comm[, colnames(V)]
  }
  if (!is.null(SSii) & length(SSii) != max(rowSums(comm))) {
    stop("The length of PSVmncd is not equal to the species richness of the community with the highest species richness. Set PSVmncd=NULL, PSVpool=NULL and run pcd again.")
  }
  m <- dim(comm)[1]
  n <- dim(comm)[2]
  nsr <- max(rowSums(comm))
  if (is.null(SSii) & is.null(SCii)) {
    SSii <- array(0, nsr)
    n1 <- 2
    for (n2 in 1:nsr) {
      temp <- array(0, reps)
      for (t in 1:reps) {
        # rp <- sample(n)
        # pick1 <- rp[1:n1]
        # rp <- sample(n)
        # pick2 <- rp[1:n2]

        pick1 <- sample(n, size = n1)
        pick2 <- sample(n, size = n2)

        C11 <- V[pick1, pick1]
        C22 <- V[pick2, pick2]
        C12 <- V[pick1, pick2]
        invC22 <- solvex(C22)
        S11 <- C11 - C12 %*% invC22 %*% t(C12)
        SS11 <- (n1 * sum(diag(S11)) - sum(S11)) / (n1 * (n1 - 1))
        temp[t] <- SS11
      }
      SSii[n2] <- mean(temp)
    }
    SCii <- 1 - (sum(V) - sum(diag(V))) / (n * (n - 1))
  }
  PCD <- array(NA, c(m, m))
  PCDc <- array(NA, c(m, m))
  PCDp <- array(NA, c(m, m))
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      pick1 <- (1:n)[comm[i, ] == 1]
      pick2 <- (1:n)[comm[j, ] == 1]
      n1 <- length(pick1)
      n2 <- length(pick2)
      C <- V[c(pick1, pick2), c(pick1, pick2)]
      C11 <- C[1:n1, 1:n1]
      C22 <- C[(n1 + 1):(n1 + n2), (n1 + 1):(n1 + n2)]
      C12 <- C[1:n1, (n1 + 1):(n1 + n2)]
      if (is.null(dim(C12))) {
        if (is.null(dim(C22))) {
          C12 <- as.matrix(C12)
        } else {
          C12 <- t(as.matrix(C12))
        }
      }
      invC11 <- solvex(C11)
      S22 <- C22 - t(C12) %*% invC11 %*% C12
      invC22 <- solvex(C22)
      S11 <- C11 - C12 %*% invC22 %*% t(C12)
      if (n1 > 1) {
        SC11 <- (n1 * sum(diag(C11)) - sum(C11)) / (n1 *
          (n1 - 1))
        SS11 <- (n1 * sum(diag(S11)) - sum(S11)) / (n1 *
          (n1 - 1))
      } else {
        SC11 <- (n1 * sum(diag(C11)) - sum(C11)) / (n1 *
          (n1))
        SS11 <- (n1 * sum(diag(S11)) - sum(S11)) / (n1 *
          (n1))
      }
      if (n2 > 1) {
        SC22 <- (n2 * sum(diag(C22)) - sum(C22)) / (n2 *
          (n2 - 1))
        SS22 <- (n2 * sum(diag(S22)) - sum(S22)) / (n2 *
          (n2 - 1))
      } else {
        SC22 <- (n2 * sum(diag(C22)) - sum(C22)) / (n2 *
          (n2))
        SS22 <- (n2 * sum(diag(S22)) - sum(S22)) / (n2 *
          (n2))
      }
      D <- (n1 * SS11 + n2 * SS22) / (n1 * SC11 + n2 * SC22)
      a <- length(unique(c(pick1, pick2)))
      b <- length(pick1) - a
      cc <- length(pick2) - a
      dsor <- 2 * a / (2 * a + b + cc) - 1
      pred.D <- (n1 * SSii[n2] + n2 * SSii[n1]) / (n1 * SCii +
        n2 * SCii)
      pred.dsor <- 1 - 2 * n1 * n2 / ((n1 + n2) * n)
      PCD[i, j] <- D / pred.D
      PCDc[i, j] <- dsor / pred.dsor
      PCDp[i, j] <- PCD[i, j] / PCDc[i, j]
    }
  }
  colnames(PCD) <- rownames(comm)
  rownames(PCD) <- rownames(comm)
  colnames(PCDc) <- rownames(comm)
  rownames(PCDc) <- rownames(comm)
  colnames(PCDp) <- rownames(comm)
  rownames(PCDp) <- rownames(comm)
  return(list(
    PCD = as.dist(t(PCD)), PCDc = as.dist(t(PCDc)),
    PCDp = as.dist(t(PCDp)), PSVmncd = SSii, PSVpool = SCii
  ))
}




# ratio of number of speciation events / age of clade; code from Title and Rabosky 2019
nodeDensity <- function(tree) {
  maxBT <- max(branching.times(tree))
  nodeCounts <- sapply(1:length(tree$tip.label), function(x) {
    n <- 0
    childnode <- x
    parentNode <- tree$edge[which(tree$edge[, 2] == childnode), 1]
    while (parentNode != (length(tree$tip.label) + 1)) {
      n <- n + 1
      childnode <- parentNode
      parentNode <- tree$edge[which(tree$edge[, 2] == childnode), 1]
    }
    return(n)
  })

  # avoid rates of 0 by adding the root node
  nodeCounts <- nodeCounts + 1

  return(setNames(nodeCounts / maxBT, tree$tip.label))
}


# DR metric / inverse equal splits; code from Title and Rabosky 2019
DRstat <- function(tree) {
  spRate <- function(sp, tree) {
    # get branch lengths from root to tip
    edges <- vector()
    daughterNode <- match(sp, tree$tip.label)
    while (daughterNode != (length(tree$tip.label) + 1)) {
      parentNode <- tree$edge[which(tree$edge[, 2] == daughterNode), 1]
      edges <- c(edges, tree$edge.length[which(tree$edge[, 1] == parentNode & tree$edge[, 2] == daughterNode)])
      daughterNode <- parentNode
    }

    res <- sum(sapply(1:length(edges), function(x) edges[x] * (1 / (2^(x - 1)))))
    res <- res^(-1)

    return(res)
  }

  rates <- unlist(lapply(tree$tip.label, function(x) spRate(x, tree)))
  names(rates) <- tree$tip.label

  return(rates)
}


get_mapped_edge <- function(x) {
  data.frame(x$mapped.edge) %>%
    rownames_to_column("edges") %>%
    as_tibble() %>%
    separate(edges, c("node_1", "node_2"), sep = ",")
}



num_trans <- function(mrca, simtree) {
  trans_count <- function(desc) {
    edges <- get_mapped_edge(simtree) %>%
      mutate_at(vars(node_1:node_2), as.numeric)

    edges %>%
      filter(node_1 %in% desc) %>%
      pivot_longer(-c(node_1, node_2), names_to = "state", values_to = "value") %>%
      filter(value != 0) %>%
      group_by(node_1, node_2) %>%
      add_count() %>%
      filter(n > 1) %>%
      select(node_1, node_2) %>%
      unique() %>%
      nrow(.)
  }
  tibble(mrca) %>%
    mutate(
      desc = map(mrca, ~ c(.x, getDescendants(simtree, .x))),
      trans_num = map_dbl(desc, trans_count)
    )
}



gg_color <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



color_clades <- function(tree, nodes, col.ordered = TRUE) {
  n <- length(nodes)

  if (col.ordered == TRUE) {
    cols <- viridis(n)
  } else {
    gg <- gg_color(n)
    cols <- sample(gg)
  }

  br_col <- rep("grey85", nrow(tree$edge))

  for (i in 1:n) {
    clade_desc <- getDescendants(tree, nodes[[i]])

    if (all(!clade_desc %in% nodes)) {
      br_col[which.edge(tree, clade_desc)] <- cols[[i]]
    } else {
      overlap_clades <- clade_desc[which(clade_desc %in% nodes)]

      if (length(overlap_clades) > 1) {
        clade_edges <- c()
        for (j in 1:length(overlap_clades)) {
          clade_edges <- append(
            clade_edges,
            which.edge(tree, getDescendants(tree, overlap_clades[[j]]))
          )
        }
      } else {
        clade_edges <- which.edge(tree, getDescendants(tree, overlap_clades))
      }

      overlap_edges <- which.edge(tree, clade_desc)
      br_col[overlap_edges[!overlap_edges %in% clade_edges]] <- cols[[i]]
    }
  }

  clade_cols <- br_col
  return(clade_cols)
}



fishbase_synonyms <- function(species) {
  tax_fb <- load_taxa() %>%
    as_tibble() %>%
    dplyr::select(Species, Family, Order, SpecCode) %>%
    dplyr::mutate(Genus_species = str_replace(Species, " ", "_"))

  if (class(species) == "phylo") {
    x <- tibble(Genus_species = species$tip.label) %>%
      left_join(., tax_fb, by = "Genus_species") %>%
      mutate(fishbase = ifelse(!is.na(Family), Genus_species, NA)) %>%
      select(Genus_species, fishbase, Family, Order, SpecCode)
  } else {
    x <- tibble(Genus_species = species) %>%
      # mutate(Genus_species = gsub("_", " ", Genus_species)) %>%
      left_join(., tax_fb, by = "Genus_species") %>%
      mutate(fishbase = ifelse(!is.na(Family), Genus_species, NA)) %>%
      select(Genus_species, fishbase, Family, Order, SpecCode)
  }

  safe_species_list <- function(...) {
    out <- species_list(...)
    if (!length(out)) {
      return(NA_character_)
    }
    out
  }


  if (all(!is.na(x$fishbase))) {
    cat("All species names are up to date.")
  } else {


    # get synonyms for species that aren't recognized
    syns <- x %>%
      filter(is.na(fishbase)) %>%
      mutate(code = purrr:::map(Genus_species, rfishbase:::synonyms)) %>%
      unnest() %>%
      filter(Status != "misapplied name") %>%
      mutate(fishbase_name = map_chr(SpecCode, ~ safe_species_list(SpecCode = .))) %>%
      select(SpecCode, Genus_species, fishbase_name = synonym)

    out <- left_join(x, syns, by = "Genus_species") %>%
      mutate(fishbase = ifelse(is.na(fishbase), fishbase_name, fishbase)) %>%
      left_join(., tax_fb, by = c("fishbase" = "Genus_species")) %>%
      mutate(
        Genus_species = gsub(" ", "_", Genus_species),
        fishbase = gsub(" ", "_", fishbase)
      ) %>%
      select(current_name = Genus_species, valid_name = fishbase, Family = Family.y, Order = Order.y) %>%
      unique()
  }
  out
}


# code tweaked from geiger function
subset.phylo <- function(x, taxonomy, rank = "", ncores = 1, drop.tax = TRUE, ...) {
  phy <- x
  if (!rank %in% colnames(taxonomy)) {
    stop(paste(sQuote(rank), " does not appear as a column name in 'taxonomy'",
      sep = ""
    ))
  }
  xx <- match(phy$tip.label, rownames(taxonomy))
  new <- as.matrix(cbind(tip = phy$tip.label, rank = taxonomy[
    xx,
    rank
  ]))
  drop <- apply(new, 1, function(x) {
    if (any(is.na(x)) | any(x ==
      "")) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  if (any(drop) & drop.tax == FALSE) {
    warning(paste("Information for some tips is missing from 'taxonomy'; offending tips will be left unpruned:\n\t",
      paste(phy$tip.label[drop], collapse = "\n\t"),
      sep = ""
    ))
  } else {
    warning(paste("Information for some tips is missing from 'taxonomy'; offending tips WILL BE pruned:\n\t",
      paste(phy$tip.label[drop], collapse = "\n\t"),
      sep = ""
    ))
    drop.sp <-  phy$tip.label[drop]
    phy <- drop.tip(phy, drop.sp)
    new <- new[!rownames(new) %in% drop.sp, ]
  }
  tips <- phy$tip.label
  hphy <- geiger::hashes.phylo(phy, tips = tips, ncores = ncores)
  tax <- as.data.frame(new, stringsAsFactors = FALSE)
  stax <- split(tax$tip, tax$rank)
  rank_hashes <- sapply(stax, function(ss) geiger:::.hash.tip(ss, tips = tips))
  pruned <- hphy
  if(drop.tax){
    pruned$tip.label <- new[,2]
  } else {
    pruned$tip.label <- ifelse(drop == TRUE, tax$tip, tax$rank)
  }
  if (!all(zz <- rank_hashes %in% hphy$hash) & drop.tax == FALSE) {
    warning(paste(paste("non-monophyletic at level of ",
      rank,
      sep = ""
    ), ":\n\t", paste(sort(nonmon <- names(rank_hashes)[!zz]),
      collapse = "\n\t"
    ), sep = ""))
    vv <- which(pruned$tip.label %in% nonmon)
    pruned$tip.label[vv] <- tax$tip[vv]
  } else {
    warning(
      paste(paste("non-monophyletic at level of ",
        rank,
        sep = ""
      ), ":\n\t", paste(sort(nonmon <- names(rank_hashes)[!zz]),
        collapse = "\n\t"
      ), sep = ""),
      "\n taxa in families will be pruned from tree."
    )
    vv <- which(pruned$tip.label %in% nonmon)
    drop.sp <- c(drop.sp, tax$tip[vv])
    pruned <- drop.tip(phy, drop.sp)
    tax <- tax[!tax$tip %in% drop.sp, ]
    pruned$tip.label <- tax$rank
    warning(paste("\nA total of", length(drop.sp), "tips dropped from tree."))
  }
  rank_phy <- geiger:::unique.phylo(pruned)
  rank_phy$tip.label <- as.character(rank_phy$tip.label)
  return(rank_phy)
}


# code from phytools blog to extract Q matrix
Q.mat <- function(x){
  Q<-matrix(NA,length(x$states),length(x$states))
  Q[]<-c(0,x$rates)[x$index.matrix+1]
  diag(Q)<-0
  diag(Q)<--rowSums(Q)
  colnames(Q)<-rownames(Q)<-x$states
  Q
}


areColors <- function(x) {
  sapply(x, function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  })
}


# code modified from phytools function plotTree.wBars
barplot_phylo <- function(x, bar.col, tree, bar.width, bar.offset, scale, 
                          tip.labels, tip.cex, tip.col, tip.offset, legend, legend.cex) {
  x <- if(hasArg(scale)) x * scale else x
  col <- if(hasArg(bar.col)) {
    if (all(bar.col %in% colors()) |
        all(areColors(bar.col))) 
      bar.col 
    else {
      l <- levels(as.factor(bar.col))
      l <- l[!is.na(l) & l != ""]
      gg_color <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
      lcol <- setNames(sapply(gg_color(length(l)), color.id), l)
      col <- lcol[bar.col]
      col[is.na(col) | col == ""] <- "grey30"
      col
    }
  } else "grey80"
  tree <- tree
  obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  if(obj$type != "fan") 
    stop("This function only works for a fan phylogeny. Re-plot phylogeny as a fan and then try this function again.")
  
  w <- if(hasArg(bar.width)) bar.width else 0.2
  h <- max(nodeHeights(tree))
  sw <- if(hasArg(bar.offset)) bar.offset else 1
  if(length(w) == 1){
    for(i in 1:length(x)){
      theta<-atan(obj$yy[i]/obj$xx[i])
      s<-if(obj$xx[i]>0) 1 else -1
      dx<-s*h*cos(theta)+s*cos(theta)*sw
      dy<-s*h*sin(theta)+s*sin(theta)*sw
      x1<-dx+(w/2)*cos(pi/2-theta)-s*min(0,min(x))*cos(theta)
      y1<-dy-(w/2)*sin(pi/2-theta)-s*min(0,min(x))*sin(theta)
      x2<-dx-(w/2)*cos(pi/2-theta)-s*min(0,min(x))*cos(theta)
      y2<-dy+(w/2)*sin(pi/2-theta)-s*min(0,min(x))*sin(theta)
      x3<-s*x[i]*cos(theta)+x2
      y3<-s*x[i]*sin(theta)+y2
      x4<-s*x[i]*cos(theta)+x1
      y4<-s*x[i]*sin(theta)+y1
      polygon(c(x1,x2,x3,x4),c(y1,y2,y3,y4),col=col[i], border="transparent")
    }
  } else {
    for(i in 1:length(x)){
      theta<-atan(obj$yy[i]/obj$xx[i])
      s<-if(obj$xx[i]>0) 1 else -1
      dx<-s*h*cos(theta)+s*cos(theta)*sw
      dy<-s*h*sin(theta)+s*sin(theta)*sw
      x1<-dx+(w[i]/2)*cos(pi/2-theta)-s*min(0,min(x))*cos(theta)
      y1<-dy-(w[i]/2)*sin(pi/2-theta)-s*min(0,min(x))*sin(theta)
      x2<-dx-(w[i]/2)*cos(pi/2-theta)-s*min(0,min(x))*cos(theta)
      y2<-dy+(w[i]/2)*sin(pi/2-theta)-s*min(0,min(x))*sin(theta)
      x3<-s*x[i]*cos(theta)+x2
      y3<-s*x[i]*sin(theta)+y2
      x4<-s*x[i]*cos(theta)+x1
      y4<-s*x[i]*sin(theta)+y1
      polygon(c(x1,x2,x3,x4),c(y1,y2,y3,y4),col=col[i], border="transparent")
    }
    
  }
  if(hasArg(tip.labels)){
    if(tip.labels == TRUE) {
      tips <- tree$tip.label
      tip.cex <- if(hasArg(tip.cex)) tip.cex else 1
      tip.col <- if(hasArg(tip.col)) tip.col else "black"
      tip.offset <- if(hasArg(tip.offset)) tip.offset else 0
      for(i in 1:length(tips)){
        tmp <- rect2polar(obj$xx[i], obj$yy[i])
        angle <- atan(obj$yy[i]/obj$xx[i])*180/pi
        tmp <- polar2rect(tmp$r + tip.offset, tmp$angle)
        
        if(obj$xx[i] < 0) {
          text(tmp$x,tmp$y, paste(gsub("_"," ",tips[i])," ", sep=""),
               pos=2,srt=angle, offset=0, cex=tip.cex, col=tip.col)
        } else {
          text(tmp$x,tmp$y, paste(" ",gsub("_"," ",tips[i]),sep=""),
               pos=4,srt=angle, offset=0, cex=tip.cex, col=tip.col)
        }
        
      }
      
    }
    
    if(hasArg(legend) & legend == TRUE & all(!bar.col %in% colors()))  {
      legend.cex <- if(hasArg(legend.cex)) legend.cex else 1
      legend("topleft", names(lcol), col=lcol, bty = "n", pch=16, cex = legend.cex)
    }
  }
}























# this is the code copied directly from the geiger-v2 github repo here: https://github.com/mwpennell/geiger-v2/blob/master/R/congruify.R; was having issues getting code to work with geiger package loaded on computer

.build_classification=function(species){
  .split_lab=function(label){
    lab=gsub(".", "_", label, fixed=TRUE)
    lab=gsub(" ", "_", lab, fixed=TRUE)
    lab=unlist(strsplit(lab, "_", fixed=TRUE))
    lab=lab[lab!=""]
    lab
  }
  
  data.frame(genus=sapply(species, function(s) .split_lab(s)[1]), species=species, stringsAsFactors=FALSE)
}


.build_calibrations=function(dat, scion, scion_desc=NULL, tol=0){
  #	dat: branching times from stock; rows are 1:(Ntip(stock)+Nnode(stock))
  #			time                             hash
  #	1001 352.234677 3a4adb7cc0d4a51b9012dfb5615b3d71
  #	1002 243.269677 33757769ee61bde8dd5574ae35b47053
  
  #	scion: phylo tree with 'hash' object -- to be scaled from stock
  
  fetch_spanning=function(phy, nd, desc){
    # desc: a list from 1:(Ntip(phy)+Nnode(phy)) of tips descended from 'nd'
    if(nd<=Ntip(phy)) return(NULL)
    dd=geiger:::.get.desc.of.node(nd,phy)[1:2]
    tt=sapply(dd, function(x) return(desc[[x]][1]))
    return(phy$tip.label[sort(tt)])
  }
  
  if(is.null(scion_desc)) scion_desc=geiger:::.cache.descendants(scion)$tips
  
  N=Ntip(scion)
  stock_times=dat
  scion_hash=scion$hash
  df=data.frame(MRCA=scion_hash[(N+1):length(scion_hash)], MaxAge=NA, MinAge=NA, taxonA=NA, taxonB=NA, valid=FALSE, stringsAsFactors=FALSE)
  for(i in 1:nrow(df)){
    if(!is.na(hash.cur<-df$MRCA[i])){
      if(hash.cur%in%stock_times$hash){
        node.idx=i+N
        df[i,c("MaxAge","MinAge")]<-age.idx<-stock_times$time[match(hash.cur, stock_times$hash)]
        df[i,c("taxonA","taxonB")]<-taxa.idx<-fetch_spanning(scion, node.idx, scion_desc)
        if(age.idx>tol & all(!is.na(taxa.idx))) df[i,"valid"]=TRUE
      }
    }
  }
  df=df[df$valid,]
  return(df[,-which(names(df)=="valid")])
}

congruify.phylo=function(reference, target, taxonomy=NULL, tol=0, scale=c(NA, "PATHd8", "treePL"), ncores=NULL){
  ## adding requirement for ncbit
  ## require(ncbit)
  
  stock=reference
  scion=target
  #	stock: a time-calibrated phylogeny with tip-labels that can be treated as an exemplar for clades in 'scion'
  #		-- e.g., tip.label in 'stock' might be "Pinus" while in 'scion' we might have "Pinus_cembra"
  #		-- tips in 'stock' can be contained in 'scion' (FIXME: is this true?)
  #	megaphylogeny: a rooted phylogeny that is to be time-scaled based on 'stock'
  #	taxonomy: linkage table between tipsets for 'stock' and 'scion'; if empty, one is attempted to be built by 'scion' labels
  #		-- if NULL, 'stock' tips must correspond to tips in 'scion'... e.g., A, B, C in 'stock'; A_sp1, B_sp2, C_sp3 in 'scion'
  #		-- rownames of taxonomy must be tips in megaphylogeny
  
  ## functions
  method=match.arg(unname(sapply(scale, toString)), c("NA", "PATHd8", "treePL"))
  
  hashes.mapping <- function (phy, taxa, mapping){
    ## GENERAL FUNCTION: find hash tag for every edge in tree (using phy$tip.label or predefined set of 'taxa')
    # returns list of hash tags from node 1:(Nnode(phy)+Ntip(phy))
    # taxa: set of species used to create hash keys
    # mapping: named list -- names are tips in 'phy'; elements are tips represented by each tip in 'phy' (and should also be present in 'taxa')
    mapping=mapping[names(mapping)%in%phy$tip.label]
    if(is.null(taxa)) stop("Must supply 'tips'.")
    if(!all(names(mapping)%in%phy$tip.label)) stop("'mapping' must be named list with names in tip labels of 'phy'.")
    
    mapping=mapping[match(names(mapping), phy$tip.label)]
    descendants <- geiger:::.cache.descendants(phy)$tips
    hashes <- sapply(descendants, function(desc) geiger:::.hash.tip(unlist(mapping[desc]), taxa))
    empty=geiger:::.md5(integer(length(taxa)))
    hashes[hashes==empty]=NA
    phy$hash=hashes
    phy=geiger:::.uniquify_hashes(phy)
    
    return(phy)
  }
  
  times.mapping=function(phy, taxa, mapping){
    #	mapping: named list -- names are tips in 'phy'; elements are tips represented by each tip in 'phy' (and should also be present in 'taxa')
    #	taxa: species that are represented by tips in 'stock'
    
    # find hash tags for stock 'phy'
    stock=hashes.mapping(phy, taxa, mapping)
    tmp=geiger:::heights.phylo(stock)
    tmp$hash=stock$hash
    dat=data.frame(time=tmp[,"end"], hash=stock$hash, stringsAsFactors=FALSE)
    dat$hash[1:Ntip(stock)]=NA 	# destroy keys that are associated with tips
    
    return(list(stock=stock,dat=dat))
  }
  
  smooth_scion=function(stock, scion, scion_desc, taxa, spp, tol=0.01, method=c("PATHd8", NA, "treePL")){
    method=match.arg(toString(method), c("NA", "PATHd8", "treePL"))
    if(!is.ultrametric(stock, tol=tol)) warning("Supplied 'stock' is non-ultrametric.")
    stock_tmp=times.mapping(stock, taxa, spp)
    stock=stock_tmp$stock
    stock_dat=stock_tmp$dat
    calibration=geiger:::.build_calibrations(stock_dat, scion, scion_desc, tol=tol)
    if(!nrow(calibration)) {
      warning("No concordant branches reconciled between 'stock' and 'scion'; ensure that 'tax' involves rownames found as tip labels in 'scion'")
      return(NA)
    }
    if(method=="PATHd8") {
      phy=PATHd8.phylo(scion, calibration, base=".tmp_PATHd8", rm=FALSE)
      phy$hash=c(rep("", Ntip(phy)), phy$node.label)
      phy$hash[phy$hash==""]=NA
    } else if(method=="treePL") {
      phy=treePL.phylo(scion, calibration, base=".tmp_treePL", rm=FALSE)
      phy$hash=c(rep("", Ntip(phy)), phy$node.label)
      phy$hash[phy$hash==""]=NA
    } else if(method=="NA"){
      phy=NULL
    }
    
    stock$node.label=stock$hash[(Ntip(stock)+1):max(stock$edge)]
    stock$node.label[is.na(stock$node.label)]=""
    return(list(phy=phy, calibrations=calibration, reference=stock, target=scion))
  }
  
  ## end functions
  
  ## PROCESSING ##
  classification=taxonomy
  unfurl=FALSE
  
  if(inherits(stock,"phylo")) {
    stock=list(stock)
    unfurl=TRUE
  }
  if(is.null(classification)) {
    classification=as.data.frame(unique(as.matrix(.build_classification(scion$tip.label)),MARGIN=2))
  }
  
  tips=unique(unlist(lapply(stock, function(x) x$tip.label)))
  spp=lapply(tips, function(o) {
    x=rownames(classification)[which(classification==o, arr.ind=TRUE)[,1]]
    x=x[x%in%scion$tip.label]
  })
  names(spp)=tips
  taxa=unique(unlist(spp))
  
  scion=geiger:::hashes.phylo(scion, taxa, ncores)
  scion_desc=geiger:::.cache.descendants(scion)$tips
  if(is.null(scion$edge.length)) scion$edge.length=numeric(nrow(scion$edge)) ## JME 01302013
  
  f=lapply
  results=f(1:length(stock), function(i) {
    phy=stock[[i]]
    smooth_scion(phy, scion, scion_desc, taxa, spp, tol=tol, method=method)
  })
  
  if(unfurl) results=results[[1]]
  return(results)
}

## END CONGRUIFICATION FUNCTIONS ##


write.treePL=function(phy, calibrations, nsites=10000, min=0.0001, base="", opts=list(smooth=100, nthreads=8, optad=0, opt=1, cvstart=1000, cviter=3, cvend=0.1, thorough=TRUE)){
  #	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB' from .build_calibrations
  #	MRCA							MaxAge     MinAge                                  taxonA                                  taxonB
  #	c65bacdf65aa29635bec90f3f0447c6e 352.234677 352.234677                          Inga_chartacea             Encephalartos_umbeluziensis
  #	d4bc6557dccbd4e8e18b867979f34f8e 243.269677 243.269677                          Inga_chartacea                     Nuphar_sagittifolia
  #	5e30c16531694aff7d94da3b35806677 217.632627 217.632627                          Inga_chartacea                  Schisandra_glaucescens
  
  if(file.exists(inp<-paste(base,"infile",sep="."))) unlink(inp)
  if(file.exists(int<-paste(base,"intree",sep="."))) unlink(inp)
  
  poss=list(
    cv="numeric",
    collapse="boolean",
    checkconstraints="boolean",
    cvstart="numeric",
    cvstop="numeric",
    cvmultstep="numeric",
    verbose="boolean",
    lftemp="numeric",
    pltemp="numeric",
    plcool="numeric",
    lfstoptemp="numeric",
    plstoptemp="numeric",
    lfrtstep="numeric",
    plrtstep="numeric",
    thorough="boolean",
    lfiter="integer",
    pliter="integer",
    cviter="integer",
    ldfsimaniter="integer",
    plsimaniter="integer",
    cvsimaniter="integer",
    calcgrad="numeric",
    paramverbose="boolean",
    prime="boolean",
    opt="boolean",
    optad="boolean",
    optcvad="boolean",
    moredetail="boolean",
    moredetailad="boolean",
    moredetailcvad="boolean",
    randomcv="boolean",
    ftol="numeric",
    xtol="numeric",
    mapspace="boolean",
    nthreads="integer"
  )
  if(length(opts)==0) {
    print(poss)
    stop("No 'opts' specified")
  }
  
  # correct small branch lengths
  z=phy$edge.length[which(phy$edge.length>0)]
  if(any(z<min)){
    scl=min/min(z)
    phy$edge.length=phy$edge.length*scl
  }
  write.tree(phy, file=int)
  
  ##	check appropriateness of constraints ##
  #	check for 'calibrations' and 'phy' congruence
  #	if(!is.null(phy)){
  #		check=function(t, phy) all(t%in%phy$tip.label)
  #		a=check(calibrations$taxonA, phy)
  #		b=check(calibrations$taxonB, phy)
  
  #		if(!all(c(a,b))) stop("Some calibrations not encountered in tree.")
  #	}
  
  ##	build r8s file
  #	calibrations$fixage=ifelse(calibrations$MinAge==calibrations$MaxAge, TRUE, FALSE)
  constraints<-constraintnames<-character(nrow(calibrations))
  for(i in 1:nrow(calibrations)){
    cal=calibrations[i,]
    taxon=cal$MRCA
    desc=c(cal$taxonA, cal$taxonB)
    
    txt1=ifelse(!is.na(cal$MinAge), paste("min =", taxon, cal$MinAge, sep=" "), "")
    txt2=ifelse(!is.na(cal$MaxAge), paste("max =", taxon, cal$MaxAge, sep=" "), "")
    txt=paste(txt1,txt2,sep="\n")
    
    constraints[i]=txt
    constraintnames[i]=paste("mrca =", taxon, desc[1], desc[2], sep=" ")
  }
  infile=list(
    tree=paste("treefile = ", int, sep=""),
    ns=paste("numsites = ", nsites, sep=""),
    names=paste(unlist(constraintnames), collapse="\n"),
    mrca=paste(unlist(constraints), collapse="\n"),
    out=paste("outfile = ", paste(base, "dated", "tre", sep="."), sep=""),
    opt=paste(names(opts), opts, sep="=", collapse="\n")
  )
  
  inp=paste(base,"infile",sep=".")
  writeLines(paste(infile,collapse="\n\n"), con=inp)
  attr(inp, "method")="treePL"
  return(inp)
}

write.r8s=function(phy=NULL, calibrations, base="", blformat=c(lengths="persite", nsites=10000, ultrametric="no", round="yes"), divtime=c(method="NPRS", algorithm="POWELL"), describe=c(plot="chrono_description"), cv=c(cvStart=0, cvInc=0.5, cvNum=8), do.cv=FALSE){
  #	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB' from .build_calibrations
  #	MRCA							MaxAge     MinAge                                  taxonA                                  taxonB
  #	c65bacdf65aa29635bec90f3f0447c6e 352.234677 352.234677                          Inga_chartacea             Encephalartos_umbeluziensis
  #	d4bc6557dccbd4e8e18b867979f34f8e 243.269677 243.269677                          Inga_chartacea                     Nuphar_sagittifolia
  #	5e30c16531694aff7d94da3b35806677 217.632627 217.632627                          Inga_chartacea                  Schisandra_glaucescens
  
  if(file.exists(inp<-paste(base,"infile",sep="."))) unlink(inp)
  
  ##	check appropriateness of constraints ##
  #	check for 'calibrations' and 'phy' congruence
  #	if(!is.null(phy)){
  #		check=function(t, phy) all(t%in%phy$tip.label)
  #		a=check(calibrations$taxonA, phy)
  #		b=check(calibrations$taxonB, phy)
  
  #		if(!all(c(a,b))) stop("Some calibrations not encountered in tree.")
  #	}
  
  ##	build r8s file
  #	calibrations$fixage=ifelse(calibrations$MinAge==calibrations$MaxAge, TRUE, FALSE)
  constraints<-constraintnames<-character(nrow(calibrations))
  for(i in 1:nrow(calibrations)){
    cal=calibrations[i,]
    taxon=cal$MRCA
    desc=c(cal$taxonA, cal$taxonB)
    
    txt=paste(paste("\tfixage taxon=", taxon,sep=""), paste("age=", cal$MinAge, ";\n", sep=""), sep=" ")
    
    constraints[i]=txt
    constraintnames[i]=paste("\tMRCA", taxon, desc[1], paste(desc[2], ";\n", sep=""), sep=" ")
  }
  #	phy$node.label=NULL
  cv.code <- ""
  if(do.cv) {
    cv.code <- 	paste(paste(names(cv), cv, sep="="), ";\n", sep="", collapse="")
  }
  infile=paste(c(
    "#nexus\n",
    "begin trees;\n",
    paste("tree r8s = ", write.tree(phy), "\n", sep=""),
    "end;\n",
    "begin r8s;\n",
    paste("\tblformat ", paste(names(blformat), blformat, collapse=" ", sep="="), ";\n", sep=""),
    names=paste(unlist(constraintnames), collapse=""),
    mrca=paste(unlist(constraints), collapse=""),
    "\tcollapse;\n",
    paste("\tdivtime ", paste(names(divtime), divtime, collapse=" ", sep="="), cv.code, ";\n", sep=""),
    paste("\tdescribe ", paste(names(describe), describe, sep="="), ";\n", sep="", collapse=""),
    "end;"
  ),collapse="")
  
  inp=paste(base,"infile",sep=".")
  writeLines(paste(infile,collapse="\n\n"), con=inp)
  attr(inp, "method")="r8s"
  return(inp)
}


write.pathd8=function(phy, calibrations, base=""){
  #	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB'
  
  if(file.exists(inp<-paste(base,"infile",sep="."))) unlink(inp)
  
  ##	check appropriateness of constraints ##
  #	check for 'calibrations' and 'phy' congruence
  check=function(t, phy) all(t%in%phy$tip.label)
  a=check(calibrations$taxonA, phy)
  b=check(calibrations$taxonB, phy)
  
  if(!all(c(a,b))) stop("Some calibrations not encountered in tree.")
  
  ##	build PATHd8 file
  calibrations$fixage=ifelse(calibrations$MinAge==calibrations$MaxAge, TRUE, FALSE)
  constraints<-constraintnames<-character(nrow(calibrations))
  for(i in 1:nrow(calibrations)){
    cal=calibrations[i,]
    taxon=cal$MRCA
    desc=c(cal$taxonA, cal$taxonB)
    if(cal$fixage) {
      txt=paste("mrca:", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("fixage=", cal$MinAge, ";", sep=""), sep=" ")
    } else {
      txt1=paste("mrca:", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("minage=", cal$MinAge, ";", sep=""), sep=" ")
      txt2=paste("mrca:", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("maxage=", cal$MaxAge, ";", sep=""), sep=" ")
      txt=paste(txt1,txt2,sep="\n")
    }
    constraints[i]=txt
    constraintnames[i]=paste("name of mrca: ", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("name=", cal$MRCA, ";", sep=""), sep=" ")
  }
  phy$node.label=NULL
  infile=list(tree=write.tree(phy),
              mrca=paste(unlist(constraints), collapse="\n"),
              names=paste(unlist(constraintnames), collapse="\n")
  )
  
  inp=paste(base,"infile",sep=".")
  writeLines(paste(infile,collapse="\n\n"), con=inp)
  attr(inp, "method")="pathd8"
  return(inp)
}

PATHd8.phylo=function(phy, calibrations=NULL, base="", rm=TRUE){
  #	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB'
  #		-- if NULL, simple ultrametricization of 'phy' is performed
  
  phy$node.label=NULL
  if(!is.null(calibrations)){
    infile=write.pathd8(phy, calibrations, base)
  } else {
    infile=paste(base, "infile", sep=".")
    write.tree(phy, infile)
  }
  smooth.file=paste(base, "smoothed.tre", sep=".")
  parsed.outfile=paste(base, "pathd8.out", sep=".")
  outfile=paste(base, "pathd8.orig.out", sep=".")
  if(file.exists(outfile)) unlink(outfile)
  if(!system("which PATHd8", ignore.stdout=TRUE)==0) stop("Install 'PATHd8' before proceeding.")
  system(paste("PATHd8 -n", infile, "-pn >", outfile, sep=" "))
  system(paste("grep \"d8 tree\" ", outfile, ">", parsed.outfile, sep=" "))
  smoothed=read.tree(parsed.outfile)
  if(rm & base=="") {
    unlink(parsed.outfile)
    unlink(smooth.file)
    unlink(outfile)
    unlink(infile)
  }
  return(smoothed)
}

treePL.phylo=function(phy, calibrations=NULL, base="", rm=TRUE){
  phy$node.label=NULL
  if(!is.null(calibrations)){
    infile=write.treePL(phy=phy, calibrations=calibrations, base=base)
  } else {
    infile=paste(base, "infile", sep=".")
    write.tree(phy, infile)
  }
  smooth.file=paste(base, "dated.tre", sep=".")
  outfile=paste(base, "treePL.orig.out", sep=".")
  if(file.exists(outfile)) unlink(outfile)
  if(!system("which treePL", ignore.stdout=TRUE)==0) stop("Install 'treePL' before proceeding.")
  system(paste("treePL ", infile, " >", outfile, sep=" "))
  #system(paste("grep \"tree\" ", outfile, ">", parsed.outfile, sep=" "))
  smoothed=read.tree(smooth.file)
  if(rm & base=="") {
    unlink(smooth.file)
    unlink(outfile)
    unlink(infile)
  }
  return(smoothed)
}

r8s.phylo=function(phy, calibrations=NULL, base="r8srun", ez.run="none", rm=TRUE,  blformat=c(lengths="persite", nsites=10000, ultrametric="no", round="yes"), divtime=c(method="NPRS", algorithm="POWELL"), cv=c(cvStart=0, cvInc=0.5, cvNum=8), do.cv=FALSE){
  if(grepl("nprs",ez.run, ignore.case=TRUE)) {
    divtime <- c(method="NPRS", algorithm="POWELL")
    do.cv <- FALSE
  }
  if(grepl("pl",ez.run, ignore.case=TRUE)) {
    divtime <- c(method="PL", algorithm="qnewt")
    cv <- c(cvStart=0, cvInc=0.5, cvNum=8)
    do.cv <- TRUE
  }
  #	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB'
  #		-- if NULL, simple ultrametricization of 'phy' is performed
  
  phy$node.label=NULL
  if(!is.null(calibrations)){
    infile=write.r8s(phy, calibrations, base, blformat=blformat, divtime=divtime, cv=cv, do.cv=do.cv)
  } else {
    infile=paste(base, "infile", sep=".")
    write.tree(phy, infile)
  }
  smooth.file=paste(base, "smoothed.tre", sep=".")
  parsed.outfile=paste(base, "r8s.out", sep=".")
  outfile=paste(base, "r8s.orig.out", sep=".")
  if(file.exists(outfile)) unlink(outfile)
  if(!system("which r8s", ignore.stdout=TRUE)==0) stop("Install 'r8s' before proceeding.")
  system(paste("r8s -b -f", infile, " >", outfile, sep=" "))
  system(paste("grep \"tree r8s\" ", outfile, ">", parsed.outfile, sep=" "))
  smoothed=read.tree(parsed.outfile)
  if(rm & base=="") {
    unlink(parsed.outfile)
    unlink(smooth.file)
    unlink(outfile)
    unlink(infile)
  }
  return(smoothed)
}
