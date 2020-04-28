# test if a value is continuous
is_continuous <- function(x) {
  if (length(x) == 1) {TRUE} else {identical(as.numeric(x), seq(min(x), max(x), 1))}
}

# reconstructions for multiple variables
recs_all_process <- function(tree, character_table, method = c("ML", "MP"), 
                             variable_names = NULL, save_plot = FALSE, 
                             write_results = TRUE, tip_offset = 0.1, 
                             width = 5, height = 0.7, no.margin = TRUE, 
                             label.offset = 0.2, use.edge.length = FALSE, 
                             save_binary = FALSE,
                             output_directory = "Reconstruction_results") {
  
  len <- ifelse(class(tree)[1] == "multiPhylo", length(tree), 1)
  lec <- ifelse(class(character_table)[1] == "list", length(character_table), 1)
  
  if (is.null(variable_names)) {
    variable_names <- names(character_table)
    if (is.null(variable_names)) {
      stop("'variable_names' must be defined if 'character_table' is not named.")
    }
  } else {
    if (lec != length(variable_names)) {
      stop("length of 'variable_names' and 'character_table' must be the same.")
    }
  }
  
  if (len > 1) {
    if (write_results == TRUE | save_binary == TRUE) {
      dir.create(output_directory)
    }
    output_directory <- paste0(output_directory, "/Results_for_tree_", 1:len)
  }
  
  message("\nAncestral reconstructions have starter, please wait:")
  results <- lapply(1:len, function(z) {
    if (write_results == TRUE) {
      ## files to be created
      rec_files <- paste0(output_directory[z], "/", variable_names, "_rec_table.csv")
      pdf_files <- paste0(output_directory[z], "/", variable_names, "_niche_evol.pdf")
    }
    res <- lapply(1:length(variable_names), function(x) {
      ## tree and data
      if (len == 1) {tr <- tree} else {tr <- tree[[z]]}
      if (lec == 1) {ct <- character_table} else {ct <- character_table[[x]]}
      tree_data <- geiger::treedata(tr, ct)
      
      ## reconstructions
      if (method[1] %in% c("ML", "MP")) {
        if (method[1] == "ML") {
          recons <- bin_ml_rec(tree_data)
        } else {
          recons <- bin_par_rec(tree_data)
        }
      } else {
        stop("Argument 'method' is not valid.")
      }
      
      ## smoothing
      rec_smooth <- smooth_rec(recons)
      
      ## saving results
      if (write_results == TRUE) {
        dir.create(output_directory[z])
        write.csv(rec_smooth, rec_files[x], row.names = T)
        
        ## plot
        if (save_plot == TRUE) {
          pdf(pdf_files[x], width = 7, height = 3)
          par(mfrow = c(1:2))
          par(cex = 0.85)
          plot(tree_data$phy, no.margin = no.margin, label.offset = label.offset, 
               use.edge.length = use.edge.length)
          niche_labels(tree_data$phy, rec_smooth, tip_offset = tip_offset, 
                       width = width, height = height)
          niche_legend("topright", cex = 0.7)
          plot(tree_data$phy, no.margin = no.margin, label.offset = label.offset, 
               use.edge.length = use.edge.length)
          niche_labels(tree_data$phy, rec_smooth, tip_offset = tip_offset, 
                       label_type = "tip", width = width, height = height)
          nichevol_labels(tree_data$phy, rec_smooth, width = width, height = height)
          nichevol_legend("topright", cex = 0.7)
          dev.off()
        } 
      }
      
      return(rec_smooth)
    })
    
    if (save_binary == TRUE) {
      if (write_results == FALSE) {dir.create(output_directory[z])}
      save(res, file = paste0(output_directory[z], "/reconstruction_list.Rdata"))
    }
    message("\tProcess ", z, " of ", len, " finished")
    return(res)
  })
  
  if (len == 1) {return(results[[1]])} else {return(results)}
}

# variability among multiple reconstructions obtained with trees of distinct branch length
rec_variability <- function(rec_list) {
  if (missing(rec_list)) {
    stop("Argument 'rec_list' needs to be defined.")
  }
  
  # replacing ? for 0.5 and changing to numeric
  for (h in 1:length(rec_list)) {
    rec_list[[h]][rec_list[[h]] == "?"] <- "0.5"
    if ("LogLik" %in% rownames(rec_list[[h]])) {
      rec_list[[h]] <- rec_list[[h]][-c((nrow(rec_list[[h]]) - 2):nrow(rec_list[[h]])), ]
    }
    class(rec_list[[h]]) <- "numeric" 
  }
  
  # variability
  tab <- do.call(rbind, lapply(rec_list, function(y) {c(y)}))
  vars <- apply(tab, 2, mean)
  
  # filling matrix
  rec_varp <- rec_list[[1]]
  rec_varp[] <- vars
  
  return(rec_varp)
}

# labels to show of variability in reconstructions
niche_labels_var <- function(tree, whole_var_table, label_type = "tip_node",
                             tip_offset = 0.015, present_col = "#e41a1c", 
                             unknown_col = "#969696", absent_col = "#377eb8", 
                             width = 1, height = 1) {
  if (missing(tree)) {stop("Argument tree needs to be defined.")}
  if (missing(whole_var_table)) {stop("Argument whole_var_table needs to be defined.")}
  
  # preparing colors
  my_pal <- colorRampPalette(c(absent_col, unknown_col, present_col))
  cols <- my_pal(100)
  
  # reorganizing table
  tlab <- tree$tip.label
  nrt <- nrow(whole_var_table)
  rns <- c(tlab, rownames(whole_var_table)[(length(tlab) + 1):nrt])
  whole_var_table <- rbind(whole_var_table[tlab, ],
                           whole_var_table[(length(tlab) + 1):nrt, ])
  rownames(whole_var_table) <- rns
  
  # getting info from plot
  tp_info <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  if (tp_info$type != "phylogram") {
    stop("niche_labels can be used only for plots of type phylogram.")
  }
  if (!tp_info$direction %in% c("rightwards", "leftwards")) {
    stop("niche_labels can be used only for rightwards or leftwards phylograms.")
  }
  if (tp_info$direction == "leftwards") {tip_offset <- -tip_offset}
  
  xx <- tp_info$xx
  yy <- tp_info$yy
  edges <- tp_info$edge
  tpos <- 1:tp_info$Ntip
  npos <- (tp_info$Ntip + 1):(tp_info$Ntip + tp_info$Nnode)
  otips <- xx[tpos]
  rtip <- range(otips)
  if ((rtip[2] - rtip[1]) <= 0.00001) {
    xx[tpos] <- rep(max(otips), length(otips))
  }
  
  # organizing data
  tpol <- ncol(whole_var_table)
  wt <- ((max(yy) / tp_info$Ntip) / 10) * (height * 6)
  wpol <- wt / tpol
  
  h_vertices <- seq(0, wt, wpol)
  
  # plotting bars
  if (label_type %in% c("tip", "node", "tip_node")) {
    if (label_type %in% c("tip", "tip_node")) {
      barss <- sapply(tpos, function(j) {
        ys <- yy[j] - (wt / 2)
        hver <- ys + h_vertices
        wdt <- 0.01 * width
        xs <- xx[j] + tip_offset; xs1 <- xs - (wdt / 2); xs2 <- xs + (wdt / 2)
        wver <- rep(c(xs1, xs2), each = 2)
        
        polys <- sapply(1:(length(h_vertices) - 1), function(x) {
          pcolor <- cols[ifelse(round(whole_var_table[j, x] * 100) == 0, 1, 
                                round(whole_var_table[j, x] * 100))]
          
          yss <- c(hver[x], hver[x + 1], hver[x + 1], hver[x])
          
          polygon(x = wver, y = yss, col = pcolor, border = NA)
        })
      })
    }
    
    if (label_type %in% c("node", "tip_node")) {
      barss <- sapply(npos, function(j) {
        ys <- yy[j] - (wt / 2)
        hver <- ys + h_vertices
        wdt <- 0.01 * width
        xs <- xx[j]; xs1 <- xs - (wdt / 2); xs2 <- xs + (wdt / 2)
        wver <- rep(c(xs1, xs2), each = 2)
        
        polys <- sapply(1:(length(h_vertices) - 1), function(x) {
          pcolor <- cols[ifelse(round(whole_var_table[j, x] * 100) == 0, 1, 
                                round(whole_var_table[j, x] * 100))]
          
          yss <- c(hver[x], hver[x + 1], hver[x + 1], hver[x])
          
          polygon(x = wver, y = yss, col = pcolor, border = NA)
        })
      })
    }
    
  } else {
    stop("Argument 'label_type' is not correct, see help(niche_labels).")
  }
}

# find limits of niche or places where niche evolution occurs 
all_limits <- function(value_vector, split_list, type = c("niche", "change"),
                       changes_list = NULL) {
  if (missing(value_vector)) {stop("Argument 'value_vector' needs to be defined.")}
  if (missing(split_list)) {stop("Argument 'split_list' needs to be defined.")}
  if (type[1] == "change" & is.null(changes_list)) {
    stop("Argument 'changes_list' needs to be defined for detecting changes.")
  }
  
  l <- lapply(1:length(split_list), function(x) {
    if (type[1] == "change") {value_vector[] <- changes_list[[x]]}
    j <- split_list[[x]]
    if (length(j) > 0) {
      if (is_continuous(j)) {
        c(j[j == min(j)], j[j == max(j)])
      } else {
        ls <- vector()
        ad <- 1
        while (ad < (length(value_vector) - 1)) {
          ocon <- sum(c(ad, ad + 1) %in% j) >= 1
          if (value_vector[ad] != value_vector[ad + 1] & ocon) {
            if (sum(j %in% ad) >= 1) {
              ls[ad] <- j[j == ad]
            } else {
              ls[ad] <- j[j == (ad + 1)]
            }
          }
          ad <- ad + 1
        }
        gx <- j[j %in% na.omit(ls)]
        gt <- na.omit(ls)[na.omit(ls) %in% j]
        names(gt) <- sapply(gt, function(x) {names(gx[gx == x])})
        
        if (max(gt) != max(j)) {gt <- c(gt, j[j == max(j)])}
        if (1 %in% j) {gt <- c(j[1], gt)}
        return(gt)
      }
    } else {
      return(integer())
    }
  })
  if (type[1] == "change") {
    names(l) <- names(changes_list)
  } else {
    names(l) <- names(split_list)
  }
  return(l)
}

# find values of previous limits
critical_limits <- function(all_limits, return = c("values", "bins")) {
  if (return[1] == "values") {
    lapply(all_limits, function(x) {
      if (length(x) > 0) {
        v <- rep(1:2, ceiling(length(x) / 2))
        vals <- sapply(1:length(x), function (y) {
          strsplit(names(x)[y], " to ")[[1]][v[y]]
        })
        mat <- matrix(as.numeric(vals), ncol = 2, byrow = TRUE)
        colnames(mat) <- c("From", "To"); return(mat)
      } else {
        mat <- matrix(ncol = 2); colnames(mat) <- c("From", "To") 
        return(mat[!is.na(mat[, 1]), ])
      }
    })
  } else {
    lapply(all_limits, function(x) {
      mat <- matrix(x, ncol = 2, byrow = TRUE)
      colnames(mat) <- c("From", "To"); return(mat)
    })
  }
}

# find values of niche evolution using previous functions
find_evol_lims <- function(whole_rec_table, from, to, return = c("values", "bins"),
                           present = "1", absent = "0", unknown = "?") {
  if (missing(whole_rec_table)) {
    stop("Argument 'whole_rec_table' needs to be defined.")
  }
  if (missing(from)) {stop("Argument 'from' needs to be defined.")}
  if (missing(to)) {stop("Argument 'to' needs to be defined.")}
  
  # finding
  to_find <- whole_rec_table[c(from, to), ]
  tn <- whole_rec_table[1, ]
  
  ## changes
  expansion <- sapply(1:ncol(to_find), function(z) {
    from <- to_find[1, z]; to <- to_find[2, z]
    if (from == absent & to == present) {TRUE} else {FALSE}
  })
  
  retraction <- sapply(1:ncol(to_find), function(z) {
    from <- to_find[1, z]; to <- to_find[2, z]
    if (from == present & to == absent) {TRUE} else {FALSE}
  })
  
  stable <- !expansion & !retraction
  
  uns <- list(expansion = expansion, retraction = retraction, stable = stable)
  
  ## limits of changes
  lims <- lapply(uns, function(x) {tn[] <- x; which(tn == TRUE)})
  alims <- all_limits(value_vector = tn, split_list = lims, type = "change", 
                      changes_list = uns)
  
  ## critical values and bin numbers
  avals <- critical_limits(alims, return = "values")
  
  abins <- critical_limits(alims, return = "bins")
  
  # results
  if (sum(stable) < length(tn)) {
    if (return[1] == "values") {return(avals)} else {return(abins)}
  } else {
    warning("No changes detected in comparisons.")
    if (return[1] == "values") {return(avals)} else {return(abins)}
  }
}

# find values of niche limits using previous functions
find_niche_lims <- function(whole_rec_table, tip_or_node, return = c("values", "bins"),
                            present = "1", absent = "0", unknown = "?") {
  if (missing(whole_rec_table)) {
    stop("Argument 'whole_rec_table' needs to be defined.")
  }
  if (missing(tip_or_node)) {stop("Argument 'tip_or_node' needs to be defined.")}
  if (length(tip_or_node) > 1) {stop("Argument 'tip_or_node' must be of length 1.")}
  
  # finding
  tn <- whole_rec_table[tip_or_node, ]
  uns <- unique(tn)
  
  ns <- ifelse(uns == present, "present", ifelse(uns == absent, "absent", "unknown"))
  
  ## limits
  lims <- lapply(uns, function(x) {which(tn == x)})
  names(lims) <- ns
  
  alims <- all_limits(value_vector = tn, split_list = lims, type = "niche")
  
  ## critical values and bin numbers
  avals <- critical_limits(alims, return = "values")
  
  abins <- critical_limits(alims, return = "bins")
  
  # results
  if (return[1] == "values") {return(avals)} else {return(abins)}
}

# creates a raster layer to show areas of niche and niche evolution
nichevol_layer <- function(variable, niche_list, 
                           return = c("niche", "evolution", "nichevol"),
                           evol_list = NULL, id_unknown = TRUE) {
  if (missing(variable)) {stop("Argument 'variable' needs to be defined.")}
  if (class(variable)[1] != "RasterLayer") {
    stop("Argument 'variable' needs to be an object of class 'RasterLayer'.")
  }
  if (missing(niche_list)) {stop("Argument 'niche_list' needs to be defined.")}
  if (!return[1] %in% c("nichevol", "niche", "evolution")) {
    stop("Arguments 'return' is not valid.")
  }
  if (return[1] %in% c("nichevol", "evolution") & is.null(evol_list)) {
    stop("Arguments 'evol_list' needs to be defined to identify changes in niche.")
  }
  
  # preparing objects
  vals <- variable[]
  unk <- ifelse(id_unknown == TRUE, 10, 0)
  
  # reclassidication and details
  if (return[1] %in% c("niche", "nichevol")) {
    ns <- names(niche_list)
    nvs <- ifelse(ns == "present", 100, ifelse(ns == "absent", 0, unk))
    rec_l <- lapply(1:length(niche_list), function(x) {
      cbind(niche_list[[x]], becomes = nvs[x])
    })
    recm <- do.call(rbind, rec_l)
    recm <- recm[order(recm[, 1]), ]
    for (i in 2:nrow(recm)) {recm[i, 1] <- recm[i - 1, 2]}
    nnas <- range(c(recm[, 1:2]))
    
    ni <- variable
    ni <- raster::reclassify(ni, recm)
    ni[vals <= nnas[1] | vals > nnas[2]] <- NA
  }
  
  if (return[1] %in% c("evolution", "nichevol")) {
    pres <- niche_list$present
    ns <- names(evol_list)
    nvs <- ifelse(ns == "expansion", 1, ifelse(ns == "retraction", 2, 0))
    rec_l <- lapply(1:length(evol_list), function(x) {
      etab <- evol_list[[x]]
      if (nrow(etab) > 0) {
        if(nrow(etab) > 1) {
          vs <- rep(nvs[x], nrow(etab))
          if (nvs[x] == 1) {vs <- seq(1, 10, 2)[1:length(vs)]}
          if (nvs[x] == 2) {vs <- seq(2, 10, 2)[1:length(vs)]}
          cbind(etab, becomes = vs)
        } else {
          nvst <- ifelse(etab[1, 2] <= pres[1, 1], 0, ifelse(nvs[x] == 0, 0, 2))
          vs <- nvs[x] + nvst
        }
        cbind(etab, becomes = vs)
      } else {
        cbind(etab, becomes = numeric())
      }
    })
    recm <- do.call(rbind, rec_l)
    recm <- recm[order(recm[, 1]), ]
    for (i in 2:nrow(recm)) {recm[i, 1] <- recm[i - 1, 2]}
    nnas <- range(c(recm[, 1:2]))
    
    ev <- variable
    ev <- raster::reclassify(ev, recm)
    ev[vals <= nnas[1] | vals > nnas[2]] <- NA
  }
  
  # results
  if (return[1] == "nichevol") {return(raster::ratify(ni + ev))}
  if (return[1] == "niche") {return(raster::ratify(ni))}
  if (return[1] == "evolution") {return(raster::ratify(ev))}
}