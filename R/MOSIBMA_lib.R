## Create a list of designs
create_particle <- function(balance, full_design, factor_level, unit, particle_number,
                            structure_matrix, structure_factor){
  if(balance == TRUE){
    # Return a list. Each element is an integer representing each level's occurrence
    # number of a factor.
    replicate <- lapply(seq_along(factor_level), function(i){
      return(unit / length(factor_level[[i]]))
    })
    # Return a data.frame representing a design before randomized.
    standard_design <- sapply(seq_along(factor_level), function(i){
      return(rep(factor_level[[i]], each=replicate[[i]]))
    })
    standard_design <- data.frame(standard_design)
    # Return a list of data.frame. Each represents a design, which is also a particle.
    # Each particle is the result of randomization of standard_design by factor.
    particle <- lapply(1:particle_number, function(i){
      standard_design[] <- lapply(standard_design, sample)
      return(standard_design)
    })
  }else{
    particle <- lapply(1:particle_number, function(i){
      design <- full_design[sample(1:nrow(full_design), unit),]
      return(design)
    })
  }

  # When the design has a structure.
  if(is.null(structure_matrix) == FALSE){
    idx_same_treatment <- lapply(structure_matrix, function(i){
      return(apply(i, 2, function(j) which(j == 1)))
    })
    # Return a list.
    particle <- lapply(particle, function(i){
      unique_treatment <- lapply(1:length(idx_same_treatment), function(j){
        return(sapply(structure_factor[[j]], function(k){
          return(sample(rep(factor_level[[k]],
                            each = ncol(idx_same_treatment[[j]])/length(factor_level[[k]])),
                        ncol(idx_same_treatment[[j]])))
        }))
      })
      for(j in 1:length(idx_same_treatment)){
        for(k in 1:ncol(idx_same_treatment[[j]])){
          for(l in 1:length(structure_factor[[j]])){
            i[idx_same_treatment[[j]][,k], structure_factor[[j]][l]] <- unique_treatment[[j]][k,l]
          }
        }
      }
      return(i)
    })
  }
  return(particle)
}

# create wordlength pattern for a design under a P_w
wordlength <- function(m, P_w, factor_number){
  result <- rep(0, factor_number)
  colnames_m <- colnames(m)
  ind_set <- stringr::str_count(colnames_m, pattern = ":")[-1]
  for(i in 1:factor_number){
    ind <- which(ind_set == (i-1)) + 1
    result[i] <- psych::tr(t(m[,ind]) %*% P_w %*% m[,ind])
  }
  result <- round(result, digits = 3)
  return(result)
}
# compare designs under a SIB_time
#' @importFrom psych tr
compare <- function(particle, G_MA_list, factor_number, factor_level_number,
                    model.matrix_text, weight, G_MA_list_number){
  # create the model matrix
  M <- lapply(particle, function(i){
    namei <- names(i)
    i[,namei] <- lapply(i[,namei], factor)
    for(j in 1:factor_number){
      #contrasts(i[,j]) <- contr.poly(factor_level_number[j])
      contrasts(i[,j]) <- contr.poly(as.numeric(levels(i[,j])))
    }
    m <- model.matrix(eval(parse(text = model.matrix_text)), i)
    m <- apply(m, 2, function(x){x/sqrt(t(x)%*%x)})
    return(m)
  })
  # find the wordlength pattern and the best designs under each G_MA
  WL <- lapply(G_MA_list, function(i){
    a_WL <- lapply(M, function(j){  # under a G_MA, generate designs' WL
      a_WL_list <- lapply(1:length(i), function(k){
        return(wordlength(m = j, P_w = i[[k]], factor_number = factor_number))
      })
      return(Reduce("+", a_WL_list))
    })
    weight_a_WL <- lapply(a_WL, function(j){
      sum(j*weight)
    })
    weight_a_WL <- unlist(weight_a_WL)
    min_index <- which(weight_a_WL == min(weight_a_WL)) # choose all that reach the min
    result <- list("min_index" = min_index, "min_a_WL" = a_WL[min_index])
    return(result)
  })
  # find the G_MA design(s)
  if(G_MA_list_number == 1){
    G_MA_info <- list("G_MA_index" = WL[[1]]$min_index[1], "G_MA_WL" = WL[[1]]$min_a_WL[[1]])
  }else{
    unlist_WL <- unlist(WL, recursive = FALSE)
    min_index_list <- unlist_WL[which(names(unlist_WL) == "min_index")]
    intersect_index <- Reduce(intersect, min_index_list)
    if(length(intersect_index) == 0){ # find intersection, if no, admissible designs
      r <- 1:G_MA_list_number
      G_MA_index <- lapply(1:G_MA_list_number, function(i){
        final_index <- min_index_list[[i]]
        for(j in r[-i]){
          if(length(final_index) == 1){break}
          a_WL2 <- lapply(M[final_index], function(k){
            a_WL_list2 <- lapply(1:length(G_MA_list[[j]]), function(l){
              return(wordlength(m = k, P_w = G_MA_list[[j]][[l]], factor_number = factor_number))
            })
            return(Reduce("+", a_WL_list2))
          })
          weight_a_WL2 <- lapply(a_WL2, function(k){
            sum(k*weight)
          })
          weight_a_WL2 <- unlist(weight_a_WL2)
          final_index <- final_index[which(weight_a_WL2 == min(weight_a_WL2))]
        }
        return(final_index[1])
      })
      G_MA_index <- unlist(G_MA_index)
      G_MA_WL <- lapply(1:G_MA_list_number, function(i){WL[[i]]$min_a_WL[[ which(WL[[i]]$min_index==G_MA_index[i]) ]]})
      G_MA_info <- list("G_MA_index" = G_MA_index, "G_MA_WL" = G_MA_WL) # G_MA_index refers to the index of G_MA designs in the particle under the G_MA
    }else{ # if yes, one design
      G_MA_WL <- lapply(1:G_MA_list_number, function(i){WL[[i]]$min_a_WL[[ which(WL[[i]]$min_index==intersect_index[1]) ]]})
      G_MA_info <- list("G_MA_index" = intersect_index[1], "G_MA_WL" = G_MA_WL)
    }
  }
  return(G_MA_info)
}
# replace q columns of X_design with Y_design, where X_design and Y_design are good under the P_w.
# particle are designs of a SIB_time: e.g. all_particle[[1]]
# should report designs, not index
mix_operation <- function(X_design, Y_design, factor_number, G_MA_list, q, factor_level_number,
                          model.matrix_text, weight, G_MA_list_number){
  for(i in 1:q){
    if(length(X_design)==1 & length(Y_design)==1){
      RX <- lapply(1:factor_number, function(j){ # RX: reduced X
        X_design[[1]][,j] <- Y_design[[1]][,j]
        return(X_design)
      })
      RX <- unlist(RX, recursive = FALSE)
      X_info <- compare(particle = RX, G_MA_list = G_MA_list, factor_number = factor_number,
                               factor_level_number = factor_level_number,
                               model.matrix_text = model.matrix_text, weight = weight,
                               G_MA_list_number = G_MA_list_number)
      X_design <- RX[X_info$G_MA_index]
    }else if(length(X_design)==1 & length(Y_design)!=1){
      X_design <- lapply(1:G_MA_list_number, function(j){
        RX <- lapply(1:factor_number, function(k){
          X_design[[1]][,k] <- Y_design[[j]][,k]
          return(X_design)
        })
        RX <- unlist(RX, recursive = FALSE)
        X_info <- compare(particle = RX, G_MA_list = G_MA_list[j], factor_number = factor_number,
                                 factor_level_number = factor_level_number,
                                 model.matrix_text = model.matrix_text, weight = weight,
                                 G_MA_list_number = 1)
        return(RX[X_info$G_MA_index])
      })
      X_design <- unlist(X_design, recursive = FALSE)
    }else if(length(X_design)!=1 & length(Y_design)==1){
      X_design <- lapply(1:G_MA_list_number, function(j){
        RX <- lapply(1:factor_number, function(k){
          X_design[[j]][,k] <- Y_design[[1]][,k]
          return(X_design)
        })
        RX <- unlist(RX, recursive = FALSE)
        X_info <- compare(particle = RX, G_MA_list = G_MA_list[j], factor_number = factor_number,
                                 factor_level_number = factor_level_number,
                                 model.matrix_text = model.matrix_text, weight = weight,
                                 G_MA_list_number = 1)
        return(RX[X_info$G_MA_index])
      })
      X_design <- unlist(X_design, recursive = FALSE)
    }else{
      X_design <- lapply(1:G_MA_list_number, function(j){
        RX <- lapply(1:factor_number, function(k){
          X_design[[j]][,k] <- Y_design[[j]][,k]
          return(X_design)
        })
        RX <- unlist(RX, recursive = FALSE)
        X_info <- compare(particle = RX, G_MA_list = G_MA_list[j], factor_number = factor_number,
                                 factor_level_number = factor_level_number,
                                 model.matrix_text = model.matrix_text, weight = weight,
                                 G_MA_list_number = 1)
        return(RX[X_info$G_MA_index])
      })
      X_design <- unlist(X_design, recursive = FALSE)
    }
  }
  return(X_design)
}
# mix with a new design
mix_with_new <- function(balance, full_design, factor_level, unit, factor_number,
                         structure_matrix, structure_factor, X_design, q_new,
                         factor_level_number, model.matrix_text,
                         weight, G_MA_list_number, G_MA_list){
  new_design <- create_particle(balance = balance, full_design = full_design,
                                factor_level = factor_level,
                                unit = unit, particle_number = 1,
                                structure_matrix = structure_matrix,
                                structure_factor = structure_factor)
  for(i in 1:q_new){
    if(length(X_design)==1){
      RX <- lapply(1:factor_number, function(j){ # RX: reduced X
        X_design[[1]][,j] <- new_design[[1]][,j]
        return(X_design)
      })
      RX <- unlist(RX, recursive = FALSE)
      X_info <- compare(particle = RX, G_MA_list = G_MA_list, factor_number = factor_number,
                        factor_level_number = factor_level_number,
                        model.matrix_text = model.matrix_text, weight = weight,
                        G_MA_list_number = G_MA_list_number)
      X_design <- RX[X_info$G_MA_index]
    }else{
      X_design <- lapply(1:G_MA_list_number, function(j){
        RX <- lapply(1:factor_number, function(k){
          X_design[[j]][,k] <- new_design[[1]][,k]
          return(X_design)
        })
        RX <- unlist(RX, recursive = FALSE)
        X_info <- compare(particle = RX, G_MA_list = G_MA_list[j], factor_number = factor_number,
                          factor_level_number = factor_level_number,
                          model.matrix_text = model.matrix_text, weight = weight,
                          G_MA_list_number = 1)
        return(RX[X_info$G_MA_index])
      })
      X_design <- unlist(X_design, recursive = FALSE)
    }
  }
  return(X_design)
}


# move X
move_X <- function(balance, full_design, factor_level, unit, factor_number,
                   structure_matrix, structure_factor, q_new,
                   factor_level_number, model.matrix_text, G_MA_list,
                   weight, G_MA_list_number, X_design, mixwGB, mixwLB){
  if(is.null(mixwLB)){ # no mixwLB
    if(length(X_design)==1 & length(mixwGB)==1){ # length(X_design) and length(mixwGB) == 1
      candidate <- Reduce(append, list(X_design, mixwGB))
      candidate.info <- compare(particle = candidate, G_MA_list = G_MA_list,
                                factor_number = factor_number,
                                factor_level_number = factor_level_number,
                                model.matrix_text = model.matrix_text, weight = weight,
                                G_MA_list_number = G_MA_list_number)
      candidate_temp <- candidate[candidate.info$G_MA_index]
      which_G_MA_index <- which(candidate.info$G_MA_index == 1)
      for(i in which_G_MA_index){
        candidate_temp[i] <- mix_with_new(balance = balance, full_design = full_design,
                                          factor_level = factor_level,
                                          unit = unit, factor_number = factor_number,
                                          structure_matrix = structure_matrix,
                                          structure_factor = structure_factor,
                                          X_design = X_design, q_new = q_new,
                                          factor_level_number = factor_level_number,
                                          model.matrix_text = model.matrix_text,
                                          weight = weight, G_MA_list_number = 1,
                                          G_MA_list = G_MA_list[i])
      }
      return(candidate_temp)
    }else{ # length(X_design) or length(mixwGB) != 1
      if(length(X_design)==1){
        X_design <- X_design[rep(1,G_MA_list_number)]
      }
      if(length(mixwGB)==1){
        mixwGB <- mixwGB[rep(1,G_MA_list_number)]
      }
      candidate.info <- lapply(1:G_MA_list_number, function(i){
        candidate <- Reduce(append, list(X_design[i], mixwGB[i]))
        candidate.info2 <- compare(particle = candidate, G_MA_list = G_MA_list[i],
                                   factor_number = factor_number,
                                   factor_level_number = factor_level_number,
                                   model.matrix_text = model.matrix_text, weight = weight,
                                   G_MA_list_number = 1)
        if(candidate.info2$G_MA_index == 1){
          return(mix_with_new(balance = balance, full_design = full_design, factor_level = factor_level,
                              unit = unit, factor_number = factor_number,
                              structure_matrix = structure_matrix,
                              structure_factor = structure_factor, X_design = X_design[i], q_new = q_new,
                              factor_level_number = factor_level_number, model.matrix_text = model.matrix_text,
                              weight = weight, G_MA_list_number = 1,
                              G_MA_list = G_MA_list[i]))
        }else{return(candidate[candidate.info2$G_MA_index])}
      })
      candidate.info <- unlist(candidate.info, recursive = FALSE)
      return(candidate.info)
    }
  }else{ #  mixwLB exists
    if(length(X_design)==1 & length(mixwGB)==1 & length(mixwLB)==1){ # length(X_design) and length(mixwGB) and length(mixwLB) == 1
      candidate <- Reduce(append, list(X_design, mixwGB, mixwLB))
      candidate.info <- compare(particle = candidate, G_MA_list = G_MA_list,
                                factor_number = factor_number,
                                factor_level_number = factor_level_number,
                                model.matrix_text = model.matrix_text, weight = weight,
                                G_MA_list_number = G_MA_list_number)
      candidate_temp <- candidate[candidate.info$G_MA_index]
      which_G_MA_index <- which(candidate.info$G_MA_index == 1)
      for(i in which_G_MA_index){
        candidate_temp[i] <- mix_with_new(balance = balance, full_design = full_design,
                                          factor_level = factor_level,
                                          unit = unit, factor_number = factor_number,
                                          structure_matrix = structure_matrix,
                                          structure_factor = structure_factor,
                                          X_design = X_design, q_new = q_new,
                                          factor_level_number = factor_level_number,
                                          model.matrix_text = model.matrix_text,
                                          weight = weight, G_MA_list_number = 1,
                                          G_MA_list = G_MA_list[i])
      }
      return(candidate_temp)
    }else{ # length(X_design) or length(mixwGB) or length(mixwLB) != 1
      if(length(X_design)==1){
        X_design <- X_design[rep(1,G_MA_list_number)]
      }
      if(length(mixwGB)==1){
        mixwGB <- mixwGB[rep(1,G_MA_list_number)]
      }
      if(length(mixwLB)==1){
        mixwLB <- mixwLB[rep(1,G_MA_list_number)]
      }
      candidate.info <- lapply(1:G_MA_list_number, function(i){
        candidate <- Reduce(append, list(X_design[i], mixwGB[i], mixwLB[i]))
        candidate.info2 <- compare(particle = candidate, G_MA_list = G_MA_list[i],
                                   factor_number = factor_number,
                                   factor_level_number = factor_level_number,
                                   model.matrix_text = model.matrix_text, weight = weight,
                                   G_MA_list_number = 1)
        if(candidate.info2$G_MA_index == 1){
          return(mix_with_new(balance = balance, full_design = full_design, factor_level = factor_level,
                              unit = unit, factor_number = factor_number,
                              structure_matrix = structure_matrix,
                              structure_factor = structure_factor, X_design = X_design[i], q_new = q_new,
                              factor_level_number = factor_level_number, model.matrix_text = model.matrix_text,
                              weight = weight, G_MA_list_number = 1,
                              G_MA_list = G_MA_list[i]))
        }else{return(candidate[candidate.info2$G_MA_index])}
      })
      candidate.info <- unlist(candidate.info, recursive = FALSE)
      return(candidate.info)
    }
  }
}
# move GB or LB
move_GBLB <- function(LB, X_design, factor_number,factor_level_number,
                      model.matrix_text, G_MA_list, weight, G_MA_list_number){
  if(is.null(X_design)){ # move GB
    if(all(lengths(LB) == 1)){
      LB <- unlist(LB, recursive = FALSE)
      candidate.info <- compare(particle = LB, G_MA_list = G_MA_list,
                                factor_number = factor_number,
                                factor_level_number = factor_level_number,
                                model.matrix_text = model.matrix_text, weight = weight,
                                G_MA_list_number = G_MA_list_number)
    }else{
      LB <- lapply(1:length(LB), function(i){
        if(lengths(LB[i]) == 1){
          LB[[i]] <- LB[[i]][rep(1,G_MA_list_number)]
        }else{
          LB[[i]] <- LB[[i]]
        }
      })
      candidate.info <- lapply(1:G_MA_list_number, function(i){
        candidate <- purrr::map(LB, i)
        candidate.info2 <- compare(particle = candidate, G_MA_list = G_MA_list[i],
                                   factor_number = factor_number,
                                   factor_level_number = factor_level_number,
                                   model.matrix_text = model.matrix_text, weight = weight,
                                   G_MA_list_number = 1)
        return(candidate.info2)
      })
      x <- unlist(candidate.info, recursive = FALSE)
      G_MA_index <- unlist(x[which(names(x)=="G_MA_index")])
      if(length(Reduce(intersect, G_MA_index)) != 0){
        G_MA_index <- Reduce(intersect, G_MA_index)
      }
      G_MA_WL <- unname(x[which(names(x)=="G_MA_WL")])
      candidate.info <- list("G_MA_index" = G_MA_index, "G_MA_WL" = G_MA_WL)
    }
    return(candidate.info)
  }else{ # move LB
    if(length(LB)==1 & length(X_design)==1){
      candidate <- append(LB, X_design)
      candidate.info2 <- compare(particle = candidate, G_MA_list = G_MA_list,
                                 factor_number = factor_number,
                                 factor_level_number = factor_level_number,
                                 model.matrix_text = model.matrix_text, weight = weight,
                                 G_MA_list_number = G_MA_list_number)
      candidate.info <- candidate[candidate.info2$G_MA_index]
    }else{
      single_X <- 0
      single_LB <- 0
      if(length(X_design)==1){
        X_design <- X_design[rep(1,G_MA_list_number)]
        single_X <- 1
      }
      if(length(LB)==1){
        LB <- LB[rep(1,G_MA_list_number)]
        single_LB <- 1
      }
      candidate.info <- lapply(1:G_MA_list_number, function(i){
        candidate <- append(LB[i], X_design[i])
        candidate.info2 <- compare(particle = candidate, G_MA_list = G_MA_list[i],
                                   factor_number = factor_number,
                                   factor_level_number = factor_level_number,
                                   model.matrix_text = model.matrix_text, weight = weight,
                                   G_MA_list_number = 1)
        result <- list("candidate_design" = candidate[candidate.info2$G_MA_index],
                       "index" = candidate.info2$G_MA_index)
        return(result)
      })
      x <- unlist(candidate.info, recursive = FALSE)
      if(all(x[which(names(x)=="index")]==1)&&single_LB==1){
        candidate.info <- LB[1]
      }else if(all(x[which(names(x)=="index")]==2)&single_X==1){
        candidate.info <- X_design[1]
      }else{
        candidate.info <- unname(unlist(x[which(names(x)=="candidate_design")], recursive = FALSE))
      }
    }
    return(candidate.info)
  }
}
# check whether the G_MA_index under several SIB_time are the same or not
same_WL <- function(GB, SIB_time){
  x <- purrr::map(GB, 2)
  y <- paste0("x[[",1:SIB_time,"]]", collapse = ",")
  y <- paste0("DescTools::AllIdentical(",y,")")
  return(eval(parse(text = y)))
}

no_improvement_WL <- function(current_GB, new_GB){
  ans <- lapply(1:length(current_GB), function(i){
    DescTools::AllIdentical(current_GB[[i]]$G_MA_WL, new_GB[[i]]$G_MA_WL)
  })
  return(all(unlist(ans)))
}


