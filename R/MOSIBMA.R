#' MOSIBMA
#'
#' Use swarm intelligence-based algorithm to search for designs that satisfy the
#' corresponding minimum aberration criterion.
#'
#'
#'
#'
#' @param balance If set to \code{TRUE}, balanced designs are considered and
#' the argument \code{full_design} is not required.
#' If the user considers unbalanced designs, the argument should set to \code{TRUE}.
#' @param full_design A data.frame stands for the complete factorial design.
#' Required only when \code{balance} is set to \code{FALSE}.
#' @param factor_level A list where each element is a numeric vector specifying
#' levels of a factor.
#' @param unit An integer representing the number of experimental units in a design.
#' @param particle_number An integer indicating the initial particle number.
#' @param particle_increase An integer indicating the number of particles increased sequentially.
#' @param SIB_time An integer indicating the number of MOSIBMA processes that will run simultaneously.
#' Should larger than 1.
#' @param G_MA_list A list where each element is a list of matrices standing for the orthogonal
#' projection matrices corresponding to a minimum aberration criterion.
#' @param q_GB An integer indicating the number of columns of a design that
#' should be mixed with the corresponding columns of a design in the global best particle.
#' @param q_LB An integer indicating the number of columns of a design that
#' should be mixed with the corresponding columns of a design in the local best particle.
#' @param q_new An integer indicating the number of columns of a design that
#' should be mixed with the corresponding columns of a random design.
#' @param t An integer indicating the initial iteration number.
#' @param t_increase An integer indicating the increased iteration number sequentially.
#' @param structure_matrix A list where each element is an incidence matrix.
#' Required only when groups of units need to have the same value of some factors, e.g., split-plot
#' or strip-plot design, and so on.
#' @param structure_factor A list where each element is a numeric vector specifying
#' which factors should have the same value in terms of an incidence matrix.
#' @param initial_design A list of data.frames representing designs that are promising to be minimum
#' aberration designs.
#' @param iteration An integer indicating the total iteration number of MOSIBMA.
#'
#' @return
#' A list containing two elements:
#' \describe{
#' \item{GMA_design_candidate}{A list of candidate designs.}
#' \item{GMA_design.info}{A list consists \code{G_MA_index} and \code{G_MA_WL}.
#' \code{G_MA_index} represents the index of minimum aberration designs in \code{GMA_design_candidate}.
#' \code{G_MA_WL} represents the wordlength patterns of the minimum aberration designs.
#' The order of \code{G_MA_index} and \code{G_MA_WL} follows the order of \code{G_MA_list}.}
#' }
#'
#' @export

MOSIBMA <- function(balance = TRUE, full_design = NULL, factor_level, unit,
                    particle_number = 10, particle_increase = 10,
                    SIB_time = 3, G_MA_list,
                    q_GB = 1, q_LB = 1, q_new = 1, t = 10,
                    t_increase = 10, structure_matrix = NULL,
                    structure_factor = NULL, initial_design = NULL,
                    iteration = 10){
  a <- 0 ## the number of convergence occurred
  ################## step1: create particles
  all_particle <- lapply(1:SIB_time, function(i){
    p <- create_particle(balance = balance, full_design = full_design,
                         factor_level = factor_level,
                         unit = unit,
                         particle_number = particle_number,
                         structure_matrix = structure_matrix,
                         structure_factor = structure_factor)
    return(append(p, initial_design))
  })
  #names(all_particle) <- paste0("SIB_time", 1:SIB_time)
  ################## global information
  factor_number <- length(factor_level)
  factor_level_number <- lengths(factor_level)
  model.matrix_text <- paste0("~.^", factor_number)
  weight <- sort(100^rep(1:factor_number), decreasing = TRUE)
  G_MA_list_number <- length(G_MA_list)
  ################## step2: MOSIBMA
  #### t=1
  # get GB in index
  GB <- lapply(all_particle, function(i){
    compare(particle = i, G_MA_list = G_MA_list, factor_number = factor_number,
            factor_level_number = factor_level_number,
            model.matrix_text = model.matrix_text, weight = weight,
            G_MA_list_number = G_MA_list_number)
  })
  # mix with GB
  mixwGB <- lapply(1:SIB_time, function(i){
    mix_op <- lapply(1:length(all_particle[[i]]), function(j){
      mix_operation(X_design = all_particle[[i]][j],
                    Y_design = all_particle[[i]][GB[[i]]$G_MA_index],
                    factor_number = factor_number, G_MA_list = G_MA_list, q = q_GB,
                    factor_level_number = factor_level_number,
                    model.matrix_text = model.matrix_text, weight = weight,
                    G_MA_list_number = G_MA_list_number)
    })
  })
  # move X, X is designs not index
  X <- lapply(1:SIB_time, function(i){
    X_op <- lapply(1:length(all_particle[[i]]), function(j){
      move_X(balance = balance, full_design = full_design, factor_level = factor_level,
             unit = unit, factor_number = factor_number,
             structure_matrix = structure_matrix, structure_factor = structure_factor,
             q_new = q_new, factor_level_number = factor_level_number,
             model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
             weight = weight, G_MA_list_number = G_MA_list_number,
             X_design = all_particle[[i]][j], mixwGB = mixwGB[[i]][[j]], mixwLB = NULL)
    })
  })
  # move LB, LB is designs
  LB <- lapply(1:SIB_time, function(i){
    LB_op <- lapply(1:length(all_particle[[i]]), function(j){
      move_GBLB(LB = all_particle[[i]][j], X_design = X[[i]][[j]],
                factor_number = factor_number, factor_level_number = factor_level_number,
                model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                weight = weight, G_MA_list_number = G_MA_list_number)
    })
  })
  # move GB, GB is index
  GB <- lapply(1:SIB_time, function(i){
    move_GBLB(LB = LB[[i]], X_design = NULL,
              factor_number = factor_number, factor_level_number = factor_level_number,
              model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
              weight = weight, G_MA_list_number = G_MA_list_number)
  })
  #### t>1
  for(i in 2:t){
    # mix with GB
    mixwGB <- lapply(1:SIB_time, function(i){
      if(length(GB[[i]]$G_MA_index)==1){
        Y_design <- LB[[i]][[GB[[i]]$G_MA_index]]
      }else{
        Y_design <- unlist(lapply(1:G_MA_list_number, function(k){
          if(length(LB[[i]][[GB[[i]]$G_MA_index[k]]]) == 1){
            return(LB[[i]][[GB[[i]]$G_MA_index[k]]][1])
          }else{
            return(LB[[i]][[GB[[i]]$G_MA_index[k]]][k])
          }
          }), recursive = FALSE)
      }
      mix_op <- lapply(1:length(X[[i]]), function(j){
        mix_operation(X_design = X[[i]][[j]],
                      Y_design = Y_design,
                      factor_number = factor_number, G_MA_list = G_MA_list, q = q_GB,
                      factor_level_number = factor_level_number,
                      model.matrix_text = model.matrix_text, weight = weight,
                      G_MA_list_number = G_MA_list_number)
      })
    })
    # mix with LB
    mixwLB <- lapply(1:SIB_time, function(i){
      mix_op <- lapply(1:length(X[[i]]), function(j){
        mix_operation(X_design = X[[i]][[j]],
                      Y_design = LB[[i]][[j]],
                      factor_number = factor_number, G_MA_list = G_MA_list, q = q_GB,
                      factor_level_number = factor_level_number,
                      model.matrix_text = model.matrix_text, weight = weight,
                      G_MA_list_number = G_MA_list_number)
      })
    })
    # move X, X is designs not index
    X <- lapply(1:SIB_time, function(i){
      X_op <- lapply(1:length(X[[i]]), function(j){
        move_X(balance = balance, full_design = full_design, factor_level = factor_level,
               unit = unit, factor_number = factor_number,
               structure_matrix = structure_matrix, structure_factor = structure_factor,
               q_new = q_new, factor_level_number = factor_level_number,
               model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
               weight = weight, G_MA_list_number = G_MA_list_number,
               X_design = X[[i]][[j]], mixwGB = mixwGB[[i]][[j]],
               mixwLB = mixwLB[[i]][[j]])
      })
    })
    # move LB, LB is designs
    LB <- lapply(1:SIB_time, function(i){
      LB_op <- lapply(1:length(LB[[i]]), function(j){
        move_GBLB(LB = LB[[i]][[j]], X_design = X[[i]][[j]],
                  factor_number = factor_number, factor_level_number = factor_level_number,
                  model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                  weight = weight, G_MA_list_number = G_MA_list_number)
      })
    })
    # move GB, GB is index
    GB <- lapply(1:SIB_time, function(i){
      move_GBLB(LB = LB[[i]], X_design = NULL,
                factor_number = factor_number, factor_level_number = factor_level_number,
                model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                weight = weight, G_MA_list_number = G_MA_list_number)
    })
  }
  ################## step3: MOSIBMA--increase t or particle number
  ii <- 1
  repeat{
    # make sure each SIB's WL are the same
    while((same_WL(GB = GB, SIB_time = SIB_time) == FALSE) & (ii < iteration)){ # each SIB's WL are not the same--increase t
      current_GB <- GB
      for(i in 1:t_increase){
        # mix with GB
        mixwGB <- lapply(1:SIB_time, function(i){
          if(length(GB[[i]]$G_MA_index)==1){
            Y_design <- LB[[i]][[GB[[i]]$G_MA_index]]
          }else{
            Y_design <- unlist(lapply(1:G_MA_list_number, function(k){
              if(length(LB[[i]][[GB[[i]]$G_MA_index[k]]]) == 1){
                return(LB[[i]][[GB[[i]]$G_MA_index[k]]][1])
              }else{
                return(LB[[i]][[GB[[i]]$G_MA_index[k]]][k])
              }
            }), recursive = FALSE)
          }
          mix_op <- lapply(1:length(X[[i]]), function(j){
            mix_operation(X_design = X[[i]][[j]],
                          Y_design = Y_design,
                          factor_number = factor_number, G_MA_list = G_MA_list, q = q_GB,
                          factor_level_number = factor_level_number,
                          model.matrix_text = model.matrix_text, weight = weight,
                          G_MA_list_number = G_MA_list_number)
          })
        })
        # mix with LB
        mixwLB <- lapply(1:SIB_time, function(i){
          mix_op <- lapply(1:length(X[[i]]), function(j){
            mix_operation(X_design = X[[i]][[j]],
                          Y_design = LB[[i]][[j]],
                          factor_number = factor_number, G_MA_list = G_MA_list, q = q_GB,
                          factor_level_number = factor_level_number,
                          model.matrix_text = model.matrix_text, weight = weight,
                          G_MA_list_number = G_MA_list_number)
          })
        })
        # move X, X is designs not index
        X <- lapply(1:SIB_time, function(i){
          X_op <- lapply(1:length(X[[i]]), function(j){
            move_X(balance = balance, full_design = full_design, factor_level = factor_level,
                   unit = unit, factor_number = factor_number,
                   structure_matrix = structure_matrix, structure_factor = structure_factor,
                   q_new = q_new, factor_level_number = factor_level_number,
                   model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                   weight = weight, G_MA_list_number = G_MA_list_number,
                   X_design = X[[i]][[j]], mixwGB = mixwGB[[i]][[j]],
                   mixwLB = mixwLB[[i]][[j]])
          })
        })
        # move LB, LB is designs
        LB <- lapply(1:SIB_time, function(i){
          LB_op <- lapply(1:length(LB[[i]]), function(j){
            move_GBLB(LB = LB[[i]][[j]], X_design = X[[i]][[j]],
                      factor_number = factor_number, factor_level_number = factor_level_number,
                      model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                      weight = weight, G_MA_list_number = G_MA_list_number)
          })
        })
        # move GB, GB is index
        GB <- lapply(1:SIB_time, function(i){
          move_GBLB(LB = LB[[i]], X_design = NULL,
                    factor_number = factor_number, factor_level_number = factor_level_number,
                    model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                    weight = weight, G_MA_list_number = G_MA_list_number)
        })
      }
      ii <- ii+1
      if((no_improvement_WL(current_GB = current_GB, new_GB = GB) == TRUE) & (ii < iteration)){ # no improvement for increasing t--increase particle number
        # increase particle number
        added_particle <- lapply(1:SIB_time, function(i){
          x <- create_particle(balance = balance, full_design = full_design,
                               factor_level = factor_level, unit = unit,
                               particle_number = particle_increase,
                               structure_matrix = structure_matrix,
                               structure_factor = structure_factor)
          return(lapply(x, function(j){list(j)}))
        })
        # update X
        X <- lapply(1:SIB_time, function(i){
          return(append(X[[i]], added_particle[[i]]))
        })
        # update LB
        LB <- lapply(1:SIB_time, function(i){
          return(append(LB[[i]], added_particle[[i]]))
        })
        # update GB in index
        GB <- lapply(1:SIB_time, function(i){
          move_GBLB(LB = LB[[i]], X_design = NULL,
                    factor_number = factor_number, factor_level_number = factor_level_number,
                    model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                    weight = weight, G_MA_list_number = G_MA_list_number)
        })
        for(i in 2:t_increase){
          # mix with GB
          mixwGB <- lapply(1:SIB_time, function(i){
            if(length(GB[[i]]$G_MA_index)==1){
              Y_design <- LB[[i]][[GB[[i]]$G_MA_index]]
            }else{
              Y_design <- unlist(lapply(1:G_MA_list_number, function(k){
                if(length(LB[[i]][[GB[[i]]$G_MA_index[k]]]) == 1){
                  return(LB[[i]][[GB[[i]]$G_MA_index[k]]][1])
                }else{
                  return(LB[[i]][[GB[[i]]$G_MA_index[k]]][k])
                }
              }), recursive = FALSE)
            }
            mix_op <- lapply(1:length(X[[i]]), function(j){
              mix_operation(X_design = X[[i]][[j]],
                            Y_design = Y_design,
                            factor_number = factor_number, G_MA_list = G_MA_list, q = q_GB,
                            factor_level_number = factor_level_number,
                            model.matrix_text = model.matrix_text, weight = weight,
                            G_MA_list_number = G_MA_list_number)
            })
          })
          # mix with LB
          mixwLB <- lapply(1:SIB_time, function(i){
            mix_op <- lapply(1:length(X[[i]]), function(j){
              mix_operation(X_design = X[[i]][[j]],
                            Y_design = LB[[i]][[j]],
                            factor_number = factor_number, G_MA_list = G_MA_list, q = q_GB,
                            factor_level_number = factor_level_number,
                            model.matrix_text = model.matrix_text, weight = weight,
                            G_MA_list_number = G_MA_list_number)
            })
          })
          # move X, X is designs not index
          X <- lapply(1:SIB_time, function(i){
            X_op <- lapply(1:length(X[[i]]), function(j){
              move_X(balance = balance, full_design = full_design, factor_level = factor_level,
                     unit = unit, factor_number = factor_number,
                     structure_matrix = structure_matrix, structure_factor = structure_factor,
                     q_new = q_new, factor_level_number = factor_level_number,
                     model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                     weight = weight, G_MA_list_number = G_MA_list_number,
                     X_design = X[[i]][[j]], mixwGB = mixwGB[[i]][[j]],
                     mixwLB = mixwLB[[i]][[j]])
            })
          })
          # move LB, LB is designs
          LB <- lapply(1:SIB_time, function(i){
            LB_op <- lapply(1:length(LB[[i]]), function(j){
              move_GBLB(LB = LB[[i]][[j]], X_design = X[[i]][[j]],
                        factor_number = factor_number, factor_level_number = factor_level_number,
                        model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                        weight = weight, G_MA_list_number = G_MA_list_number)
            })
          })
          # move GB, GB is index
          GB <- lapply(1:SIB_time, function(i){
            move_GBLB(LB = LB[[i]], X_design = NULL,
                      factor_number = factor_number, factor_level_number = factor_level_number,
                      model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                      weight = weight, G_MA_list_number = G_MA_list_number)
          })
        }
        ii <- ii+1
      }
    }
    if(!(ii<iteration)){
      print("not converged")
      break
    } # first condition for breaking the repeat loop
    a <- a+1
    if(a==1){
      pre_finalized_GB <- GB
    }else{
      post_finalized_GB <- GB
    }
    if(a>1){
      if(no_improvement_WL(current_GB = pre_finalized_GB, new_GB = post_finalized_GB)==TRUE){ # second condition for breaking the repeat loop
        break
      }else{
        pre_finalized_GB <- post_finalized_GB
      }
    }
    # increase particle number
    added_particle <- lapply(1:SIB_time, function(i){
      x <- create_particle(balance = balance, full_design = full_design,
                           factor_level = factor_level, unit = unit,
                           particle_number = particle_increase,
                           structure_matrix = structure_matrix,
                           structure_factor = structure_factor)
      return(lapply(x, function(j){list(j)}))
    })
    # update X
    X <- lapply(1:SIB_time, function(i){
      return(append(X[[i]], added_particle[[i]]))
    })
    # update LB
    LB <- lapply(1:SIB_time, function(i){
      return(append(LB[[i]], added_particle[[i]]))
    })
    # update GB in index
    GB <- lapply(1:SIB_time, function(i){
      move_GBLB(LB = LB[[i]], X_design = NULL,
                factor_number = factor_number, factor_level_number = factor_level_number,
                model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                weight = weight, G_MA_list_number = G_MA_list_number)
    })
    for(i in 2:t_increase){
      # mix with GB
      mixwGB <- lapply(1:SIB_time, function(i){
        if(length(GB[[i]]$G_MA_index)==1){
          Y_design <- LB[[i]][[GB[[i]]$G_MA_index]]
        }else{
          Y_design <- unlist(lapply(1:G_MA_list_number, function(k){
            if(length(LB[[i]][[GB[[i]]$G_MA_index[k]]]) == 1){
              return(LB[[i]][[GB[[i]]$G_MA_index[k]]][1])
            }else{
              return(LB[[i]][[GB[[i]]$G_MA_index[k]]][k])
            }
          }), recursive = FALSE)
        }
        mix_op <- lapply(1:length(X[[i]]), function(j){
          mix_operation(X_design = X[[i]][[j]],
                        Y_design = Y_design,
                        factor_number = factor_number, G_MA_list = G_MA_list, q = q_GB,
                        factor_level_number = factor_level_number,
                        model.matrix_text = model.matrix_text, weight = weight,
                        G_MA_list_number = G_MA_list_number)
        })
      })
      # mix with LB
      mixwLB <- lapply(1:SIB_time, function(i){
        mix_op <- lapply(1:length(X[[i]]), function(j){
          mix_operation(X_design = X[[i]][[j]],
                        Y_design = LB[[i]][[j]],
                        factor_number = factor_number, G_MA_list = G_MA_list, q = q_GB,
                        factor_level_number = factor_level_number,
                        model.matrix_text = model.matrix_text, weight = weight,
                        G_MA_list_number = G_MA_list_number)
        })
      })
      # move X, X is designs not index
      X <- lapply(1:SIB_time, function(i){
        X_op <- lapply(1:length(X[[i]]), function(j){
          move_X(balance = balance, full_design = full_design, factor_level = factor_level,
                 unit = unit, factor_number = factor_number,
                 structure_matrix = structure_matrix, structure_factor = structure_factor,
                 q_new = q_new, factor_level_number = factor_level_number,
                 model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                 weight = weight, G_MA_list_number = G_MA_list_number,
                 X_design = X[[i]][[j]], mixwGB = mixwGB[[i]][[j]],
                 mixwLB = mixwLB[[i]][[j]])
        })
      })
      # move LB, LB is designs
      LB <- lapply(1:SIB_time, function(i){
        LB_op <- lapply(1:length(LB[[i]]), function(j){
          move_GBLB(LB = LB[[i]][[j]], X_design = X[[i]][[j]],
                    factor_number = factor_number, factor_level_number = factor_level_number,
                    model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                    weight = weight, G_MA_list_number = G_MA_list_number)
        })
      })
      # move GB, GB is index
      GB <- lapply(1:SIB_time, function(i){
        move_GBLB(LB = LB[[i]], X_design = NULL,
                  factor_number = factor_number, factor_level_number = factor_level_number,
                  model.matrix_text = model.matrix_text, G_MA_list = G_MA_list,
                  weight = weight, G_MA_list_number = G_MA_list_number)
      })
    }
    ii <- ii+1
  }
  ################## step4: return the G-MA design(s)
  GMA_design_candidate <- unlist(unlist(LB, recursive = FALSE), recursive = FALSE)
  GMA_design.info <- compare(particle = GMA_design_candidate,
                             G_MA_list = G_MA_list, factor_number = factor_number,
                             factor_level_number = factor_level_number,
                             model.matrix_text = model.matrix_text,
                             weight = weight, G_MA_list_number = G_MA_list_number)
  ans <- list("GMA_design_candidate" = GMA_design_candidate, "GMA_design.info" = GMA_design.info)
  return(ans)
}
