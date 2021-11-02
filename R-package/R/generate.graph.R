gen_rr_adjmat <- function(n_node, degree, beta, alpha, type = c("ferro", "glass")) {
    g <- igraph::sample_k_regular(n_node, degree)
    adj <- as.matrix(igraph::as_adjacency_matrix(g, type = "both"))
    
    ind_nonzero <- which(adj != 0, arr.ind = TRUE)
    ind_nonzero <- ind_nonzero[ind_nonzero[, 1] > ind_nonzero[, 2],]
    num_nonzero <- nrow(ind_nonzero)
    
    if (type == "ferro") {
      value <- c(alpha, rep(beta, num_nonzero - 1))
    } else if (type == "glass") {
      len_pos <- round(num_nonzero / 2)
      len_neg <- num_nonzero - len_pos
      value_pos <- c(alpha, rep(beta, len_pos - 1))
      value_neg <- c(-alpha, rep(-beta, len_neg - 1))
      value <- c(value_pos, value_neg)
    }
    
    value <- sample(value, size = num_nonzero)
    for (i in 1:num_nonzero) {
      adj[ind_nonzero[i, 1], ind_nonzero[i, 2]] <- value[i]
      adj[ind_nonzero[i, 2], ind_nonzero[i, 1]] <- value[i]
    }
    adj
  }

gen_4nn_cyc <-
  function(n_node,
           degree,
           beta,
           alpha,
           type = c("ferro", "glass", "glass_weak")) {
    adj <- matrix(0, n_node, n_node)
    index_lr <-
      (row(adj) - col(adj) == 1 & col(adj) %% sqrt(n_node) != 0)
    index_ud <- (row(adj) - col(adj) == sqrt(n_node))
    cycle_lr <-
      (col(adj) %% sqrt(n_node) == 0 &
         col(adj) - row(adj) == sqrt(n_node) - 1)
    cycle_ud <-
      (col(adj) <= sqrt(n_node) &
         row(adj) - col(adj) == sqrt(n_node) * (sqrt(n_node) - 1))
    index <- index_lr | index_ud | cycle_lr | cycle_ud
    
    adj[index] <- 1
    adj <- adj + t(adj)
    
    ind_nonzero <- which(adj != 0, arr.ind = TRUE)
    ind_nonzero <- ind_nonzero[ind_nonzero[, 1] > ind_nonzero[, 2],]
    num_nonzero <- nrow(ind_nonzero)
    
    if (type == "ferro") {
      value <- c(alpha, rep(beta, num_nonzero - 1))
    } else if (type == "glass") {
      len_pos <- round(num_nonzero / 2)
      len_neg <- num_nonzero - len_pos
      value_pos <- c(alpha, rep(beta, len_pos - 1))
      value_neg <- c(-alpha, rep(-beta, len_neg - 1))
      value <- c(value_pos, value_neg)
    } else if (type == "glass_weak") {
      value <- c(-alpha, rep(beta, num_nonzero - 1))
    }
    
    value <- sample(value, size = num_nonzero, replace = FALSE)
    for (i in 1:num_nonzero) {
      adj[ind_nonzero[i, 1], ind_nonzero[i, 2]] <- value[i]
      adj[ind_nonzero[i, 2], ind_nonzero[i, 1]] <- value[i]
    }
    adj
  }

sim_theta <- function(p, type = 1, seed, beta, degree, alpha) {
  # chain structure: size = p - 1
  if (type == 1) {
    theta <- matrix(0, p, p)
    # for(i in 1: p) {
    #   theta[i, i] <- sample(c(-0.5, 0, 0.5), 1)
    # }
    for (i in 1:(p - 1)) {
      theta[i, i + 1] <- sample(c(0.5,-0.5), 1)
      theta[i + 1, i] <- theta[i, i + 1]
    }
  }
  # random graph
  if (type == 2) {
    # set.seed(seed)
    sigma_inv <-
      fastclime::fastclime.generator(
        n = 20,
        d = p,
        graph = "random",
        verbose = FALSE
      )
    sparse <- as.matrix(sigma_inv$theta)
    theta <- matrix(0, p, p)
    # for(i in 1:p) {
    #   theta[i, i] <- sample(c(-0.5, 0, 0.5), 1)
    # }
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        if (sparse[i, j] != 0) {
          theta[i, j] <- sample(c(0.5,-0.5), 1)
          theta[j, i] <- theta[i, j]
        }
      }
    }
  }
  # 4 nearest neighbor: size = 2 * p
  if (type == 3) {
    theta <- matrix(0, p, p)
    # for(i in 1:p) {
    #   theta[i, i] <- sample(c(-0.5, 0, 0.5), 1)
    # }
    sqr <- sqrt(p)
    for (i in 0:(sqr - 1)) {
      for (j in 1:(sqr - 1)) {
        theta[i * sqr + j, i * sqr + j + 1] <- sample(c(-0.5, 0.5), 1)
        theta[i * sqr + j + 1, i * sqr + j] <-
          theta[i * sqr + j, i * sqr + j + 1]
      }
    }
    for (i in 0:(sqr - 2)) {
      for (j in 1:sqr) {
        theta[i * sqr + j, i * sqr + j + sqr] <- sample(c(-0.5, 0.5), 1)
        theta[i * sqr + j + sqr, i * sqr + j] <-
          theta[i * sqr + j, i * sqr + j + sqr]
      }
    }
  }
  # 8 nearest neighbor: size = 4 * p
  if (type == 4) {
    theta <- matrix(0, p, p)
    # for(i in 1:p) {
    #   theta[i, i] <- sample(c(-0.5, 0, 0.5), 1)
    # }
    sqr <- sqrt(p)
    for (i in 0:(sqr - 1)) {
      for (j in 1:(sqr - 1)) {
        theta[i * sqr + j, i * sqr + j + 1] <- sample(c(-0.5, 0.5), 1)
        theta[i * sqr + j + 1, i * sqr + j] <-
          theta[i * sqr + j, i * sqr + j + 1]
      }
    }
    for (i in 0:(sqr - 2)) {
      for (j in 1:sqr) {
        theta[i * sqr + j, i * sqr + j + sqr] <- sample(c(-0.5, 0.5), 1)
        theta[i * sqr + j + sqr, i * sqr + j] <-
          theta[i * sqr + j, i * sqr + j + sqr]
      }
    }
    for (i in 0:(sqr - 2)) {
      for (j in 1:(sqr - 1)) {
        theta[i * sqr + j, i * sqr + j + sqr + 1] <- sample(c(-0.5, 0.5), 1)
        theta[i * sqr + j + sqr + 1, i * sqr + j] <-
          theta[i * sqr + j, i * sqr + j + sqr + 1]
      }
    }
    for (i in 0:(sqr - 2)) {
      for (j in 2:sqr) {
        theta[i * sqr + j, i * sqr + j + sqr - 1] <- sample(c(-0.5, 0.5), 1)
        theta[i * sqr + j + sqr - 1, i * sqr + j] <-
          theta[i * sqr + j, i * sqr + j + sqr - 1]
      }
    }
  }
  # star shape: size = p - 1
  if (type == 5) {
    theta <- matrix(0, p, p)
    # for(i in 1:p) {
    #   theta[i, i] <- sample(c(-0.5, 0, 0.5), 1)
    # }
    for (i in 2:p) {
      theta[1, i] <- sample(c(-0.5, 0.5), 1)
      theta[i, 1] <- theta[1, i]
    }
  }
  # local dense graph: 6*6
  if (type == 6) {
    theta <- matrix(0, p, p)
    # for(i in 1:p) {
    #   theta[i, i] <- sample(c(-0.5, 0, 0.5), 1)
    # }
    for (i in 1:5) {
      for (j in (i + 1):6) {
        theta[i, j] <- sample(c(-0.5, 0.5), 1)
        theta[j, i] <- theta[i, j]
      }
    }
  }
  # local dense graph: 20*20
  if (type == 7) {
    theta <- matrix(0, p, p)
    # for(i in 1:p) {
    #   theta[i, i] <- sample(c(-0.5, 0, 0.5), 1)
    # }
    for (i in 1:19) {
      theta[i, i + 1] <- sample(c(-0.5, 0.5), 1)
      theta[i + 1, i] <- theta[i, i + 1]
    }
  }
  
  if (type == 8)
    theta <- gen_rr_adjmat(p, degree, beta, alpha, type = "ferro")
  if (type == 9)
    theta <- gen_rr_adjmat(p, degree, beta, alpha, type = "glass")
  if (type == 10)
    theta <- gen_4nn_cyc(p, degree, beta, alpha, type = "ferro")
  if (type == 11)
    theta <- gen_4nn_cyc(p, degree, beta, alpha, type = "glass")
  if (type == 12)
    theta <- gen_4nn_cyc(p, degree, beta, alpha, type = "glass_weak")
  
  return(theta)
}

#' Generate a frequency dataset from binary markov network
#'
#' @param n
#' @param p
#' @param type
#' @param seed
#' @param graph.seed
#' @param beta
#' @param degree
#' @param alpha
#'
#' @return
#' @export
#'
#' @examples
generate.bmn.freq.data <-
  function(n,
           p,
           type = 1,
           seed = NULL,
           graph.seed = NULL,
           beta = 0.7,
           degree = 3,
           alpha = 0.4) {
    if (is.null(graph.seed)) {
      graph_seed <- round(runif(1 , 0, .Machine$integer.max))
    }
    theta <-
      sim_theta(
        p,
        type = type,
        seed = graph_seed,
        beta = beta,
        degree = degree,
        alpha = alpha
      )
    
    if (is.null(seed)) {
      seed <- round(runif(1 , 0, .Machine$integer.max))
    }
    
    data <- sample_by_conf(n = n, theta = theta, seed = seed)
    return(list(data = data, theta = theta))
  }

#' 
#' Generate a dataset with independent observations from binary markov network
#' 
#' @param n 
#' @param p 
#' @param type 
#' @param seed 
#' 
#' 
#' @return
#' @export
#'
#' @examples
generate.bmn.data <- function(n, p, type = 1, seed = NULL, graph.seed = NULL, beta = 0.7, degree = 3, alpha = 0.4) {
  if (is.null(graph.seed)) {
    graph_seed <- round(runif(1 , 0, .Machine$integer.max))
  }
  theta <- sim_theta(p, type, graph_seed, beta, degree, alpha)
  Ising_Gibbs(theta, n, )
}
