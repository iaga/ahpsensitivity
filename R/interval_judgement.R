random.varying <- function(lower, upper, nmatrices){
  interval <- list()
  for (i in c(1:length(lower))){
    interval[[i]] <- runif(nmatrices, lower[i], upper[i])
  }
  return(interval)
}

#' A function to obtain the consistency index
#'
#' @param lambda_max largest eigenvalue
#' @param n total number of criteria
#'
#' @section Details:
#' \deqn{C.I. = \frac{\lambda_{max}-n}{n-1}} where
#' \eqn{\lambda_{max}} is the largest eigenvalue and \eqn{n} represents the total number of
#' criteria
#' @keywords consistency index, largest eigenvalue
#' @examples
#' consistency.index(lambda_max = 1, n = 6)
#' @return
#' @export
consistency.index <- function(lambda_max, n){
  ci <-  (lambda_max-n)/(n-1)
}


#' A function to obtain the consistency ratio
#'
#'
#'
#' @param ci Consistency index
#' @param ri Random inconsistency
#'
#' @section Details:
#' \deqn{C.R. = \frac{C.I.}{R.I.}}
#' @note
#' The random inconsistency number is a function of the number of criteria
#' @keywords
#' @examples
#' consistency.ratio(.05, 0.90)
#' @return
#' @export
consistency.ratio <- function(ci, ri){
  cr <- ci/ri
}


#' A function to create a pairwise comparison matrix
#'
#'
#'
#' @param x Vector, it contains the upper triangular elements of the matrix
#' @param n Integer, total number of criteria
#'
#' @section Details:
#' Given a vector with the upper triangular elements of a pairwise comparison matrix,
#' this function returns a matrix in the pairwise comparison form
#' \deqn{C = \left[ {\begin{array}{cccc}
#'  1 & a_{12} & \dots & a_{1m}\\
#'  1/a_{12} & 1 & \dots & a_{2m} \\
#'  \dots & \dots & 1 & \dots\\
#'  1/a_{1m} & 1/a_{2m} & \dots &  1
#'  \end{array} } \right]}
#' @keywords pairwise comparison matrix
#' @examples
#' pcmatrix(c(1,2,3,4,5,6), 4)
#' @return
#' @export
pcmatrix <- function(x, n){
  # a matrix with zeros in the diagonal
  mat_mesa1 <- matrix(rep(0,n*n), nrow=n)
  # upper triangular matrix elements (by column)
  mat_mesa1[upper.tri(mat_mesa1)] <- x
  # transpose matrix
  mat_mesa1_t <- t(mat_mesa1)
  # lower triangular elements
  mat_mesa1[lower.tri(mat_mesa1)] <- 1/mat_mesa1_t[lower.tri(mat_mesa1_t)]
  # diagonal elements to one
  diag(mat_mesa1) <- rep(1,n)
  return(mat_mesa1)
}

lambda <- function(matrices, w){
  lambda <- (matrices %*% matrix(w, ncol=1))/w
}


kolmogorov <- function(w, promedio, desviacion){
  end <- ncol(w)-1
  c <- 0
  ks_test <- list()
  for (i in c(1:end)){
    ks_test[[i]] <- ks.test(w[i], pnorm, promedio[i], desviacion[i])
    if (ks_test[[i]]$p.value < 0.05){
      c <- c + 1
    }
  }
  if (c > 0){
    print("Kolmogorov-Smirnov test failed...")
  }
  else {
    print("Normality test passed")
  }
  return(ks_test)
}


lower.bound <- function(promedio, s){
  lower_limit <- promedio-(2.58*s)
}


upper.bound <- function(promedio, s){
  upper_limit <- promedio+(2.58*s)
}

interval.judgement <- function(data_lower, data_upper, nmatrices, norm_test){
  obj_interval <- list()
  # upper and lower intervals
  obj_interval$lower <- data_lower[upper.tri(data_lower)]
  obj_interval$upper <- data_upper[upper.tri(data_upper)]



  # random matrices (list)
  interval <- random.varying(obj_interval$lower, obj_interval$upper, nmatrices)

  # list to array
  interval_data <- t(simplify2array(interval))

  # total number of criteria
  n= ncol(data_lower)
  matrices_list <- list()
  obj_interval$matrices_list <- lapply(as.data.frame(interval_data), pcmatrix, n)

  # Weight calculation
  obj_interval$matrices_w <- lapply(obj_interval$matrices_list, weights_matrix, n)
  # Normalized weight
  obj_interval$norm_matrices_w <- lapply(obj_interval$matrices_w, norm_fun)

  # Inconsistency
  lambda <- mapply(lambda, obj_interval$matrices_list, obj_interval$norm_matrices_w)
  obj_interval$lambda_max <- apply(lambda, 2, mean)
  obj_interval$ci <- consistency.index(obj_interval$lambda_max,n)
  ri <- ri_table$ri[n]
  obj_interval$cr <- consistency.ratio(obj_interval$ci, ri)

  # the inconsistency matrix are transformed to NA
  df <- as.data.frame.list(obj_interval$norm_matrices_w)
  df <- as.data.frame(t(df))
  df$cr <- obj_interval$cr
  df$cr[which(df$cr > .1)] <- NA
  # Na rows are eliminated
  obj_interval$w_consistent <-  df[complete.cases(df), ]

  # Interval calculation: stadistics
  obj_interval$maximum = apply(obj_interval$w_consistent[,1:n], 2, max)
  obj_interval$minimum = apply(obj_interval$w_consistent[,1:n], 2, min)
  obj_interval$promedio = apply(obj_interval$w_consistent[,1:n], 2, mean)
  obj_interval$desviacion = apply(obj_interval$w_consistent[,1:n], 2, sd)


  # Kolomgorov-Smirnoff test

  if (norm_test == TRUE){
    obj_interval$normality <- kolmogorov(obj_interval$w_consistent, obj_interval$promedio, obj_interval$desviacion)
  }



  # Auxiliar dataframe
  cat = 1:n
  data_row <- data.frame(cat, obj_interval$promedio, obj_interval$desviacion) %>%
    set_names(c("cat_row", "promedio_row", "desviacion_row"))
  data_row$lower_limit_row <- lower.bound(obj_interval$promedio, obj_interval$desviacion)
  data_row$upper_limit_row <- upper.bound(obj_interval$promedio, obj_interval$desviacion)
  data_col <- data.frame(cat, obj_interval$promedio, obj_interval$desviacion) %>%
    set_names(c("cat_col", "promedio_col", "desviacion_col"))
  data_col$lower_limit_col <- lower.bound(obj_interval$promedio, obj_interval$desviacion)
  data_col$upper_limit_col <- upper.bound(obj_interval$promedio, obj_interval$desviacion)


  # Elements in the upper triangular
  m <- matrix(1:(n*n),n,n)
  m_triang <- upper.tri(m)
  indices <- data.frame (ind = which (m_triang == TRUE, arr.ind = TRUE))%>%
    set_names(c("cat_row", "cat_col"))

  # Df that relates the indices with the input data
  data_row <- merge(indices, data_row, by= "cat_row")
  data_col <- merge(indices, data_col, by= "cat_col")
  data_row <- data_row[order(data_row$cat_row, data_row$cat_col),]
  data_col <- data_col[order(data_col$cat_row, data_col$cat_col),]
  data <- merge(data_row, data_col, by= c("cat_row", "cat_col"))

  # Interval judgment condition check (4 conditions)
  data$cond_1 <- ifelse((data$lower_limit_col <= data$lower_limit_row)
                        & (data$upper_limit_row<=data$upper_limit_col), 1, 0)
  data$cond_2 <- ifelse((data$lower_limit_row < data$lower_limit_col)
                        & (data$upper_limit_col < data$upper_limit_row), 2, 0)
  data$cond_3 <- ifelse((data$lower_limit_row < data$lower_limit_col)
                        & (data$lower_limit_col < data$upper_limit_row)
                        & (data$upper_limit_row < data$upper_limit_col), 3, 0)
  data$cond_4 <- ifelse((data$lower_limit_col < data$lower_limit_row)
                        & (data$lower_limit_row < data$upper_limit_col)
                        & (data$upper_limit_col < data$upper_limit_row), 4, 0)
  data$cond <- round(rowSums(data[11:14]))

  # Auxiliar df with zeros
  data$one_minus_p1 <- zeros(nrow(data))
  data$one_minus_p2 <- zeros(nrow(data))
  data$one_minus_p3 <- zeros(nrow(data))
  data$one_minus_p4 <- zeros(nrow(data))

  # Df with probabilities of rr
  by_out <- data %>% group_by(cond)
  data_prob <- data.frame(data$cat_row, data$cat_col, data$cond)%>%
    set_names(c("cat_row", "cat_col", "cond"))

  # **Probability of rr condition 1:**
  data_1 <- subset(data, cond == 1)
  if (sum(data$cond_1) != 0){
    data_1$F_row_upper <- zeros(nrow(data_1))
    data_1$F_row_lower <- zeros(nrow(data_1))
    #pnorm, grado de libertad
    data_1$F_col_upper <- mapply(pnorm, q=data_1$upper_limit_row, mean=data_1$promedio_col,
                                 sd = data_1$desviacion_col, lower.tail = TRUE)
    data_1$F_col_lower <- mapply(pnorm, q=data_1$lower_limit_row, mean=data_1$promedio_col,
                                 sd = data_1$desviacion_col, lower.tail = TRUE)
    data_1$dif_1 <- data_1$F_col_upper - data_1$F_col_lower
    data_prob <- merge(x = data_prob, y = data_1[ , c("cat_row", "cat_col", "dif_1")],
                       by=c("cat_row", "cat_col"), all.x=TRUE)
  }



  # **Probability of rr condition 2:**
  data_2 <- subset(data, cond == 2)
  if (sum(data$cond_2) != 0){
    data_2$F_row_upper <- zeros(nrow(data_2))
    data_2$F_row_lower <- zeros(nrow(data_2))
    data_2$F_col_upper <- mapply(pnorm, q=data_2$upper_limit_col, mean=data_2$promedio_row,
                                 sd = data_2$desviacion_row, lower.tail = TRUE)
    data_2$F_col_lower <- mapply(pnorm, q=data_2$lower_limit_col, mean=data_2$promedio_row,
                                 sd = data_2$desviacion_row, lower.tail = TRUE)
    data_2$dif_2 <- data_2$F_col_upper- data_2$F_col_lower
    data_prob <- merge(x = data_prob, y = data_2[ , c("cat_row", "cat_col", "dif_2")],
                       by=c("cat_row", "cat_col"), all.x=TRUE)
  }



  # **Probability of rr condition 3:**
  data_3 <- subset(data, cond == 3)
  if (sum(data$cond_3) != 0){
    data_3$F_row_upper <- mapply(pnorm, q=data_3$upper_limit_row, mean=data_3$promedio_row,
                                 sd = data_3$desviacion_row, lower.tail = TRUE)
    data_3$F_row_lower <- mapply(pnorm, q=data_3$lower_limit_col, mean=data_3$promedio_row,
                                 sd = data_3$desviacion_row, lower.tail = TRUE)
    data_3$F_col_upper <- mapply(pnorm, q=data_3$upper_limit_row, mean=data_3$promedio_col,
                                 sd = data_3$desviacion_col, lower.tail = TRUE)
    data_3$F_col_lower <- mapply(pnorm, q=data_3$lower_limit_col, mean=data_3$promedio_col,
                                 sd = data_3$desviacion_col, lower.tail = TRUE)
    data_3$dif_3 <- (data_3$F_col_upper
                     -data_3$F_col_lower)*(data_3$F_col_upper- data_3$F_col_lower)
    data_prob <- merge(x = data_prob, y = data_3[ , c("cat_row", "cat_col", "dif_3")],
                       by=c("cat_row", "cat_col"), all.x=TRUE)
  }



  # **Probability of rr condition 4:**
  data_4 <- subset(data, cond == 4)
  if (sum(data$cond_4) != 0){
    data_4$F_row_upper <- mapply(pnorm, q=data_4$upper_limit_col, mean=data_4$promedio_row,
                                 sd = data_4$desviacion_row, lower.tail = TRUE)
    data_4$F_row_lower <- mapply(pnorm, q=data_4$lower_limit_row, mean=data_4$promedio_row,
                                 sd = data_4$desviacion_row, lower.tail = TRUE)
    data_4$F_col_upper <- mapply(pnorm, q=data_4$upper_limit_col, mean=data_4$promedio_col,
                                 sd = data_4$desviacion_col, lower.tail = TRUE)
    data_4$F_col_lower <- mapply(pnorm, q=data_4$lower_limit_row, mean=data_4$promedio_col,
                                 sd = data_4$desviacion_col, lower.tail = TRUE)
    data_4$dif_4 <- (data_4$F_row_upper
                     -data_4$F_row_lower) * (data_4$F_col_upper- data_4$F_col_lower)
    data_prob <- merge(x = data_prob, y = data_4[ , c("cat_row", "cat_col", "dif_4")],
                       by=c("cat_row", "cat_col"), all.x=TRUE)
  }


  # All the probabilites in one df
  data_prob[is.na(data_prob)] <- 0
  obj_interval$data <- data
  obj_interval$data_prob <- data_prob
  end <- length(data_prob)
  data_prob$p_ij <- rowSums(data_prob[4:end])
  obj_interval$p_ij <- round(data_prob$p_ij, digits = 3)
  data_prob$one_minus_p_ij <- 1-data_prob$p_ij
  obj_interval$ones_minus_pij <- round(data_prob$one_minus_p_ij, digits = 3)

  #Probabilities
  small <- matrix(rep(1,n*n), nrow=n)
  small[lower.tri(small)] <- data_prob$one_minus_p_ij
  small <- t(small)
  small[lower.tri(small)] <- data_prob$one_minus_p_ij
  obj_interval$p <- 1-colProds(small)
  obj_interval$p <- round(obj_interval$p, digits = 3)

  class(obj_interval) <- "interval.judgement"
  return(obj_interval)
}



#' A function to obtain random matrices
#'
#' @param data_lower A list of dataframes objects with the lower bounds for each pairwise comparison
#' matrix
#' @param data_upper A list of dataframes objects with the upper bounds for each pairwise comparison
#' matrix
#' @param nmatrices Integer, total number of random matrices to generate
#'
#' @section Details:
#' Given matrices
#' \deqn{C^{L} = \left[ {\begin{array}{cccc}
#'  1 & a^L_{12} & \dots & a^L_{1m}\\
#'  1/a^L_{12} & 1 & \dots & a^L_{2m} \\
#'  \dots & \dots & 1 & \dots\\
#'  1/a^L_{1m} & 1/a^L_{2m} & \dots &  1
#'  \end{array} } \right],\hspace{.5cm}
#'  C^{U} = \left[ {\begin{array}{cccc}
#'  1 & a^U_{12} & \dots & a^U_{1m}\\
#'  1/a^U_{12} & 1 & \dots & a^U_{2m} \\
#'  \dots & \dots & 1 & \dots\\
#'  1/a^U_{1m} & 1/a^U_{2m} & \dots &  1
#'  \end{array} } \right]}
#' with the lower and upper bounds of an interval pairwise comparison
#' matrix, this function generates a set of random matrices within the interval. The random
#' numbers for each comparison are generated with an uniform distribution.
#' @seealso
#' runif about the uniform distribution
#' @references Saaty, Thomas L. y Vargas, Luis G.: Uncertainty and rank order in the analytic hierarchy
#'process. European Journal of Operational Research, 1987, 32(1), pp. 107–117. ISSN 0377-
#'  2217. doi: 10.1016/0377-2217(87)90275-X.
#' @keywords random numbers, matrix
#' @return A list of matrices
#' @examples
#' @export
random.pcmatrices <- function(data_lower, data_upper, nmatrices){
  interval_l <- list()
  # upper and lower intervals
  interval_l$lower <- data_lower[upper.tri(data_lower)]
  interval_l$upper <- data_upper[upper.tri(data_upper)]

  # random matrices (list)
  interval <- random.varying(interval_l$lower, interval_l$upper, nmatrices)

  # list to array
  interval_data <- t(simplify2array(interval))

  interval_data <- t(simplify2array(interval))
  n= ncol(data_lower)
  matrices_list <- list()
  interval_l$matrices_list <- lapply(as.data.frame(interval_data), pcmatrix, n)
  return(interval_l)
}


#' A function to obtain the probabilities of rank reversals
#'
#'
#'
#' @param lower A list of dataframes objects with the lower bounds for each pairwise comparison matrix
#' @param upper A list of dataframes objects with the upper bounds for each pairwise comparison matrix
#' @param nmatrices Integer, defines the total number of random matrices to generate
#' @param norm_test Logical, if TRUE determines if the components of the right eigenvector of all the generated random
#' pairwise comparison matrices are normally distributed, via the kolmogorov-smirnoff test
#'
#' @section Details:
#' The probabilities \eqn{p_i}, \eqn{i = 1, \dots , n} that a given alternative will reverse
#' rank with another alternative are given by \deqn{p_i=1-\prod_{j=1}^n (1-p_{ij})} where
#' \eqn{p_{ij}} is the probability of rank reversal for two alternatives, it is calculated considering
#' different cases for the lower and upper bounds on the \eqn{i} and \eqn{j} components of the right eigenvector \eqn{w} and the
#' probability cumulative distribution \eqn{F_i(x_i)}:
#' \tabular{lll}{
#' Case \tab Condition \tab \eqn{p_{ij}} \cr
#' 1 \tab \eqn{w_j^L\leq w_i^L} and \eqn{w_i^U\leq w_j^U} \tab \eqn{F_j(w_i^U)-F_j(w_i^L)} \cr
#' 2 \tab \eqn{w_i^L < w_j^L} and \eqn{w_j^U < w_i^U} \tab \eqn{F_i(w_j^U)-F_i(w_j^U)} \cr
#' 3 \tab \eqn{w_i^L < w_j^L < w_i^U < w_j^U} \tab \eqn{(F_i(w_i^U)-F_i(w_j^L))(F_j(w_i^U)-F_j(w_j^L))}\cr
#' 4 \tab \eqn{w_j^L < w_i^L < w_j^U < w_i^U} \tab \eqn{(F_i(w_j^U)-F_i(w_i^L))(F_j(w_j^U)-F_j(w_i^L))} \cr
#' }
#'
#' @references Saaty, Thomas L. y Vargas, Luis G.: Uncertainty and rank order in the analytic hierarchy
#'process. European Journal of Operational Research, 1987, 32(1), pp. 107–117. ISSN 0377-
#'  2217. doi: 10.1016/0377-2217(87)90275-X.
#' @keywords rank reversal, interval pairwise comparison matrix
#' @examples
#'
#' @return A list of list, for each interval pairwise comparison matrix
#' \describe{
#'   \item{matrix_list}{List of random matrices (rm)}
#'   \item{matrices_w}{List of weight of rm}
#'   \item{norm_matrices_w}{List of normalized weights of rm}
#'   \item{lambda_max}{List of largest eigenvalues}
#'   \item{ci}{List of consistency indices}
#'   \item{cr}{List of consistency ratios}
#'   \item{w_consistent}{Dataframe of normalized weights of consistent random matrices (crm)}
#'   \item{maximum}{Named vector, maximum value of the crm for each criteria}
#'   \item{minimum}{Named vector, minimum value of the crm for each criteria}
#'   \item{mean}{Named vector, mean value of the crm for each criteria}
#'   \item{sd}{Named vector, standard deviation value of the crm for each criteria}
#'   \item{normality}{List of list with the normality test for each criteria}
#'   \item{p_ij}{Vector, probabilities of rank reversal between two alternatives \eqn{A_i} and \eqn{A_j}}
#'   \item{p_i}{Vector, probabilities that a given alternative will reverse rank with other alternative}
#' }
#' @export
ahp.interval <- function(lower, upper, nmatrices, norm_test){
  out_list <- list()
  for (i in c(1:(length(n2_lower)))){
    out_list[[i]] <- interval.judgement(lower[[i]], upper[[i]], nmatrices, norm_test = TRUE)
    names(out_list)[i] <- paste0("matrix_", i)
  }
  return(out_list)
}

