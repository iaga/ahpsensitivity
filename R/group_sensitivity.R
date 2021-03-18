#' A function to obtain the criteria weight
#'
#' Given a pairwise comparison matrix this function allows you to obtain the
#' criteria  weight (row geometric mean)
#' for one user and any number of criteria
#'
#' @param x a pairwise comparison matrix
#' @param n Integer, the total number of criteria
#'
#' @section Details:
#' The weight criteria is obtained by the computation of the geometric mean
#' \deqn{\bar{w}_r^{P_j} = \sqrt[m]{a_{r1}a_{r2}\dots a_{rm}}}
#' where \eqn{a_{rm}} are elements in the pairwise comparisons matrix
#' \deqn{C = \left[ {\begin{array}{cccc}
#'  1 & a_{12} & \dots & a_{1m}\\
#'  1/a_{12} & 1 & \dots & a_{2m} \\
#'  \dots & \dots & 1 & \dots\\
#'  1/a_{1m} & 1/a_{2m} & \dots &  1
#'  \end{array} } \right]}
#' @references Marie Ivanco, Gene Hou, Jennifer Michaeli, Sensitivity analysis method to address user disparities in the analytic hierarchy process,
#' Expert Systems with Applications, Volume 90, 2017, Pages 111-126, ISSN 0957-4174,
#' https://doi.org/10.1016/j.eswa.2017.08.003.
#' @keywords weights matrix
#' @examples
#' pcm <- pcmatrix(c(1,2,3,4,5,6), 4)
#' weights.matrix(pcm, 4)
#' @export
weights.matrix <- function(x,n){
  nthroot(rowProds(x),n)
}

#' The normalized criteria weight
#'
#' This function allows you to obtain the normalized criteria weight
#' (normalized geometric mean) for one user
#'
#' @param x Vector of criteria weight
#'
#' @section Details: For each component \eqn{r} of the vector the normalized weight
#' is obtained by
#' \deqn{\bar{w}_r = \frac{w_r}{\sum_{q=1}^n w_q}}
#' where \eqn{n} is the total number of criteria
#' @keywords normalized weights
#' @examples
#' w <- c(.7,.2,.3)
#' norm.fun(w)
#' @return
#' A vector
#' @export
norm.fun <- function(x){
  x/sum(x)
}


to_matrix <- function(x,n,by_row){
  pert_list <- list()
  for (i in c(1:length(x))){
    pert_list[[i]] <- matrix(x[,i], nrow = n, byrow = by_row)
  }
  return(pert_list)
}

to_group_matrix <- function(x, y, n){
  # x es la matriz de grupos original
  # y es el vector que contiene la perturbación i
  group_list <- list()
  for (i in 1:length(y)){
    tmp <- x
    tmp[,n] <- y[[i]]
    group_list[[i]] <- tmp
  }
  return(group_list)
}

#' Create a group vector
#'
#' This function allows you to create the group vector
#'
#' @param x a list of matrices with the normalized weight for all criteria and users (rows: critera, columns: user)
#' @param users_total the total number of users in the group
#'
#' @section Details: The group aggregated vector is obtained by following
#' \deqn{{G}_r^{(i, j_{i-1})}= \sqrt[P]{\bar{w}_r^1\bar{w}_r^2\dots \bar{w}_r^P}}
#' where \eqn{P} represents the total number of users in the group and
#' \eqn{\bar{w}_r} is the normalized weight of the rth criterion
#' @keywords group vector
#' @examples
#' @return
#' A list of vectors
#' @references Marie Ivanco, Gene Hou, Jennifer Michaeli, Sensitivity analysis method to address user disparities in the analytic hierarchy process,
#' Expert Systems with Applications, Volume 90, 2017, Pages 111-126, ISSN 0957-4174,
#' https://doi.org/10.1016/j.eswa.2017.08.003.
#' @export
group.vector <- function (x, users_total){
  group.vector <- list()
  for (i in 1:length(x)){
    group.vector [[i]] <- rowProds(x[[i]])^(1/users_total)
  }
  return(group.vector)
}

delta_pert <- function(F, pert){
  # Data frame generation: F_orig
  colnames(F) <- c(1:nrow(F))
  rownames(F) <- c(1:nrow(F))
  data_frame <- as.data.frame(as.table(F)) %>%
    set_names(c("row", "col", "orig"))

  data_frame$perturbed <- data_frame$orig
  cont=1


  # contiene la magnitud de elementos perturbados
  perturbation_list <- list()

  for (i in c(1:nrow(F))){
    for (j in c(1:ncol(F))){
      if (i != j)  {
        # se calcula la magnitud de la perturbación para el elemento i y j
        perturbation <- (data_frame$orig[which(data_frame$row == i & data_frame$col == j)]*pert)/100
        # el valor de la perturbación se almacena en la lista perturbation_list
        perturbation_list[[cont]] <- perturbation
      }
      else if (i == j) # si son iguales,
      {
        # se asigna un valor de cero a la diagonal
        perturbation_list[[cont]] <- 0
      }
      # se incrementa el contador
      cont <- cont + 1
    }
  }
  return(perturbation_list)
}


perturbation <- function(F){
  # Data frame generation: F_orig
  colnames(F) <- c(1:nrow(F))
  rownames(F) <- c(1:nrow(F))
  data_frame <- as.data.frame(as.table(F)) %>%
    set_names(c("row", "col", "orig"))

  data_frame$perturbed <- data_frame$orig
  cont=1


  # contiene la magnitud de elementos perturbados
  perturbation_list <- list()

  for (i in c(1:nrow(F))){
    for (j in c(1:ncol(F))){
      if (i != j)  {
        # se calcula la magnitud de la perturbación para el elemento i y j
        perturbation <- (data_frame$orig[which(data_frame$row == i & data_frame$col == j)]*pert)/100
        # el valor de la perturbación se almacena en la lista perturbation_list
        perturbation_list[[cont]] <- perturbation
        # se calcula el nuevo valor perturbado para la matriz de comparaciones F
        data_frame$perturbed[which(data_frame$row == i & data_frame$col == j)] <- data_frame$orig[which(data_frame$row == i & data_frame$col == j)] + perturbation
        # se asigna el inverso para la parte de abajo de la matriz triangular
        data_frame$perturbed[which(data_frame$row == j & data_frame$col == i)] <- 1/data_frame$perturbed[which(data_frame$row == i & data_frame$col == j)]
      } else if (i == j) # si son iguales,
      {
        # se asigna un valor de cero a la diagonal
        perturbation_list[[cont]] <- 0
        # se asigna el valor original al campo perturbado
        data_frame$perturbed[which(data_frame$row == i & data_frame$col == j)] <- data_frame$orig[which(data_frame$row == i & data_frame$col == j)]
      }
      # se asigna el nombre a la columna
      colnames(data_frame)[3+cont] <- paste0("perturbed_", i, "_", j)
      # se genera otra columna que contiene al valor original de F, será de utilidad en la siguiente iteración
      data_frame$perturbed <- data_frame$orig
      # se incrementa el contador
      cont <- cont + 1
    }
  }

  # se elimina la columna con la original y la última (auxiliar)
  data_frame <- data_frame %>%
    dplyr::select(-perturbed, -orig, -row, -col)
}

diff_fin <-  function(x, delta_pert){
  diff_fin_list <- list()
  for (i in 1:length(x)){
    diff_fin_list[[i]] <- (x[[1]]-x[[i]])/delta_pert[[i]]
  }
  return(diff_fin_list)
}

sum_abs <- function(x){
  # x toma como entrada la lista que contiene el cálculo de rho (diferencias finitas) para cada usuario
  abs_diff_list <- list()
  for (i in 1:length(x)){
    abs_diff_list[[i]] <- sum(abs(x[[i]]))
  }
  # without NAs
  abs_diff_list[which(is.na(abs_diff_list) == TRUE)] <- 0
  return(abs_diff_list)
}



groupsens_data <- function(users_list, n, users_total, pert){
  group_obj <- list()
  group_obj$users_list <- users_list
  group_obj$n <- n
  group_obj$users_total <- users_total
  group_obj$pert <- pert
  class(group_obj) <- "groupsens_data"
  return(group_obj)
  #data_frames_pert <- lapply(users_list, perturbation)
  #return(diff_fin_list)
}

#' Group sensitivity
#'
#' This function allows you to identify the most influential expert on the group
#'
#' @param users_list List with the pairwise comparison matrices
#' @param pert Double, the perturbation magnitude for the finite differences method
#' @keywords ahp, group sensitivity,
#' @export
#' @examples
#' @return A list of lists for each user
#' \describe{
#'      \item{users_list}{List with the pairwise comparison matrices}
#'      \item{n}{Number of criteria}
#'      \item{users_total}{Integer, the total number of users}
#'      \item{pert}{Double, the perturbation magnitude for the finite differences method}
#'      \item{users_w}{List with weights for the pairwise comparison matrices}
#'      \item{norm_users_w}{List with the normalized weights for the pairwise comparison matrices}
#'      \item{group_matrix}{Matrix with the normalized weight for all criteria (rows: criteria, columns: users)}
#'      \item{sensitivity}{List of matrices with the sensitivity coefficients for each user}
#' }
#'
#' @references Marie Ivanco, Gene Hou, Jennifer Michaeli, Sensitivity analysis method to address user disparities in the analytic hierarchy process,
#' Expert Systems with Applications, Volume 90, 2017, Pages 111-126, ISSN 0957-4174,
#' https://doi.org/10.1016/j.eswa.2017.08.003.
group.sens <- function(users_list, pert){
  group_obj <- list()
  group_obj$users_list <- users_list
  group_obj$n <- as.numeric(nrow(users_list[[1]]))
  group_obj$users_total <- length(users_list)
  group_obj$pert <- pert
  # datos necesarios
  group_obj$users_w <- lapply(users_list, weights.matrix, group_obj$n)
  group_obj$norm_users_w <- lapply(group_obj$users_w, norm.fun)
  group_obj$data_frames_pert <- lapply(group_obj$users_list, perturbation)
  group_obj$group_matrix <- matrix(unlist(group_obj$norm_users_w), nrow = group_obj$n, byrow = FALSE)
  # Pasar de data frame a matriz, para cada uno de los usuarios
  group_obj$pert_list <- lapply(group_obj$data_frames_pert, to_matrix, group_obj$n, FALSE)
  # calcular la magnitud de la perturbación para cada elemento
  group_obj$delta_pert_list <- lapply(users_list, delta_pert, pert)
  # Lista que contiene los pesos de la matriz perturbada
  group_obj$pert_weights_list <- lapply(group_obj$pert_list, lapply, weights.matrix, group_obj$n)
  # Pesos normalizados de la matriz perturbada
  group_obj$n_pert_weights_list <- lapply(group_obj$pert_weights_list, lapply, norm.fun)
  # Matrices de grupos (perturbadas)
  for (i in c(1:group_obj$users_total)){
    group_obj$group_list[[i]] <- to_group_matrix(group_obj$group_matrix, group_obj$n_pert_weights_list[[i]], i)
    #names(group_list)[i] <- paste0("user", i, "_group_matrix_p")
  }
  # Vector de grupos
  group_obj$group_list_v <- lapply(group_obj$group_list, group.vector, group_obj$users_total)
  # Normalización del vector de grupos
  group_obj$norm_group_list_v <- lapply(group_obj$group_list_v, lapply, norm.fun)
  # diferencias finitas
  for (i in c(1:group_obj$users_total)){
    group_obj$diff_fin_list[[i]] <- diff_fin(group_obj$norm_group_list_v[[i]], group_obj$delta_pert_list[[i]])
  #   names(diff_fin_list)[i] <- paste0("user", i, "_diff_fin")
  }
  # Suma de los valores absolutos obtenidos en diferencias finitas
  group_obj$abs_diff_fin_list <- lapply(group_obj$diff_fin_list, sum_abs)
  # La lista anterior se convierte a matriz
  for (i in c(1:group_obj$users_total)){
    group_obj$sensitivity[[i]] <- matrix(unlist(group_obj$abs_diff_fin_list[[i]]), nrow = group_obj$n, byrow = TRUE)
    names(group_obj$sensitivity)[i] <- paste0("user", i, "_sensitivity")
    rownames(group_obj$sensitivity[[i]]) <- 1:group_obj$n
    colnames(group_obj$sensitivity[[i]]) <- 1:group_obj$n
  }

  class(group_obj) <- "group.sens"
  return(group_obj)
}





