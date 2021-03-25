wsmtau <- function(a_ij, w_j, indice, indice_rho){
  obj_tau <- list()
  # se calcula la combinación lineal ponderada
  # indice <- a_ij%*%(w_j)
  # se calcula el valor de referencia
  # indice_rho <- median(indice)
  names_aij <- colnames(a_ij)
  num_col_w <- length(w_j)
  #ones <- array(rep(1, (nrow(a_ij)*num_col_w)), dim=(c(nrow(a_ij[[1]]), num_col_w)))

  #*******chunk 1******************************
  # Se coloca un indicador de la fila y la columna
  colnames(a_ij) <- c(1:ncol(a_ij))
  rownames(a_ij) <- c(1:nrow(a_ij))

  # se genera el data frame
  aij_data_frame <- as.data.frame(as.table(a_ij)) %>%
    set_names(c("row", "col", "aij_orig"))

  # se genera el data frame con los valores normalizados repetidos
  list_a_ij_r <- list()
  for (i in 1:length(w_j)) {
    list_a_ij_r[[i]] <- rep(a_ij[,i],length(w_j))
  }
  a_ij_rep <- data.frame(matrix(unlist(list_a_ij_r), nrow=(length(a_ij)), byrow=F),stringsAsFactors=FALSE)

  for (i in 1:length(w_j)) {
    colnames(a_ij_rep)[i] <- paste0(names_aij[i])
  }
  #*******chunk 2******************************
  nw <- length(w_j)
  na <- nrow(a_ij)

  # se crea una lista con el peso del criterio m_i repetido n veces, es decir: una lista de m elementos con n subelementos (15 elementos de 484 valores repetidos para el censo de 1998)

  list_w <- list()
  for (i in 1:length(w_j)) {
    list_w[[i]] <- rep(w_j[i],na)
  }

  # la lista se transforma en vector
  bd_w <- unlist(list_w)

  # se crea una lista con el valor del indice n_i repetido m veces, es decir: una lista de n elementos con m subelementos (484 elementos de 15 valores repetidos para el censo de 1998)

  list_indice <- list()
  for (i in 1:length(indice)) {
    list_indice[[i]] <- rep(indice[i],nw)
  }

  # la lista se transforma en vector
  bd_indice <- unlist(list_indice)

  #*******Cálculo de tau: denominador******************************
  # Se crea una lista con la multiplicación de los pesos por los valores de la encuesta. La cual está compuesta por m sublistas de n elementos (15 listas de 484 elementos = 7260 para capacidad adaptativa)

  list_wa <- list()
  for (j in 1:length(w_j)) {
    list_wa[[j]] <- w_j[j]*ones[,j]
  }
  # la lista se transforma en vector
  bd_list_wa <- unlist(list_wa)


  # Con el objetivo de convertir la lista "lista_wa" en base de datos, se crean indicadores de columna en forma vectorial
  wa_col <- list()
  for (i in 1:length(w_j)) {
    wa_col[[i]] <- rep(i,na)
  }
  # la lista se transforma en vector
  wa_col_v <- unlist(wa_col)

  # Se crea el data frame que contiene el indicador fila, el denominador w_j a_ij y el indicador columna
  wa_data_frame <- as.data.frame(as.table(as.matrix(bd_list_wa)))%>%
    set_names("row","VA","wa")%>%
    mutate(col = wa_col_v)%>%
    dplyr::select(-VA)

  #*******Cálculo de tau: numerador******************************
  # Se crea una lista con la diferencia entre el índice para cada encuesta menos el valor de rho
  list_dif <- list()

  for (i in 1:length(indice)) {
    list_dif[[i]] <- indice[i,1]-indice_rho
  }

  # se repite para el número total de criterios
  list_dif_comp <- rep(list_dif, nw)

  # la lista se convierte en vector de dimensión mxn
  list_dif_v <- unlist(list_dif_comp)


  # se genera el data frame que contiene el número de fila y la diferencia entre los índices
  indice_dif_data_frame <- as.data.frame(as.table(unlist(list_dif_v)))%>%
    set_names("row", "indice_dif")

  bd_indice_dif <- data_frame(aij_data_frame$row, indice_dif_data_frame$indice_dif)%>%
    set_names("row", "indice_dif")
  bd_indice_dif <- as.data.frame(bd_indice_dif)
  #*************************************
  # Se genera un nuevo data frame, con las combinaciones necesarias para realizar el cálculo

  # Contiene el producto wa para cada uno de los criterios

  # m listas de n elementos (15 listas de 484 elementos)
  list_crit <- list()
  for (i in 1:length(w_j)) {
    list_crit[[i]] <- subset(wa_data_frame, wa_data_frame$col == i, select = wa)
  }
  # se transforma a vector
  crit_wa <- unlist(list_crit)

  # lista con el indicador de columna
  crit_col <- list()
  for (i in 1:length(w_j)) {
    crit_col[[i]] <- rep(i,na)
  }

  crit_col_v <- unlist(crit_col)

  # base de datos que contiene el indicador de fila, el producto wa y el indicador de fila
  bd_crit <- data.frame(aij_data_frame$row, crit_wa)%>%
    set_names("row", "wa")%>%
    mutate(col = crit_col_v)


  cont=3

  # dataframe que contiene el producto wa para cada uno de los criterios
  for (i in c(1:length(w_j))){
    bd_indice_dif$criteria <-   bd_crit$wa[which(bd_crit$col == i #& bd_crit$col == j)
    )]
    #colnames(bd_indice_dif)[cont] <- paste0("criteria_", i)
    colnames(bd_indice_dif)[cont] <- paste0(names_aij[i])
    cont <- cont + 1
  }

  bd_indice_dif <- as.data.frame(bd_indice_dif)
  # asignación de un valor pequeño, diferente de cero, a la diferencia entre V_i y V_rho
  # bd_indice_dif[which(bd_indice_dif[,2] == 0),2] <- 0.001

  #*************************************
  # Se aplica la fórmula de $\tau$ a todas las combinaciones del vector de prioridades
  bd_tau <- data.frame(bd_indice_dif$row, bd_indice_dif$indice_dif)%>%
    set_names("row", "indice_dif")

  for (i in c(3:(length(w_j)+2))) {
    bd_tau$tau <- bd_indice_dif[,2]/bd_indice_dif[,i]
    #colnames(bd_tau)[i] <- paste0("tau_", i-2)
    colnames(bd_tau)[i] <- paste0(colnames(bd_indice_dif)[i])
  }

  #*************************************
  # Generación del rango de factibilidad
  # Se calcula la cota por abajo sobre tau, a partir de $a_{ij}$
  tau_prima_abajo <- a_ij_rep-1
  tau_prima_arriba <- a_ij_rep

  #*************************************
  # Verificación del rango de factibilidad "feasibility range"
  # Si está fuera del intervalo, el valor se hace cero (para facilitar el cálculo del grado de criticalidad)
  for (i in c(3:(length(w_j)+2))){
    bd_tau[which(bd_tau[,i] < tau_prima_abajo[,i-2]),i] <- NA
    bd_tau[which(bd_tau[,i] > tau_prima_arriba[,i-2]),i] <- NA
  }

  bd_tau <- bd_tau[1:nrow(a_ij),]

  # Cálculo del valor para que exista una reversión de rango
  a_ij_rev <- a_ij-bd_tau[,3:(nrow(pesos)+2)]

  #*************************************
  # Cálculo del grado de criticalidad

  # Se crea un data frame en blanco
  criticality_degree_quart = data.frame(matrix(ncol=0,nrow=0))

  # Cálculo del grado de criticalidad
  for (i in c(3:(length(w_j)+2))) {
    criticality_degree_quart[1,i-2] <- 1/summary(abs(bd_tau[,i]), na.rm = TRUE)[2]
  }

  # Las columnas del data frame se nombran a partir de la base bd_tau
  tau_names <- names(bd_tau[c(3:(length(w_j)+2))])
  colnames(criticality_degree_quart) <- tau_names

  #*************************************
  # Se obtiene la proporción de encuestas (con indicador diferente a NA) para cada indicador
  # Se crea un data frame en blanco
  prop_encuestas = data.frame(matrix(ncol=0,nrow=0))

  # Cálculo del grado de criticalidad
  for (i in c(3:(length(w_j)+2))) {
    prop_encuestas[1,i-2] <- (length(bd_tau[!is.na(bd_tau[,i]),i])/length(bd_tau[,i]))
  }

  # Las columnas del data frame se nombran a partir de la base bd_tau
  colnames(prop_encuestas) <- tau_names

  #*************************************
  # Indicadores más criticos
  # El producto del grado de criticalidad por la proporción de encuestas
  criticality <- round(criticality_degree_quart*(prop_encuestas),2)
  colnames(criticality) <- tau_names

  # El data frame se convierte a un vector
  rank_indicator_vector <- unlist(criticality)


  #*************************************
  # Data frame de salida
  dt <- data.frame(t(round(criticality_degree_quart, 2)))
  names(dt)[1] <- "sens_coeff"
  dt$p <- t(round(prop_encuestas, 2))
  dt$result <- (rank_indicator_vector)
  dt <- dt[order(dt$result, decreasing = TRUE),]
  # se elimina la fila con NAs
  dt <- na.omit(dt)

  obj_tau$tau <- bd_tau
  obj_tau$a_ij_hat <- a_ij_rev
  obj_tau$sensitivity <- dt

  class(obj_tau) <- "wsmtau"
  return(obj_tau)
}

#' A function to obtain the most critical measure of performance in terms of some
#' reference value for a set of indices
#'
#'
#'
#' @param x_ij A list of matrices with the standarized measures of performance
#' for each index \eqn{i} and each criteria \eqn{j}
#' @param w_ij A list of matrices with the importance weights for each criteria \eqn{j}
#' @param indices A list of dataframes with the preferences for each alternative
#' @param rho A list with the references values for each index
#'
#' @section Details:
#' For each index \eqn{i}, the minimum change in a measure of performance \eqn{x^h_{ij}} for a criteria \eqn{j} that
#' generates a rank reversal is obtained with
#' \deqn{\tau_{ij}^h = \frac{V_i^{h}-V_i^{\rho}}{w_{ij}x^h_{ij}}}
#' In terms of \eqn{\tau_{ij}^h}, the critical indicator value
#' \deqn{\mathcal{C}_{ij} = \frac{1}{\Delta_{ij}}\times p_{ij}}
#' considers the sensitivity coefficient \eqn{\Delta_{ij}=|\tau_{ij}^h|^{Q_1}}
#' (first quartile \eqn{Q_1}) and the probability of rank reversals \eqn{p_{ij}}.
#' Moreover, the modified measure of performance is given by
#' \deqn{\hat{x}_{ij}^h = {x}_{ij}^h - \tau_{ij}^h}
#' @keywords
#' @return A list of lists for each index
#' \describe{
#'      \item{tau}{A dataframe with the threshold value \eqn{\tau} for each criteria}
#'      \item{x_ij_hat}{A dataframe with the modified measure of performance value \eqn{\hat x_{ij}}
#'      \item{sensitivity}{A list with components: summarizes the sensitivity coefficient,
#'      probability of rank reversal and critical indicator value for all the
#'      measures of performance with the first quartile}
#'      }
#' @examples
#' @references
#' Sensitivity analysis for household vulnerability assessment: a case of study from Brazil surveys
#'
#' Triantaphyllou, Evangelos y Sánchez, Alfonso: A Sensitivity Analysis Approach for Some Deterministic
#' Multi-Criteria Decision-Making Methods*. Decision Sciences, 1997, 28(1), pp. 151–194.
#' doi: 10.1111/j.1540-5915.1997.tb01306.x.
#' https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1540-5915.1997.tb01306.x
#' @return
#' @export
wsm.tau <- function(x_ij, w_ij, indices, rho){
  out_list <- list()
  for (i in c(1:(length(indices)))){
    out_list[[i]] <- wsmtau(x_ij[[i]], w_ij[[i]], indices[[i]], rho[[i]])
    #names(out_list)[i] <- paste0("fp_", i)
    }
    return(out_list)
 }

