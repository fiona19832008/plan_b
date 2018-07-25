# Function to estimate mu and theta
scRNA_cluster_glm <- function(Dmatrix, n) {
  require(MASS)
  # Dmatrix: single cell RNA seq (scRNA-seq) data matrix with genes as rows and cells as columnes
  # n: threshold for removing genes with less than n cells have reading
  
  # Remove genes with less than n cells having reading
  no0 <- Dmatrix>0
  genesum <- c()
  for (i in 1:nrow(no0)) {
    genesum[i] = sum(no0[i,])
  }
  N_more <- which(genesum >= n)
  Data_N_more <- Dmatrix[N_more,]
  rownames(Data_N_more) <- paste0(N_more, "_", rownames(Dmatrix[N_more,]))
  # Fit negative binomial distribution and estimate mu and overdispersion
  UMI <- colSums(Data_N_more)
  stat <- matrix(ncol = 6, nrow = nrow(Data_N_more), byrow = TRUE)
  mu <- matrix(ncol = ncol(Data_N_more), nrow = nrow(Data_N_more), byrow = TRUE)
  for (i in 1:nrow(Data_N_more)){
    try({modelNB <- glm.nb(Data_N_more[i,]~1+offset(log(UMI)))})
    stat[i,] <- c(as.numeric(modelNB$theta), as.numeric(modelNB$coefficients),
                  as.numeric(modelNB$converged), as.numeric(modelNB$SE.theta),
                  summary(modelNB)$coefficient[2], as.numeric(modelNB$twologlik))
    mu[i,] <- as.numeric(modelNB$fitted.values)
  }
  colnames(stat) <- c("theta", "intercept", "converged", "theta SE", "intercept SE", "2loglike")
  rownames(stat) <- rownames(Data_N_more)
  colnames(mu) <- colnames(Data_N_more)
  rownames(mu) <- rownames(Data_N_more)
  stat <- stat[which(stat[,3]==1), ]
  mu <- mu[which(stat[,3]==1), ]
  # Return the following matrixes
  # Data_N_more: data matrix removed genes with less than n cells having counts
  # stat: with theta and some fitting statistics from negative binomial distribution of eath gene
  # mu: fitted mu of each gene in each cell
  return(list(Data_N_more, stat, mu))
}

# Function to generate the weight matrix and Xij matrix
scRNA_cluster_taylor <- function(yij, mu_ij, theta_j, I) {
  # I is indicator for total UMI to be 1 or colSum
  denominator <- matrix(nrow = nrow(yij), ncol = ncol(yij))
  for (i in 1:ncol(yij)) {
    denominator[,i] <- exp(mu_ij)[,i]+theta_j
  }
  # alpha ij
  alpha <- exp(mu_ij)/denominator
  mu_ij_theta <- matrix(nrow = nrow(yij), ncol = ncol(yij))
  for (i in 1:ncol(yij)) {
    mu_ij_theta[,i] <- exp(mu_ij[,i])*theta_j
  }
  # r ij
  r_ij <- mu_ij_theta/((denominator)^2)
  # Get weight matrix W
  y_theta <- matrix(nrow = nrow(yij), ncol = ncol(yij))
  for (i in 1:ncol(yij)) {
    y_theta[,i] <- yij[,i]+theta_j
  }
  W <- (0.5)*y_theta*r_ij
  #W_1 <- sqrt(W)
  # X ij matrix
  if (I==1) {
    UMI <- colSums(yij)
  } else if (I==0) {
    UMI <- rep(1, ncol(yij))
  } else if (I==2) {
    UMI <- colSums(yij)/100
  }
  
  #UMI <- rep(1, ncol(yij))
  mu_ij_logUMI <- matrix(nrow = nrow(yij), ncol = ncol(yij), byrow = TRUE)
  for (i in 1:nrow(yij)) {
    mu_ij_logUMI[i,] <- mu_ij[i,]-log(UMI)
  }
  X_ij <- mu_ij_logUMI+((yij-alpha*y_theta)/(r_ij*y_theta))
  return(list(X_ij, W, UMI))
}

# Wrapper function to run the algarithm
scRNA_cluster <- function(DM, n, I) {
  # DM: data to be analyzed
  # n: threshold for removing genes with less than n cells have reading
  # k: Number of principal components to retrieve.
  # t: Number of iterations (allows convergence to the dominant principal component)
  DM <- as.matrix(DM)
  DM_1 <- scRNA_cluster_glm(DM, n)
  DM_W_X <- scRNA_cluster_taylor(DM_1[[1]], log(DM_1[[3]]), as.numeric(DM_1[[2]][,1]), I)
  #DM_Center <- PC_center(DM_W_X[[1]])
  #DM_WPCA <- wPCA(DM_W_X[[1]], DM_W_X[[2]])
  return(list(DM_1[[1]], DM_1[[2]], DM_1[[3]], 
              DM_W_X[[1]], DM_W_X[[2]], DM_W_X[[3]])) #DM_WPCA[[1]], DM_WPCA[[2]], DM_WPCA[[4]], DM_WPCA[[5]],
  # return [[1]]: Data_n_more
  # return [[2]]: fitting stat from Data_n_more: dispersion is [,1]
  # return [[3]]: fitted mu from Data_n_more
  # return [[4]]: Xij matrix
  # return [[5]]: Weight matrix
  # return [[6]]: UMI
} 

# Order theta after glm fitting before calculate Xij and Weight matrix
scRNA_cluster_OrderTheta <- function(DM, n, I) {
  # DM: data to be analyzed
  # n: threshold for removing genes with less than n cells have reading
  # k: Number of principal components to retrieve.
  # t: Number of iterations (allows convergence to the dominant principal component)
  DM <- as.matrix(DM)
  DM_1 <- scRNA_cluster_glm(DM, n)
  DM_input <- DM_1[[1]][order(as.numeric(DM_1[[2]][,1])), ]
  mu_input <- DM_1[[3]][order(as.numeric(DM_1[[2]][,1])), ]
  theta_input <- DM_1[[2]][order(as.numeric(DM_1[[2]][,1])), ]
  DM_W_X_original <- scRNA_cluster_taylor(DM_1[[1]], log(DM_1[[3]]), as.numeric(DM_1[[2]][,1]), I)
  DM_W_X <- scRNA_cluster_taylor(DM_input, log(mu_input), as.numeric(theta_input[,1]), I)
  return(list(DM_1[[1]], DM_1[[2]], DM_1[[3]], 
              DM_W_X[[1]], DM_W_X[[2]], DM_W_X[[3]], DM_W_X_original[[1]], DM_W_X_original[[2]])) #DM_WPCA[[1]], DM_WPCA[[2]], DM_WPCA[[4]], DM_WPCA[[5]],
  # return [[1]]: Data_n_more
  # return [[2]]: fitting stat from Data_n_more: dispersion is [,1]
  # return [[3]]: fitted mu from Data_n_more
  # return [[4]]: Xij matrix
  # return [[5]]: Weight matrix
  # return [[6]]: UMI
} 

scRNA_cluster_fixTheta <- function(DM, t, n, I) {
  # DM: data to be analyzed
  # n: threshold for removing genes with less than n cells have reading
  # k: Number of principal components to retrieve.
  # t: Number of iterations (allows convergence to the dominant principal component)
  DM <- as.matrix(DM)
  DM_1 <- scRNA_cluster_glm(DM, n)
  DM_W_X <- scRNA_cluster_taylor(DM_1[[1]], log(DM_1[[3]]), as.numeric(DM_1[[2]][,1]), I)
  theta <- rep(t, nrow(DM_1[[1]]))
  v0 <- matrix(rep(-1, nrow(DM_1[[1]])*ncol(DM_1[[1]])), nrow = nrow(DM_1[[1]]), ncol = ncol(DM_1[[1]]))
  DM_W_X_fix <- scRNA_cluster_taylor(DM_1[[1]], v0, theta, I)
  DM_W_X_fixtheta <- scRNA_cluster_taylor(DM_1[[1]], log(DM_1[[3]]), theta, I)
  #DM_Center <- PC_center(DM_W_X[[1]])
  #DM_WPCA <- wPCA(DM_W_X[[1]], DM_W_X[[2]])
  return(list(DM_1[[1]], DM_1[[2]], DM_1[[3]], 
              DM_W_X[[1]], DM_W_X[[2]], DM_W_X[[3]], DM_W_X_fix[[1]], DM_W_X_fix[[2]], DM_W_X_fixtheta[[1]], DM_W_X_fixtheta[[2]])) #DM_WPCA[[1]], DM_WPCA[[2]], DM_WPCA[[4]], DM_WPCA[[5]],
  # return [[1]]: Data_n_more
  # return [[2]]: fitting stat from Data_n_more: dispersion is [,1]
  # return [[3]]: fitted mu from Data_n_more
  # return [[4]]: Xij matrix
  # return [[5]]: Weight matrix
  # return [[6]]: UMI
  # return [[7]]: Xij matrix with fixed v0 (-1) and theta
  # return [[8]]: Weight matrixwith fixed v0 (-1) and theta
} 

scRNA_cluster_interation <- function(svdResult, yij, theta, UMI, I) {
  U = as.matrix(svdResult$U)
  S = as.matrix(svdResult$S)
  V = as.matrix(svdResult$V)
  AZ = V%*%S%*%t(U)
  new_mu0 <- matrix(nrow = nrow(AZ), ncol = ncol(AZ), byrow = F)
  for (i in 1:dim(AZ)[2]) {
    new_mu0[,i] <- AZ[,i]+log(UMI[i])
  }
  D_WPCA_new <- scRNA_cluster_taylor(yij, new_mu0, theta, I)
  return(list(D_WPCA_new[[1]], D_WPCA_new[[2]]))
}

# Restrict theta to >10 (dispersion < 0.1)
scRNA_cluster_theta_correctAll <- function(DM, n, I) {
  # DM: data to be analyzed
  # n: threshold for removing genes with less than n cells have reading
  # k: Number of principal components to retrieve.
  # t: Number of iterations (allows convergence to the dominant principal component)
  DM <- as.matrix(DM)
  DM_1 <- scRNA_cluster_glm(DM, n)
  theta <- c()
  for (i in 1:nrow(DM_1[[2]])) {
    if(DM_1[[2]][i,1]<10) {
      theta[i] = 10
    } else if (DM_1[[2]][i,1]>=10) {
      theta[i] = DM_1[[2]][i,1]
    }
  }
  fit_stat <- data.frame(DM_1[[2]], theta)
  v0 <- matrix(rep(-2, nrow(DM_1[[3]])*ncol(DM_1[[3]])), nrow = nrow(DM_1[[3]]), ncol = ncol(DM_1[[3]]))
  DM_W_X <- scRNA_cluster_taylor(DM_1[[1]], v0, theta, I)
  #DM_Center <- PC_center(DM_W_X[[1]])
  #DM_WPCA <- wPCA(DM_W_X[[1]], DM_W_X[[2]])
  return(list(DM_1[[1]], fit_stat, DM_1[[3]], 
              DM_W_X[[1]], DM_W_X[[2]], DM_W_X[[3]]))
  # return [[1]]: Data_n_more
  # return [[2]]: fitting stat from Data_n_more: dispersion is [,c(1,7)]
  # return [[3]]: fitted mu from Data_n_more
  # return [[4]]: Xij matrix
  # return [[5]]: Weight matrix
  # return [[6]]: UMI
} 

scRNA_cluster_theta_correct_OnlySwap <- function(DM, n, I, gene_index) {
  # DM: data to be analyzed
  # n: threshold for removing genes with less than n cells have reading
  # k: Number of principal components to retrieve.
  # t: Number of iterations (allows convergence to the dominant principal component)
  DM <- as.matrix(DM)
  DM_1 <- scRNA_cluster_glm(DM, n)
  theta_change <- DM_1[[2]][gene_index, 1]
  theta_tobe <- c()
  for (i in 1:length(theta_change)) {
    if(theta_change[i]<10) {
      theta_tobe[i] = 10
    } else if (theta_change[i]>=10) {
      theta_tobe[i] = theta_change[i]
    }
  }
  theta <- DM_1[[2]][,1]
  theta[gene_index] <- theta_tobe
  theta <- as.vector(theta)
  v0 <- matrix(rep(-2, nrow(DM_1[[3]])*ncol(DM_1[[3]])), nrow = nrow(DM_1[[3]]), ncol = ncol(DM_1[[3]]))
  DM_W_X <- scRNA_cluster_taylor(DM_1[[1]], v0, theta, I)
  #DM_Center <- PC_center(DM_W_X[[1]])
  #DM_WPCA <- wPCA(DM_W_X[[1]], DM_W_X[[2]])
  return(list(DM_1[[1]], DM_1[[2]], DM_1[[3]], 
              DM_W_X[[1]], DM_W_X[[2]], DM_W_X[[3]]))
  # return [[1]]: Data_n_more
  # return [[2]]: fitting stat from Data_n_more: dispersion is [,1]
  # return [[3]]: fitted mu from Data_n_more
  # return [[4]]: Xij matrix
  # return [[5]]: Weight matrix
  # return [[6]]: UMI
} 

# Function to divide each gene with its cell's total UMI
UMI_frac <- function(DM) {
  UMI <- colSums(DM)
  DM <- as.matrix(DM)
  DM_new <- matrix(nrow = nrow(DM), ncol = ncol(DM), byrow = TRUE)
  for (i in 1:nrow(DM)) {
    DM_new[i,] <- DM[i,]/UMI
  }
  return(DM_new)
}

# Wrapper function to run the algarithm using PCA as initial v_0_ij
scRNA_cluster_initByPCA <- function(DM,v0, n, I) {
  # DM: data to be analyzed
  # n: threshold for removing genes with less than n cells have reading (use 0 for not removing any gene)
  # v0: initiate v_0_ij from PCA
  # I: Even UMI (I=0) or not (I=1)
  DM <- as.matrix(DM)
  DM_1 <- scRNA_cluster_glm(DM, n)
  DM_W_X <- scRNA_cluster_taylor(DM, v0, as.numeric(DM_1[[2]][,1]), I)
  #DM_Center <- PC_center(DM_W_X[[1]])
  #DM_WPCA <- wPCA(DM_W_X[[1]], DM_W_X[[2]])
  return(list(DM_1[[1]], DM_1[[2]], DM_1[[3]], 
              DM_W_X[[1]], DM_W_X[[2]], DM_W_X[[3]])) #DM_WPCA[[1]], DM_WPCA[[2]], DM_WPCA[[4]], DM_WPCA[[5]],
  # return [[1]]: Data_n_more
  # return [[2]]: fitting stat from Data_n_more: dispersion is [,1]
  # return [[3]]: fitted mu from Data_n_more
  # return [[4]]: Xij matrix
  # return [[5]]: Weight matrix
  # return [[6]]: UMI
} 

  






