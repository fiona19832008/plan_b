# Generate two cluster
Two_Cluster_new <- function(DM, cell, seed, p_gene, q_high, g_fold) {
  # DM: data matrix
  # cell: number of cells to be sampled
  # seed: seed
  # p_gene: percent gene to be swapped
  # q_high: quantile of high expressed genes
  # g_fold: fold change
  
  set.seed(seed)
  # sample a number of cells
  R <- sample(ncol(DM), cell)
  DM_S <- DM[, R]
  DM_S <- Remove0counts(DM_S, 4)
  Gene_order <- order(rowMeans(DM_S))
  Gene_exp <- sort(rowMeans(DM_S),decreasing = F)
  Gene_name <- names(Gene_exp)
  names(Gene_exp) <- NULL
  Gene_in_order <- data.frame(Gene_order, Gene_name, Gene_exp)
  # split data to two part
  R3 <- sample(cell, round(cell/2))
  DM_S_1 <- DM_S[, R3]
  DM_S_2 <- DM_S[, -R3]
  rownames(DM_S_2) <- NULL
  # Calculate the start and end of gene to be swapped
  No_gene <- round(nrow(Gene_in_order)*p_gene)
  H_exp <- round(nrow(Gene_in_order)*q_high)
  L_exp <- sum(Gene_in_order[,3]<=Gene_in_order[H_exp,3]/g_fold)
  p <- L_exp-No_gene-1
  q <- L_exp
  s <- H_exp-No_gene-1
  t <- H_exp
  Gene_used <- Gene_in_order[c(p, q, s, t),]
  rownames(Gene_used) <- c(p,q,s,t)
  Low_Gene <- DM_S_2[Gene_order[p:q], ]
  High_Gene <- DM_S_2[Gene_order[s:t], ]
  DM_S_2new <- DM_S_2
  DM_S_2new[Gene_order[p:q],] <- High_Gene
  DM_S_2new[Gene_order[s:t],] <- Low_Gene
  rownames(DM_S_2new) <- rownames(DM_S_1)
  DM_two <- data.frame(DM_S_1, DM_S_2new)
  ID <- c(rep(1,ncol(DM_S_1)), rep(2, ncol(DM_S_2)))
  gene_index <- Gene_order[c(p:q, s:t)]
  return(list(DM_two, ID, Gene_used, Gene_in_order, gene_index))
}

Two_Cluster_DropNoChangeGene <- function(DM, cell, seed, p_gene, q_high, g_fold, Gene_drop, Low_p, Up_p) {
  # DM: data matrix
  # cell: number of cells to be sampled
  # seed: seed
  # p_gene: percent gene to be swapped
  # q_high: quantile of high expressed genes
  # g_fold: fold change
  
  set.seed(seed)
  # sample a number of cells
  R <- sample(ncol(DM), cell)
  DM_S <- DM[, R]
  DM_S <- Remove0counts(DM_S, 4)
  Gene_order <- order(rowMeans(DM_S))
  Gene_exp <- sort(rowMeans(DM_S),decreasing = F)
  Gene_name <- names(Gene_exp)
  names(Gene_exp) <- NULL
  Gene_in_order <- data.frame(Gene_order, Gene_name, Gene_exp)
  # split data to two part
  R3 <- sample(cell, round(cell/2))
  DM_S_1 <- DM_S[, R3]
  DM_S_2 <- DM_S[, -R3]
  rownames(DM_S_2) <- NULL
  # Calculate the start and end of gene to be swapped
  No_gene <- round(nrow(Gene_in_order)*p_gene)
  H_exp <- round(nrow(Gene_in_order)*q_high)
  L_exp <- sum(Gene_in_order[,3]<=Gene_in_order[H_exp,3]/g_fold)
  p <- L_exp-No_gene-1
  q <- L_exp
  s <- H_exp-No_gene-1
  t <- H_exp
  Gene_used <- Gene_in_order[c(p, q, s, t),]
  rownames(Gene_used) <- c(p,q,s,t)
  Low_Gene <- DM_S_2[Gene_order[p:q], ]
  High_Gene <- DM_S_2[Gene_order[s:t], ]
  set.seed(seed)
  if (Gene_drop!=0) {
    Other <- sample(Gene_order[-c(p:q, s:t)], round(length(Gene_order[-c(p:q, s:t)])*Gene_drop))
  } else if (Gene_drop==0) {
    Other <- 0
  }
  if (Up_p != 0) {
    Upper <- Gene_order[(length(Gene_order)-round((length(Gene_order)-t)*Up_p)):length(Gene_order)]
  } else if (Up_p == 0) {
    Upper <- 0
  }
  if (Low_p!=0) {
    Low <- Gene_order[1:(round(p*Low_p)-1)]
  } else if (Low_p == 0) {
    Low <- 0
  }
  DM_S_2new <- DM_S_2
  DM_S_2new[Gene_order[p:q],] <- High_Gene
  DM_S_2new[Gene_order[s:t],] <- Low_Gene
  rownames(DM_S_2new) <- rownames(DM_S_1)
  DM_two <- data.frame(DM_S_1, DM_S_2new)
  if (Other == 0 & Up_p ==0) {
    DM_two_1 <- DM_two[-Low,]
  } else if (Other == 0 & Low_p ==0) {
    DM_two_1 <- DM_two[-Upper,]
  } else if (Low_p == 0 & Up_p ==0) {
    DM_two_1 <- DM_two[-Other,]
  } else if (Low_p != 0 & Up_p !=0) {
    DM_two_1 <- DM_two[-c(Low, Upper)]
  }
  ID <- c(rep(1,ncol(DM_S_1)), rep(2, ncol(DM_S_2)))
  return(list(DM_two_1, DM_two, ID, Gene_used, Gene_in_order, Other, Upper, Low))
}


# Random sample two data set and then downsample one with 70% total original UMI
UMI_dif <- function(DM, seed1, seed2, n) {
  set.seed(seed1)
  S1 <- sample(ncol(DM), n)
  DM_sample <- DM[,S1]
  Clean_DM_sample <- Remove0counts(DM_sample)
  set.seed(seed2)
  S2 <- sample(n, n/2)
  Data_a <- Clean_DM_sample[,S2]
  Data_b <- Clean_DM_sample[,-S2]
  print(summary(colSums(Data_a)))
  print(summary(colSums(Data_b)))
  newTotalCounts <- colSums(Data_b)
  Data_b_down = matrix(NA, dim(Data_b)[1], dim(Data_b)[2])
  for (i in 1:dim(Data_b)[2]) {
    Data_b_down[, i] = downSampleVector(Data_b[, i], 
                                        newTotalCounts = round(newTotalCounts[i] * runif(1, 0.1, 0.5)))
  }
  rownames(Data_b_down) <- rownames(Data_b)
  colnames(Data_b_down) <- colnames(Data_b)
  return(list(Data_a, Data_b_down, Data_b))
}

# function for remove gene with < n cells with counts 
Remove0counts <- function(DataMatrix, n) {
  no0 <- DataMatrix>0
  gsum <- c()
  for (i in 1:nrow(DataMatrix)) {
    gsum[i] <- sum(no0[i,])
  }
  gsumN <- which(gsum > n)
  CleanedData <- DataMatrix[gsumN,]
  return(CleanedData)
}

RemoveNcounts <- function(DataMatrix, n, m) {
  no0 <- DataMatrix>n
  gsum <- c()
  for (i in 1:nrow(DataMatrix)) {
    gsum[i] <- sum(no0[i,])
  }
  gsumN <- which(gsum > m)
  CleanedData <- DataMatrix[gsumN,]
  return(CleanedData)
}

