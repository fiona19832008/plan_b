main = function() {
  source("scRNA_cluster.R")
  source("runWSVD.R")
  source("A1A2loaddata.R")
  source("Simulate_data_new.R")
  #source("Downsampling function.R")
  # Run the script as follows:
  # Rscript <No_cell> <DataPath> <WSVDpath> <PCApath> <p_gene> <q_high> <g_fold> <k>
  
  args = commandArgs(T)
  No_cell = args[1]
  DataPath = args[2]
  WSVDpath = args[3]
  PCApath = args[4]
  p_gene = args[5]
  q_high = args[6]
  g_fold = args[7]
  k = args[8] # data set name
  I = args[9]
  
  # 4000, 4100, 6000, 6100
  DM <- loadData(paste0(DataPath, ".txt"), 
                    paste0(DataPath, ".group"))
  
  DM_swap <- Two_Cluster_new(DM[[1]], as.numeric(No_cell), 1234, 
                         as.numeric(p_gene), as.numeric(q_high), 
                         as.numeric(g_fold))
  
  rm(DM)
  # browser()
  data = t(as.matrix(DM_swap[[1]]))
  N1 = sum(DM_swap[[2]]==1)
  N2 = sum(DM_swap[[2]]==2)
  write.table(DM_swap[[3]], paste0(PCApath,k,"_",as.numeric(p_gene)*100,"_",as.numeric(q_high)*100,"_",as.numeric(g_fold)*10,"_swap.txt"), sep=" ", quote = F)
  # TPM
  data_tpm = (data*1e6)/rowSums(data)
  require(stats)
  D_PCA <- prcomp(log(data_tpm + 0.1))
  write.table(D_PCA$x, paste0(PCApath,k,"_",as.numeric(p_gene)*100,"_",as.numeric(q_high)*100,"_",as.numeric(g_fold)*10,"_PC.txt"), sep=" ", quote = F)
  write.table(D_PCA$rotation, paste0(PCApath,k,"_",as.numeric(p_gene)*100,"_",as.numeric(q_high)*100,"_",as.numeric(g_fold)*10, "_loading.txt"), sep=" ", quote = F)
  pdf(paste0(PCApath,k,"_",as.numeric(p_gene)*100,"_",as.numeric(q_high)*100,"_",as.numeric(g_fold)*10,"_PCA.pdf"))
  plot(D_PCA$x[,1:2], col=c(rep("red",N1), rep("blue",N2)))
  plot(D_PCA$x[,2:3], col=c(rep("red",N1), rep("blue",N2)))
  plot(D_PCA$x[,4:5], col=c(rep("red",N1), rep("blue",N2)))
  plot(D_PCA$x[,6:7], col=c(rep("red",N1), rep("blue",N2)))
  dev.off()
  
  D_WPCA <- scRNA_cluster(t(data), 5, as.numeric(I))
  X = D_WPCA[[4]]
  W = D_WPCA[[5]]
  
  #browser()
  svdResult = runWeightedSVD(X, W, 8)
  write.table(svdResult$U, paste0(WSVDpath,k,"_",as.numeric(p_gene)*100,"_",as.numeric(q_high)*100,"_",as.numeric(g_fold)*10, "_U.txt"), sep=" ", quote = F)
  write.table(svdResult$S, paste0(WSVDpath,k,"_",as.numeric(p_gene)*100,"_",as.numeric(q_high)*100,"_",as.numeric(g_fold)*10, "_S.txt"), sep=" ", quote = F)
  write.table(svdResult$V, paste0(WSVDpath,k,"_",as.numeric(p_gene)*100,"_",as.numeric(q_high)*100,"_",as.numeric(g_fold)*10, "_V.txt"), sep=" ", quote = F)
  pdf(paste0(WSVDpath,k,"_",as.numeric(p_gene)*100,"_",as.numeric(q_high)*100,"_",as.numeric(g_fold)*10, "_WSVD.pdf"))
  par(mfrow = c(1,1), pty = "s")
  plot(svdResult$U[, 1:2], col=c(rep("red",N1), rep("blue",N2)))
  plot(svdResult$U[, 2:3], col=c(rep("red",N1), rep("blue",N2)))
  plot(svdResult$U[, 4:5], col=c(rep("red",N1), rep("blue",N2)))
  plot(svdResult$U[, 5:6], col=c(rep("red",N1), rep("blue",N2)))
  plot(svdResult$U[, 7:8], col=c(rep("red",N1), rep("blue",N2)))
  dev.off()
}

main()
