runWeightedSVD = function(X, W, k) {
  # call matlab function weightedSVD.m 
  # assume in the directory
  
  # input
  # X: X matrix
  # W: W matrix
  
  # output
  # SVD results
  # a list of U, S, V, where U is like principal components
  
  # create a temporary folder
  temporaryFolder = paste0(Sys.info()["nodename"], "-", Sys.getpid())
  dir.create(temporaryFolder)
  XFile = paste0(temporaryFolder, "/XFile.txt")
  WFile = paste0(temporaryFolder, "/WFile.txt")
  UFile = paste0(temporaryFolder, "/UFile.txt")
  SFile = paste0(temporaryFolder, "/SFile.txt")
  VFile = paste0(temporaryFolder, "/VFile.txt")
  cmdFile = paste0(temporaryFolder, "/cmd.sh")
 
  write.table(X, XFile, quote = F, col.names = F, row.names = F)
  write.table(W, WFile, quote = F, col.names = F, row.names = F)   

  file.create(cmdFile)
  cat("module load matlab\n", file = cmdFile, append = T)
  cat("matlab -nodesktop -nosplash -r \"", file = cmdFile, append = T)
  cat("XFile = \'", XFile, "\'; WFile = \'", WFile, "\'; UFile = \'", UFile, 
      "\'; SFile = \'", SFile, "\'; VFile = \'", VFile, "\'; ",  
     "k = ", k, "; ",   
     " weightedSVD; exit;\"", file = cmdFile, append = T, sep = "")
  
  system(paste("bash", cmdFile))
  
  # read the result
  U = read.table(UFile, as.is = T, header = F)
  S = read.table(SFile, as.is = T, header = F)
  V = read.table(VFile, as.is = T, header = F)
  
  unlink(temporaryFolder)
  
  return(list = list(U = U, S = S, V = V))
}
