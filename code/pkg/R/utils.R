rdirichlet <- function(n,alpha){
  # Random Dirichlet generator
  k <- length(alpha)
  x <- lapply(seq_len(n),function(x){
    xi <- rgamma(k,shape = alpha,scale = 1)
    xi/sum(xi)
  })
  do.call(cbind,x)
}

est.moments.dirichlet <- function(p){
  # Estimates parameters of a Dirichlet distribution via method of moments
  e <- rowMeans(p)
  e2 <- rowMeans(p^2)
  k <- which.max(e)
  a0 <- (e[k] - e2[k])/(e2[k] - e[k]^2)
  e*a0
}

g <- function(n,alpha,p,epsilon = .Machine$double.xmin){
  # Gradient of a Dirichlet distribution
  n*(digamma(sum(alpha)) - digamma(alpha) + rowMeans(log(p + epsilon)))
}

invHg <- function(n,alpha,p){
  # Inverse of Hessian*Gradient of a Dirichlet distribution
  g.obs <- g(n,alpha,p)
  q.obs <- -n*trigamma(alpha)
  z.obs <- n*trigamma(sum(alpha))
  b.obs <- sum(g.obs/q.obs)/(1/z.obs + sum(1/q.obs))
  (g.obs - b.obs)/q.obs
}

est.newton.dirichlet <- function(p){
  # Estimates parameters of a Dirichlet distribution via Newton method
  n <- ncol(p)
  alpha.init <- est.moments.dirichlet(p)
  
  error <- 1
  alpha.old <- alpha.init
  while(error > 1e-6){
    alpha.new <- alpha.old - invHg(n,alpha.old,p)
    error <- max(abs((alpha.new - alpha.old)/alpha.old))
    alpha.old <- alpha.new
  }
  if(any(alpha.new < 0)){
    return(alpha.init)
  }
  alpha.new
}

catchSalmon2 <- function (paths, verbose = TRUE) 
{
  # This function is analogous to edgeR's catchSalmon but brings in the TPM
  # output from Salmon, in addition to counts.
  
  NSamples <- length(paths)
  OK <- requireNamespace("jsonlite", quietly = TRUE)
  if (!OK) 
    stop("jsonlite package required but is not installed (or can't be loaded)")
  OK <- requireNamespace("readr", quietly = TRUE)
  if (!OK) 
    stop("readr package required but is not installed (or can't be loaded)")
  for (j in 1L:NSamples) {
    if (verbose) 
      cat("Reading ", paths[j], ", ", sep = "")
    MetaFile <- file.path(paths[j], "aux_info", "meta_info.json")
    QuantFile <- file.path(paths[j], "quant.sf")
    BootFile <- file.path(paths[j], "aux_info", "bootstrap", 
                          "bootstraps.gz")
    if (!file.exists(QuantFile)) 
      stop("quant.sf file not found at specified path")
    Meta <- jsonlite::fromJSON(MetaFile)
    NTx <- Meta$num_targets
    if (is.null(NTx)) 
      NTx <- Meta$num_valid_targets
    if (is.null(NTx)) 
      stop("Can't find number of targets")
    NBoot <- Meta$num_bootstraps
    if (is.null(NBoot)) 
      stop("Can't find number of bootstraps")
    if (verbose) 
      cat(NTx, "transcripts,", NBoot, "bootstraps\n")
    if (j == 1L) {
      Counts <- matrix(0, NTx, NSamples)
      TPM <- matrix(0, NTx, NSamples)
      DF <- rep_len(0L, NTx)
      OverDisp <- rep_len(0, NTx)
      Quant1 <- suppressWarnings(readr::read_tsv(QuantFile, 
                                                 col_types = "cdddd", progress = FALSE))
      Counts[, 1L] <- Quant1$NumReads
      TPM[, 1L] <- Quant1$TPM
    } else {
      Quant <- suppressWarnings(readr::read_tsv(QuantFile, 
                                                col_types = "___dd", progress = FALSE))
      Counts[, j] <- Quant$NumReads
      TPM[, j] <- Quant$TPM
    }
    if (NBoot > 0L) {
      BootFileCon <- gzcon(file(BootFile, open = "rb"))
      Boot <- readBin(BootFileCon, what = "double", n = NTx * 
                        NBoot)
      close(BootFileCon)
      dim(Boot) <- c(NTx, NBoot)
      M <- rowMeans(Boot)
      i <- (M > 0)
      OverDisp[i] <- OverDisp[i] + rowSums((Boot[i, ] - 
                                              M[i])^2)/M[i]
      DF[i] <- DF[i] + NBoot - 1L
    }
  }
  i <- (DF > 0L)
  if (sum(i) > 0L) {
    OverDisp[i] <- OverDisp[i]/DF[i]
    DFMedian <- median(DF[i])
    DFPrior <- 3
    OverDispPrior <- median(OverDisp[i])/qf(0.5, df1 = DFMedian, 
                                            df2 = DFPrior)
    if (OverDispPrior < 1) {
      OverDispPrior <- 1
    }
    OverDisp[i] <- (DFPrior * OverDispPrior + DF[i] * OverDisp[i])/(DFPrior + DF[i])
    OverDisp <- pmax(OverDisp, 1)
    OverDisp[!i] <- OverDispPrior
  } else {
    OverDisp[] <- NA_real_
    OverDispPrior <- NA_real_
  }
  Quant1 <- as.data.frame(Quant1, stringsAsFactors = FALSE)
  dimnames(Counts) <- list(Quant1$Name, paths)
  dimnames(TPM) <- list(Quant1$Name, paths)
  row.names(Quant1) <- Quant1$Name
  Quant1$Name <- NULL
  Quant1$TPM <- Quant1$NumReads <- NULL
  Quant1$Overdispersion <- OverDisp
  list(counts = Counts, tpm = TPM, annotation = Quant1, overdispersion.prior = OverDispPrior)
}

catchKallisto2 <- function(paths,verbose = TRUE){
  NSamples <- length(paths)
  suppressPackageStartupMessages(OK <- requireNamespace("rhdf5", quietly = TRUE))
  if (!OK) 
    stop("rhdf5 package required but is not installed (or can't be loaded)")
  for (j in 1L:NSamples) {
    if (verbose) 
      cat("Reading ", paths[j], ", ", sep = "")
    h5File <- file.path(paths[j], "abundance.h5")
    tsvFile <- file.path(paths[j],'abundance.tsv')
    if (!file.exists(h5File)) 
      stop("abundance.h5 file not found at specified path")
    h5 <- rhdf5::H5Fopen(h5File)
    Quant <- suppressWarnings(readr::read_tsv(tsvFile,col_types = "____d", progress = FALSE))
    aux <- h5$aux
    NTx <- length(aux$lengths)
    NBoot <- as.integer(aux$num_bootstrap)
    if (verbose) 
      cat(NTx, "transcripts,", NBoot, "bootstraps\n")
    if (j == 1L) {
      Counts <- matrix(0, NTx, NSamples)
      TPM <- matrix(0, NTx, NSamples)
      DF <- rep_len(0L, NTx)
      OverDisp <- rep_len(0, NTx)
    }
    Counts[, j] <- h5$est_counts
    TPM[, j] <- Quant$tpm
    if (NBoot > 0L) 
      Boot <- do.call(cbind, h5$bootstrap)
    rhdf5::H5Fclose(h5)
    if (NBoot > 0L) {
      M <- rowMeans(Boot)
      i <- (M > 0)
      OverDisp[i] <- OverDisp[i] + rowSums((Boot[i, ] - M[i])^2)/M[i]
      DF[i] <- DF[i] + NBoot - 1L
    }
  }
  i <- (DF > 0L)
  if (sum(i) > 0L) {
    OverDisp[i] <- OverDisp[i]/DF[i]
    DFMedian <- median(DF[i])
    DFPrior <- 3
    OverDispPrior <- median(OverDisp[i])/qf(0.5, df1 = DFMedian, df2 = DFPrior)
    if (OverDispPrior < 1) {
      OverDispPrior <- 1
    }
    OverDisp[i] <- (DFPrior * OverDispPrior + DF[i] * OverDisp[i])/(DFPrior + DF[i])
    OverDisp <- pmax(OverDisp, 1)
    OverDisp[!i] <- OverDispPrior
  } else {
    OverDisp[] <- NA_real_
    OverDispPrior <- NA_real_
  }
  Ann <- data.frame(Length = as.integer(aux$lengths), EffectiveLength = aux$eff_lengths, 
                    Overdispersion = OverDisp, row.names = aux$ids, stringsAsFactors = FALSE)
  dimnames(Counts) <- list(aux$ids, paths)
  dimnames(TPM) <- list(aux$ids, paths)
  list(counts = Counts, tpm = TPM, annotation = Ann, overdispersion.prior = OverDispPrior)
}