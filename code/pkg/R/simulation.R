simulateExpr <- function(n.feat,n.libs,lib.sizes,top.cutoff,num.DE,fc,
                         df.bcv = 40, bcv.true = 0.2){
  # This function was written with the goal of mimicking the simulation setup
  # used in the voom paper. Specifically, we use goodTuringProportions to
  # estimate baseline abundances, different asymptotic BCV value as well as
  # Chi-square degrees of freedom, and a different dispersion trend.
  # See ?baselineAbundance_function for goodTuringProportions usage in this
  # simulation.
  
  n.groups <- length(n.libs)
  n.samples <- sum(n.libs)
  if (n.groups > 2) stop('This function does not support more than 2 groups')
  
  # Generating baseline proportions
  goodTuring <- get('baselineAbundance_function')
  baselineAbundance <- goodTuring(seq_len(n.feat) / (n.feat + 1))
  baselineAbundance <- baselineAbundance/sum(baselineAbundance)
  baselineAbundancePerGroup <- matrix(baselineAbundance,ncol = n.groups,nrow = n.feat)
  
  # Generating DE status for 2 groups
  topTx <- order(baselineAbundance,decreasing = TRUE)[seq(from = 1,length.out = min(top.cutoff,n.feat))]
  topTxDE <- sample(topTx,num.DE,replace = FALSE)
  names(topTxDE) <- as.character(sample(1:4,length(topTxDE),replace = TRUE))
  
  de.Up1.Base2 <- topTxDE[names(topTxDE) == '1']
  de.Down1.Base2 <- topTxDE[names(topTxDE) == '2']
  de.Up2.Base1 <- topTxDE[names(topTxDE) == '3']
  de.Down2.Base1 <- topTxDE[names(topTxDE) == '4']
  
  baselineAbundancePerGroup[de.Up1.Base2,1] <- fc*baselineAbundancePerGroup[de.Up1.Base2,1]
  baselineAbundancePerGroup[de.Down1.Base2,1] <- (1/fc)*baselineAbundancePerGroup[de.Down1.Base2,1]
  baselineAbundancePerGroup[de.Up2.Base1,2] <- fc*baselineAbundancePerGroup[de.Up2.Base1,2]
  baselineAbundancePerGroup[de.Down2.Base1,2] <- (1/fc)*baselineAbundancePerGroup[de.Down2.Base1,2]
  
  de <- vector('numeric',n.feat)
  de[topTxDE[names(topTxDE) %in% c('2','3')]] <- 1L # Up in group 2
  de[topTxDE[names(topTxDE) %in% c('1','4')]] <- -1L # Down in group 2
  
  # Generating expected counts
  mu0 <- lapply(seq_len(n.groups),function(x){
    size <- lib.sizes[seq(1 + n.libs[x] * (x - 1), n.libs[x] * x)]
    matrix(baselineAbundancePerGroup[,x],n.feat,1) %*% matrix(size,1,length(size))
  })
  mu0 <- do.call(cbind,mu0)
  
  # Generating random noise around dispersion trend. I am generating 1 RV per
  # transcript per group. In the voom simulation, we had 1 RV per gene per sample.
  # Here, I am arguing that the expression level and dispersion should be exactly
  # the same among libraries of the same group. Since we already generated DE
  # status in the steps above, it makes sense to have 1 random shift around the
  # dispersion trend per group and a fixed resulting dispersion per group.
  chisq <- lapply(seq_len(n.groups),function(x){
    rv <- df.bcv / rchisq(n.feat, df = df.bcv)
    matrix(rv,ncol = n.libs[x],nrow = n.feat)
  })
  chisq <- do.call(cbind,chisq)
  
  # Biological variation and Dispersion trend
  bcv0 <- bcv.true + 1/sqrt(mu0)
  disp <- bcv0 ^ 2 * chisq
  
  # Biological variation
  shape <- 1/disp
  scale <- mu0/shape
  
  expr <- matrix(rgamma(n.feat * n.samples, shape = shape, scale = scale),
                 nrow =  n.feat, ncol =  n.samples)
  
  # Shuffling matrices (to be matched with reference data)
  key <- sample(n.feat,replace = FALSE)
  expr <- expr[key,]
  mu0 <- mu0[key,]
  disp <- disp[key,]
  de <- de[key]
  
  return(list('expr' = expr,'mu' = mu0,'disp' = disp,'de' = de))
}

#' @importFrom data.table setkey copy setnames
simulateTPM <- function(contigs,contigs.subset,
                        n.libs,lib.sizes,top.cutoff,num.DE,fc){
  
  # Generating sample labels
  group <- rep(LETTERS[seq_len(length(n.libs))],times = n.libs)
  rep <- unlist(lapply(n.libs,seq_len))
  group.name <- paste(paste0('group', group), paste0('rep',rep), sep = '_')
  
  # Simulating transcript-wise expression
  trExpr <- simulateExpr(n.feat = nrow(contigs.subset),
                         n.libs = n.libs,lib.sizes = lib.sizes,
                         top.cutoff = top.cutoff,num.DE = num.DE,fc = fc)
  
  # Generating TPM values
  tpm <- trExpr$expr / contigs.subset$Length
  tpm <- 1e6 * t(t(tpm) / colSums(tpm))
  
  # Organizing contigs.subset
  dimnames(tpm) <- dimnames(trExpr$expr) <- dimnames(trExpr$mu) <- dimnames(trExpr$disp) <- list(contigs.subset$TranscriptID,group.name)
  names(trExpr$de) <- contigs.subset$TranscriptID
  
  # Merging contigs.subset values to contigs
  tpm.contigs <- expr.contigs <- mu.contigs <- disp.contigs <- matrix(NA, nrow(contigs), sum(n.libs))
  de.contigs <- rep(NA,nrow(contigs))
  
  dimnames(tpm.contigs) <- dimnames(expr.contigs) <- dimnames(mu.contigs) <- dimnames(disp.contigs) <- list(contigs$TranscriptID,group.name)
  names(de.contigs) <- contigs$TranscriptID
  
  
  key <- match(contigs.subset$TranscriptID, contigs$TranscriptID)
  tpm.contigs[key,] <- tpm
  expr.contigs[key,] <- trExpr$expr
  mu.contigs[key,] <- trExpr$mu
  disp.contigs[key,] <- trExpr$disp
  de.contigs[key] <- trExpr$de
  
  return(list('tpm' = tpm.contigs,'expr' = expr.contigs,'mu' = mu.contigs,'disp' = disp.contigs,'de' = de.contigs))
}

#' @importFrom data.table as.data.table
selectTx <- function(contigs,genome,max.tx){
  DT <- as.data.table(contigs)
  
  # Selecting reference ranking
  if (is.character(genome)) {
    if (genome == 'mm39') DT.ref <- get('GSE60450')
    if (genome == 'hg38') stop('simulation for human-like data not yet done')
  } else {
    DT.ref <- genome
  }
  
  # Merging rankings and keeping reference transcripts
  key <- match(DT$TranscriptID,DT.ref$TranscriptID)
  DT$TxRank <- DT.ref$Rank[key]
  DT <- DT[!is.na(DT$TxRank),]
  
  # Subsetting transcripts
  filtrTx <- seq(from = 1,length.out = min(max.tx,max(DT$TxRank)))
  DT.sub <- DT[TxRank %in% filtrTx,selectTx := TRUE,by = 'GeneID'][(selectTx),]
  DT.sub$selectTx <- NULL
  as.data.frame(DT.sub)
}

#' @importFrom limma strsplit2
readFasta <- function(fasta){
  # Reading input and getting tx information from fasta
  contigs <- scanFasta(fasta,quiet = TRUE)
  tx.info <- strsplit2(contigs$TranscriptID,"\\|")
  contigs[,c('TranscriptID','GeneID','TranscriptType')] <- tx.info[,c(1,2,8)]
  contigs
}

#' @importFrom Rsubread scanFasta simReads
#' @importFrom BiocParallel bplapply
#' @importFrom readr write_tsv
simulateFASTQ <- function(fasta,n.libs,lib.sizes,dest,tmpdir,paired.end,
                          fc,num.DE,top.cutoff,max.tx,genome,BPPARAM,read.length,
                          fragment.length.min){
  
  # Checking if simulation has already been run
  dir.create(dest,showWarnings = FALSE,recursive = TRUE)
  out.files <- file.path(dest,c('counts.tsv.gz','targets.tsv.gz'))
  names(out.files) <- c('counts','targets')
  if (all(file.exists(out.files))) {
    message('Reads already simulated!')
    return(invisible(NULL))
  }
  
  # Scanning fasta
  message('Preparing fasta file for simulation...')
  contigs <- readFasta(fasta)
  
  # Subsetting transcripts according to curated data and bringing metadata
  message('Subsetting transcripts...')
  contigs.subset <- selectTx(contigs,genome,max.tx)
  
  # Simulating tx-wise TPMs from genes
  message('Simulating transcript-wise TPM...')
  txTPM <- simulateTPM(contigs = contigs, contigs.subset = contigs.subset,
                       lib.sizes = lib.sizes,n.libs = n.libs, 
                       top.cutoff = top.cutoff,num.DE = num.DE,fc = fc)
  
  # Getting quality reference
  quality.source <- ifelse(read.length %in% c(75, 100),'Rsubread','rfun')
  quality.reference <- list.files(system.file(package = quality.source,'qualf'),
                                  paste0('-',read.length,'bp'),full.names = TRUE)
  
  # Running simReads
  message('Simulating reads...')
  curwd <- getwd()
  setwd(tmpdir)
  out <- bplapply(seq_len(sum(n.libs)), FUN = function(i){
    simReads(transcript.file = fasta, expression.levels = txTPM$tpm[, i], 
             output.prefix = colnames(txTPM$tpm)[i],library.size = lib.sizes[i],
             paired.end = paired.end, simplify.transcript.names = TRUE,
             fragment.length.min = fragment.length.min,read.length = read.length,
             quality.reference = quality.reference)
  },BPPARAM = BPPARAM)
  out.fastq <- normalizePath(list.files('.',"group*.*.fastq.gz",full.names = TRUE))
  setwd(curwd)
  
  # Saving metadata
  message('Saving metadata...')
  mat <- do.call(cbind,lapply(out,function(x){x$NReads}))
  colnames(mat) <- colnames(txTPM$tpm)
  
  out <- cbind(contigs[,c("TranscriptID","Length",'GeneID')],'status' = txTPM$de,mat)
  rownames(out) <- NULL
  write_tsv(x = out,file = out.files['counts'],col_names = TRUE,quote = 'none')
  
  targets <- data.frame('R1' = out.fastq[grepl("R1.fastq.gz",out.fastq)])
  if (isTRUE(paired.end)) {
    targets$R2 <- out.fastq[grepl("R2.fastq.gz",out.fastq)]
  }
  write_tsv(x = targets,file = out.files['targets'],col_names = TRUE,quote = 'none')
  
  message('Simulation of sequencing reads completed!')
}

#' @importFrom BiocParallel MulticoreParam register
simulateExperiment <- function(dest,
                               fasta,
                               genome,
                               tmpdir = tempdir(),
                               max.tx = Inf,
                               workers = 1,
                               num.DE = 3000,
                               fc = 2,
                               top.cutoff = 30000,
                               n.libs = c(3,3),
                               lib.sizes = rep(50e6,sum(n.libs)),
                               paired.end = FALSE,
                               keep.fastq = FALSE,
                               read.length = 75,
                               fragment.length.min = 150,
                               bin.salmon = "/stornext/General/data/academic/lab_smyth/baldoni.p/software/salmon-1.9.0_linux_x86_64/bin/salmon",
                               index.salmon = file.path("/stornext/General/data/academic/lab_smyth/baldoni.p/software/SalmonIndex",genome,"transcripts_index/"),
                               opts.salmon =  paste('-p',workers,'-l A --numBootstraps 100 --validateMappings'),
                               bin.kallisto = "/stornext/General/data/academic/lab_smyth/baldoni.p/software/kallisto/kallisto",
                               index.kallisto = file.path("/stornext/General/data/academic/lab_smyth/baldoni.p/software/kallistoIndex",genome,"transcripts_index"),
                               opts.kallisto = paste0('--bootstrap-samples=100 --threads=',workers),
                               run.salmon = TRUE,
                               run.kallisto = TRUE){
  
  # Setting up parallel computing
  BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE)
  register(BPPARAM = BPPARAM)
  
  # Preparing directories
  dir.create(tmpdir,recursive = TRUE,showWarnings = FALSE)
  dir.create(dest,recursive = TRUE,showWarnings = FALSE)
  
  tmpdir <- normalizePath(tmpdir)
  dest <- normalizePath(dest)
  
  # Simulating FASTQs
  simulateFASTQ(fasta = fasta,n.libs = n.libs,tmpdir = tmpdir,
                lib.sizes = lib.sizes,dest = file.path(dest,'meta'), fc = fc,
                paired.end = paired.end,num.DE = num.DE,top.cutoff = top.cutoff,
                genome = genome, max.tx = max.tx,BPPARAM = BPPARAM,
                read.length = read.length,fragment.length.min = fragment.length.min)
  
  # Quantifying FASTQs
  path.targets <- file.path(dest,'meta/targets.tsv.gz')
  quantifyReads(targets = path.targets,dest = dest,
                genome = genome,workers = workers, keep.fastq = keep.fastq,
                bin.salmon = bin.salmon,index.salmon = index.salmon,opts.salmon = opts.salmon,
                bin.kallisto = bin.kallisto,index.kallisto = index.kallisto,opts.kallisto = opts.kallisto,
                run.salmon = run.salmon, run.kallisto = run.kallisto)
  
  # Running methods
  runDTEMethods(dest = file.path(dest),run.salmon = run.salmon, run.kallisto = run.kallisto)
  
  # Organizing FASTQ files
  path.fastq <- read.delim(path.targets,header = TRUE)
  if (isTRUE(keep.fastq)) {
    message('Copying FASTQ files...')
    dir.create(file.path(dest,'fastq'))
    file.copy(from = unlist(path.fastq),to = file.path(dest,'fastq'))
  }
  file.remove(unlist(path.fastq))
}

runDTEMethods <- function(dest,run.salmon,run.kallisto){
  if(run.salmon){
    message('Running DTE methods with Salmon quantification...')
    dir.create(file.path(dest,'dte-salmon'),recursive = TRUE,showWarnings = FALSE)
    runMethods(path = file.path(dest,'quant-salmon'),dest = file.path(dest,'dte-salmon'),quantifier = 'salmon') 
    if(file.exists(file.path(dest,'dte-salmon','time.tsv'))){
      message('DTE analysis w/ Salmon completed!')
    } else{
      stop('DTE analysis w/ Salmon failed!')
    }
  }
  
  if(run.kallisto){
    message('Running DTE methods with kallisto quantification...')
    dir.create(file.path(dest,'dte-kallisto'),recursive = TRUE,showWarnings = FALSE)
    runMethods(path = file.path(dest,'quant-kallisto'),dest = file.path(dest,'dte-kallisto'),quantifier = 'kallisto') 
    if(file.exists(file.path(dest,'dte-kallisto','time.tsv'))){
      message('DTE analysis w/ kallisto completed!')
    } else{
      stop('DTE analysis w/ kallisto failed!')
    }
  }
  
  return(invisible())
}

quantifyReads <- function(targets,
                          dest,
                          genome,
                          workers,
                          keep.fastq,
                          bin.salmon,
                          index.salmon,
                          opts.salmon,
                          bin.kallisto,
                          index.kallisto,
                          opts.kallisto,
                          run.salmon,
                          run.kallisto){
  
  dir.salmon <- file.path(dest,'quant-salmon')
  dir.kallisto <- file.path(dest,'quant-kallisto')
  
  df.targets <- read.delim(targets,header = TRUE)
  
  if(run.kallisto){
    message('Quantifying reads with kallisto...')
    dir.create(dir.kallisto,recursive = TRUE,showWarnings = FALSE)
    runKallisto(bin = bin.kallisto,
                index = index.kallisto,
                options = opts.kallisto,
                targets = df.targets, 
                dest = dir.kallisto)
  }
  
  if(run.salmon){
    message('Quantifying reads with Salmon...')
    dir.create(dir.salmon,recursive = TRUE,showWarnings = FALSE)
    runSalmon(bin = bin.salmon,
              index = index.salmon,
              options = opts.salmon,
              targets = df.targets,
              dest = dir.salmon)
  }
  
  message('Quantification completed!')
}