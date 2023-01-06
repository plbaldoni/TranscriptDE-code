library(Rsubread)
library(BiocParallel)
library(devtools)
library(readr)

load_all("../pkg/")

### Parameters

dest_main <- '../../output/quasi_poisson'
fasta <- "../../data/annotation/mm39/gencode.vM27.transcripts.fa.gz"
tmpdir <- tempdir(check = TRUE)
workers <- 25
num.DE <- 3000
fc <- 2
max.tx <- Inf
genome <- 'mm39'
top.cutoff <- 30000
n.libs <- 10
lib.sizes <- rep(50e6,sum(n.libs))
paired.end <- TRUE
keep.fastq <- FALSE
read.length <- 75
fragment.length.min <- 150
bin.salmon = "/stornext/General/data/academic/lab_smyth/baldoni.p/software/salmon-1.9.0_linux_x86_64/bin/salmon"
index.salmon = file.path("/stornext/General/data/academic/lab_smyth/baldoni.p/software/SalmonIndex",genome,"transcripts_index/")
opts.salmon =  paste('-p',workers,'-l A --numBootstraps 100 --validateMappings')
bin.kallisto = "/stornext/General/data/academic/lab_smyth/baldoni.p/software/kallisto/kallisto"
index.kallisto = file.path("/stornext/General/data/academic/lab_smyth/baldoni.p/software/kallistoIndex",genome,"transcripts_index")
opts.kallisto = paste0('--bootstrap-samples=100 --threads=',workers)
run.salmon = FALSE
run.kallisto = TRUE

### Simulation

# Setting up parallel computing
BPPARAM <- MulticoreParam(workers = workers, progressbar = TRUE)
register(BPPARAM = BPPARAM)

# Preparing directories
dir.create(tmpdir,recursive = TRUE,showWarnings = FALSE)
tmpdir <- normalizePath(tmpdir)

# simulateFASTQ

# Scanning fasta
message('Preparing fasta file for simulation...')
contigs <- readFasta(fasta)

# Subsetting transcripts according to curated data and bringing metadata
message('Subsetting transcripts...')
contigs.subset <- selectTx(contigs,genome,max.tx)

# Generating sample labels
group <- rep(LETTERS[seq_len(length(n.libs))],times = n.libs)
rep <- unlist(lapply(n.libs,seq_len))
group.name <- paste(paste0('group', group), paste0('rep',rep), sep = '_')

# simulateExpr
n.groups <- length(n.libs)
n.samples <- sum(n.libs)
n.feat <- nrow(contigs.subset)
df.bcv = 40
bcv.true = 0.2

# Generating baseline proportions
goodTuring <- get('baselineAbundance_function')
baselineAbundance <- goodTuring(seq_len(n.feat) / (n.feat + 1))
baselineAbundance <- baselineAbundance/sum(baselineAbundance)
baselineAbundancePerGroup <- matrix(baselineAbundance,ncol = n.groups,nrow = n.feat)

# Generating expected counts
mu0 <- lapply(seq_len(n.groups),function(x){
  size <- lib.sizes[seq(1 + n.libs[x] * (x - 1), n.libs[x] * x)]
  matrix(baselineAbundancePerGroup[,x],n.feat,1) %*% matrix(size,1,length(size))
})
mu0 <- do.call(cbind,mu0)

# Generating random noise around dispersion trend.
chisq <- lapply(seq_len(n.groups),function(x){
  rv <- df.bcv / rchisq(n.feat, df = df.bcv)
  matrix(rv,ncol = n.libs[x],nrow = n.feat)
})
chisq <- do.call(cbind,chisq)

# Biological variation and NULL Dispersion trend on the mean
bcv0 <- bcv.true
disp <- bcv0 ^ 2 * chisq

# Biological variation
shape <- 1/disp
scale <- mu0/shape

# Bio & tech reps simulation

expr.tech <- matrix(rgamma(n.feat, shape = shape, scale = scale),
                    nrow =  n.feat, ncol =  n.samples,byrow = FALSE)

expr.bio <- matrix(rgamma(n.feat*n.samples, shape = shape, scale = scale),
                   nrow =  n.feat, ncol =  n.samples,byrow = FALSE)

for(i in c('bio', 'tech')) {
  expr <- get(paste0('expr.',i))
  dest <- file.path(dest_main,i)

  # Checking if simulation has already been run
  dir.create(dest,recursive = TRUE,showWarnings = FALSE)
  dest <- normalizePath(dest)
  dest = file.path(dest,'meta')
  dir.create(dest,showWarnings = FALSE,recursive = TRUE)
  out.files <- file.path(dest,c('counts.tsv.gz','targets.tsv.gz'))
  names(out.files) <- c('counts','targets')


  # Shuffling matrices (to be matched with reference data)
  key <- sample(n.feat,replace = FALSE)
  expr <- expr[key,]
  mu0 <- mu0[key,]
  disp <- disp[key,]

  trExpr <- list('expr' = expr,'mu' = mu0,'disp' = disp,'de' = de)

  # Generating TPM values
  tpm <- trExpr$expr / contigs.subset$Length
  tpm <- 1e6 * t(t(tpm) / colSums(tpm))

  # Organizing contigs.subset
  dimnames(tpm) <- dimnames(trExpr$expr) <- dimnames(trExpr$mu) <- dimnames(trExpr$disp) <- list(contigs.subset$TranscriptID,group.name)

  # Merging contigs.subset values to contigs
  tpm.contigs <- expr.contigs <- mu.contigs <- disp.contigs <- matrix(NA, nrow(contigs), sum(n.libs))

  dimnames(tpm.contigs) <- dimnames(expr.contigs) <- dimnames(mu.contigs) <- dimnames(disp.contigs) <- list(contigs$TranscriptID,group.name)

  key <- match(contigs.subset$TranscriptID, contigs$TranscriptID)
  tpm.contigs[key,] <- tpm
  expr.contigs[key,] <- trExpr$expr
  mu.contigs[key,] <- trExpr$mu
  disp.contigs[key,] <- trExpr$disp

  txTPM <- list('tpm' = tpm.contigs,'expr' = expr.contigs,'mu' = mu.contigs,'disp' = disp.contigs)

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

  out <- cbind(contigs[,c("TranscriptID","Length",'GeneID')],mat)
  rownames(out) <- NULL
  write_tsv(x = out,file = out.files['counts'],col_names = TRUE,quote = 'none')

  targets <- data.frame('R1' = out.fastq[grepl("R1.fastq.gz",out.fastq)])
  if (isTRUE(paired.end)) {
    targets$R2 <- out.fastq[grepl("R2.fastq.gz",out.fastq)]
  }
  write_tsv(x = targets,file = out.files['targets'],col_names = TRUE,quote = 'none')

  # Quantifying FASTQs
  dest <- file.path(dest_main,i)
  path.targets <- file.path(dest,'meta/targets.tsv.gz')
  quantifyReads(targets = path.targets,dest = dest,
                genome = genome,workers = workers, keep.fastq = keep.fastq,
                bin.salmon = bin.salmon,index.salmon = index.salmon,opts.salmon = opts.salmon,
                bin.kallisto = bin.kallisto,index.kallisto = index.kallisto,opts.kallisto = opts.kallisto,
                run.salmon = run.salmon, run.kallisto = run.kallisto)

  # Organizing FASTQ files
  path.fastq <- read.delim(path.targets,header = TRUE)
  if (isTRUE(keep.fastq)) {
    message('Copying FASTQ files...')
    dir.create(file.path(dest,'fastq'))
    file.copy(from = unlist(path.fastq),to = file.path(dest,'fastq'))
  }
  file.remove(unlist(path.fastq))
}
