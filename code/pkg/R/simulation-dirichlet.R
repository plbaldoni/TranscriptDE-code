#' zipf <- function(n, k = -0.9, a = 3e-5, b = 3e-9, scale = 1e6) {
#'   # This function was originally written by Andy Chen with the following setup:
#'   # k = -0.9, a = 3e-5, b = 3e-10, scale = 1e6
#'   
#'   r <- seq(1,n)
#'   baselineprop <- exp(k * log(r) - r * a - (r ^ 2) * b)
#'   scale * baselineprop / sum(baselineprop)
#' }
#' 
#' #' @importFrom mixtools normmix.sim
#' simulateCPM <- function(n,
#'                         p = 0.7016501, 
#'                         mu = c(-0.9351255,5.1110622), 
#'                         sigma = c(1.165072,2.056817), 
#'                         scale = 1e6){
#'   # This function generates counts per million (CPM) to be used in simulateMean
#'   # I noticed that the distribution of CPM from zipf did not match real data and
#'   # generated many more lowly expressed genes than observed in real data.
#'   
#'   # Simulating mixture
#'   mix <- normmix.sim(n = n,lambda = c(1 - p,p),mu = mu,sigma = sigma)
#'   sort(1e6*(2^(mix)/sum(2^(mix))),decreasing = TRUE)
#' }
#' 
#' simulateMean <- function(x,names,df.bcv = 10, bcv.true = 0.1){
#'   # This function was originally written by Andy Chen with the following setup:
#'   # mu0 <- zipf(n.feat)
#'   # bcv0 <- bcv.true + 1/sqrt(mu0 + 5)
#'   
#'   n.feat <- nrow(x)
#'   n.libs <- length(names)
#'   
#'   mu0 <- simulateCPM(n.feat)
#'   bcv0 <- bcv.true + 1/(mu0 + 2)
#'   
#'   disp <- bcv0 ^ 2 * df.bcv / rchisq(n.feat, df = df.bcv)
#'   
#'   shape <- 1/disp
#'   scale <- mu0/shape
#'   
#'   mu <- matrix(rgamma(n.feat * n.libs, shape = shape, scale = scale),
#'                nrow =  n.feat, ncol =  n.libs)
#'   mu <- mu[sample(1:n.feat),]
#'   
#'   colnames(mu) <- names
#'   
#'   return(cbind(x,disp = disp,mu = mu))
#' }
#' 
#' simulateMean2 <- function(x,lib.sizes,names,df.bcv = 40, bcv.true = 0.2){
#'   # This function was written with the goal of mimicking the simulation setup
#'   # used in the voom paper. Specifically, we use goodTuringProportions to
#'   # estimate baseline abundances, different asymptotic BCV value as well as
#'   # Chi-square degrees of freedom, and a different dispersion trend.
#'   # See ?baselineAbundance_function for goodTuringProportions usage in this
#'   # simulation.
#'   
#'   n.feat <- nrow(x)
#'   n.libs <- length(names)
#'   
#'   # Generating baseline proportions
#'   goodTuring <- get('baselineAbundance_function')
#'   baselineAbundance <- goodTuring(seq_len(n.feat) / (n.feat + 1))
#'   baselineAbundance <- baselineAbundance/sum(baselineAbundance)
#'   
#'   # Generating expected counts
#'   mu0 <- matrix(baselineAbundance,n.feat,1) %*% matrix(lib.sizes,1,n.libs)
#'   
#'   # Dispersion trend
#'   bcv0 <- bcv.true + 1/sqrt(mu0)
#'   
#'   disp <- bcv0 ^ 2 * df.bcv / rchisq(length(bcv0), df = df.bcv)
#'   
#'   # Biological variation
#'   shape <- 1/disp
#'   scale <- mu0/shape
#'   
#'   mu <- matrix(rgamma(n.feat * n.libs, shape = shape, scale = scale),
#'                nrow =  n.feat, ncol =  n.libs)
#'   mu <- mu[sample(1:n.feat),]
#'   
#'   colnames(mu) <- colnames(disp) <- names
#'   
#'   return(cbind(x,disp = disp,mu = mu))
#' }
#' 
#' # txUsage0 <- function(x,r,n) {
#' #   # This function was originally written when simulating reads
#' #   # with a fixed proportion of allocation across genes with the same number of
#' #   # transcripts/gene. We decided to use a Dirichlet distribution to randomly
#' #   # generate random variables to split gene-wise expression among transcripts
#' #   if(n > 4 | r > n) stop('Incompatible number of transcripts')
#' #   if(n == 1) return(x)
#' #   if(n == 2) return(c(0.8,0.2)[r]*x)
#' #   if(n == 3) return(c(0.7,0.2,0.1)[r]*x)
#' #   if(n == 4) return(c(0.6,0.2,0.1,0.1)[r]*x)
#' # }
#' 
#' txUsage <- function(x,r,n,tx) {
#'   # x: vector of gene-wise (gamma) expression 
#'   # r: rank of the transcript
#'   # n: how many transcripts the gene is expressing
#'   # tx: transcript ID
#'   n <- unique(n) # This is a gene-wise function, so all n's must be equal
#'   
#'   if (length(n) > 1) stop('n must be unique within gene')
#'   if (!(n %% 1 == 0)) stop('n must be an integer')
#'   
#'   mat <- as.matrix(x)
#'   gene.expr <- apply(mat,2,unique)
#'   tx.usage <- lapply(seq_along(gene.expr),function(igene){
#'     alpha <- getDirichletConcentration(gene.expr[igene],n = n)
#'     c(rdirichlet(1,alpha))
#'   })
#'   tx.usage <- do.call(cbind,tx.usage)
#'   
#'   mat.tx <- mat * tx.usage
#'   x.tx <- cbind(TranscriptID = tx,as.data.table(mat.tx))
#'   
#'   return(x.tx)
#' }
#' 
#' #' @importFrom data.table setkey copy setnames
#' simulateTPM <- function(x,n.libs,lib.sizes,top.cutoff,num.DE,fc){
#'   
#'   # Generating sample labels
#'   group <- rep(LETTERS[seq_len(length(n.libs))],times = n.libs)
#'   rep <- unlist(lapply(n.libs,seq_len))
#'   group.name <- paste(paste0('group', group), paste0('rep',rep), sep = '_')
#'   
#'   # Simulating gene-wise gamma expression
#'   DT <- as.data.table(x)
#'   DT.gene <- DT[,.(NTx = .N),by = 'GeneID']
#'   DT.gene <- simulateMean2(x = DT.gene,lib.sizes = lib.sizes,names = group.name)
#'   DT <- merge(DT,DT.gene,by = 'GeneID',all.x = TRUE,sort = FALSE)
#'   setkey(DT,'GeneID','TranscriptID','TxRank')
#'   
#'   # Splitting gene-wise expression into transcripts-wise expression
#'   DT.usage <- copy(DT[NTx > 1,c('GeneID','TranscriptID','TxRank','NTx',paste0('mu.',group.name)),with = FALSE])
#'   DT.usage <- DT.usage[,txUsage(.SD,r = TxRank,n = NTx,tx = TranscriptID),by = 'GeneID',.SDcols = paste0('mu.',group.name)][,-1]
#'   setnames(DT.usage,old = paste0('mu.',group.name),new = paste0('mu.',group.name,'_tx'))
#'   
#'   # Merging transcript-wise means back
#'   DT <- merge(DT,DT.usage,all.x = TRUE,by = 'TranscriptID',sort = FALSE)
#'   DT[NTx == 1,(paste0('mu.',group.name,'_tx')) := .SD,.SDcols = paste0('mu.',group.name)]
#'   
#'   # Organizing gene-wise and transcript-wise output
#'   gene.mu <- as.matrix(DT[,paste0('mu.',group.name),with = FALSE])
#'   gene.disp <- as.matrix(DT[,paste0('disp.',group.name),with = FALSE])
#'   tx.mu <- as.matrix(DT[,paste0('mu.',group.name,'_tx'),with = FALSE])
#'   tx.length <- DT$Length
#'   dimnames(gene.mu) <- dimnames(gene.disp) <- dimnames(tx.mu) <- list(DT$TranscriptID,group.name)
#'   names(tx.length) <- DT$TranscriptID
#'   
#'   DT[,(paste0('mu.',group.name)) := NULL]
#'   setnames(DT,new = paste0('mu.',group.name),old = paste0('mu.',group.name,'_tx'))
#'   
#'   # Generating DE (Selecting tx to be DE only among the top 'top.cutoff' expressed tx)
#'   topTx <- order(rowMeans(tx.mu),decreasing = TRUE)[seq(from = 1,length.out = min(top.cutoff,nrow(tx.mu)))]
#'   topTxDE <- sample(topTx,num.DE,replace = FALSE)
#'   
#'   names(topTxDE) <- as.character(sample(1:2,length(topTxDE),replace = TRUE))
#'   whichA <- grepl('A',group.name)
#'   whichB <- grepl('B',group.name)
#'   tx.mu[topTxDE[names(topTxDE) == '1'],whichB] <- fc*tx.mu[topTxDE[names(topTxDE) == '1'],whichB]
#'   tx.mu[topTxDE[names(topTxDE) == '2'],whichA] <- fc*tx.mu[topTxDE[names(topTxDE) == '2'],whichA]
#'   
#'   tpm <- tx.mu / tx.length
#'   tpm <- 1e6 * t(t(tpm) / colSums(tpm))
#'   
#'   de <- vector('numeric',nrow(tpm))
#'   names(de) <- row.names(tpm)
#'   de[topTxDE[names(topTxDE) %in% c('1')]] <- 1L
#'   de[topTxDE[names(topTxDE) %in% c('2')]] <- -1L
#'   
#'   return(list('tpm' = tpm,'de' = de,'gene.disp' = gene.disp, 'gene.mu' = gene.mu))
#' }

# getDirichletConcentration0 <- function(x,n){
#   # This function was originally created to get Dirichlet concentration para-
#   # meters from a Dirichlet regression model fitted on the GSE60450 dataset.
#   # We decided to take on a more simplistic approach to this and use a
#   # parametric function instead, while still using a function resembling the one
#   # observed from real data. Gordon suggested on June 22 2022 to use a single
#   # value of 30 for the concentration parameter. See getDirichletConcentration
#   # function for the alternative (parametric) approach.
#   
#   if(!n%%1==0) stop('n must be an integer')
#   if(n > 5) stop('n must be less than 6')
#   if(!all(x > 0)) stop('x must be a strictly positive')
#   
#   # Originally I attempted to use splines, but I ended up deciding on
#   # fitting a log-linear model.
#   # DF <- as.data.frame(ns(x = log(x),df = 3))
#   # colnames(DF) <- paste0('logTotalAveCountsSplines.',1:3)
#   DF <- data.frame(logTotalAveCounts = log(x))
#   DirichletFit <- get('DirichletFit')
#   alpha <- predict(DirichletFit[[n-1]],newdata = DF,alpha = TRUE)$alpha
#   dimnames(alpha) <- list(1:length(x),paste0('Rank.',1:n))
#   return(alpha)
# }

# DirichletConcentration <- function(x,r, a = 10, b = 40, k = 0.75){
#   # This function assumes a logistic parametric curve on the individual concen-
#   # tration parameters such that the top transcript (r = 1) always dominates
#   # the expression level of the gene and the remaining expression is equally
#   # split amont the other transcripts
#   a.rank <- a * 0.05 ^ (!r == 1)
#   a.rank + 2*(b - a.rank)/(1 + exp(k*x))
# }

# DirichletConcentration2 <- function(x, r, alpha0 = 40, ...){
#   # This function assumes a decaying function of the type f(r) = 0.9^r for the 
#   # split proportion, in which r is the rank of transcript. These proportions
#   # are, in turn, translated to individual Dirichlet concentration parameters 
#   # that sum up to alpha0 (from empirical data; Cory's RNA-seq mouse 
#   # deeply sequenced MNT WT/KO dataset from May 2022)
#   p <- 0.9^(max(r) - 1)
#   alpha0 * c(p, rep((1 - p) / (max(r) - 1), (max(r) - 1)))
# }

# getDirichletConcentration <- function(x,n){
#   if (!n %% 1 == 0)
#     stop('n must be an integer')
#   if (!all(x > 0))
#     stop('x must be strictly positive')
#   
#   alpha <- lapply(X = log(x),FUN = DirichletConcentration2,r = seq_len(n))
#   alpha <- do.call(rbind,alpha)
#   dimnames(alpha) <- list(1:length(x),paste0('Rank.',1:n))
#   return(alpha)
# }