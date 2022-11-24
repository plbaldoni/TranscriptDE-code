#' @importFrom sleuth sleuth_prep sleuth_fit sleuth_lrt sleuth_wt sleuth_results
runSleuth <- function(targets,test,quantifier){
  se <- 
    sleuth_prep(sample_to_covariates = targets,full_model = ~ group,num_cores = 1)
  
  se <- sleuth_fit(obj = se, fit_name = 'full')
  
  if (test == 'lrt') {
    test.label <- 'reduced:full'
    se <- sleuth_fit(obj = se,formula = ~ 1, fit_name = 'reduced')
    se <- sleuth_lrt(obj = se,null_model = 'reduced',alt_model = 'full')
  } else{
    test.label <- 'groupB'
    se <- sleuth_wt(obj = se,which_beta = test.label,which_model = 'full')
  }
  
  out <-
    sleuth_results(obj = se,test = test.label,test_type = test,show_all = FALSE)
  
  out <- cbind('feature' = out$target_id,out)
  out$target_id <- NULL
  
  return(out)
}

#' @importFrom tximeta tximeta
#' @importFrom fishpond scaleInfReps labelKeep swish
#' @importFrom S4Vectors mcols
runSwish <- function(targets,quantifier){
  if (grepl('salmon', quantifier)) {
    type <- 'salmon'
    targets$files <- file.path(targets$path,'quant.sf')
  } else{
    type <- 'kallisto'
    targets$files <- file.path(targets$path,'abundance.h5')
  }
  
  se <- tximeta(coldata = targets,type = type)
  se <- scaleInfReps(se)
  se <- labelKeep(se)
  se <- se[mcols(se)$keep,]
  se <- swish(y = se, x = "group")
  out <- as.data.frame(mcols(se))
  cbind('feature' = rownames(out),out)
}

#' @importFrom edgeR catchKallisto catchSalmon DGEList filterByExpr
#' @importFrom edgeR calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
runEdgeR <- function(targets,quantifier,scaled){
  if (grepl('salmon',quantifier)) {
    se <- catchSalmon(paths = targets$path) 
  } else{
    se <- catchKallisto(paths = targets$path)
  }
  
  if (scaled) {
    dge <- DGEList(counts = se$counts/se$annotation$Overdispersion, 
                   samples = targets,
                   group = targets$group,
                   genes = se$annotation)
  } else{
    dge <- DGEList(counts = se$counts, 
                   samples = targets,
                   group = targets$group,
                   genes = se$annotation)
  }
  
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  design <- model.matrix(~group,data = dge$samples)
  dge <- estimateDisp(dge,design)
  fit <- glmQLFit(dge,design)
  qlf <- glmQLFTest(fit,coef = 2)
  out <- topTags(qlf,n = Inf)
  
  cbind('feature' = rownames(out$table),out$table)
}

callMethods <- function(targets,quantifier){

  res <- list()
  time <- list()
  
  time[['sleuth-lrt']] <-
    system.time({res[['sleuth-lrt']] <- runSleuth(targets = targets, quantifier = quantifier, test = 'lrt')})
  
  time[['sleuth-wt']] <-
    system.time({res[['sleuth-wt']] <- runSleuth(targets = targets, quantifier = quantifier, test = 'wt')})
  
  time[['swish']] <-
    system.time({res[['swish']] <- runSwish(targets = targets, quantifier = quantifier)})
  
  time[['edger-scaled']] <-
    system.time({res[['edger-scaled']] <- runEdgeR(targets = targets, quantifier = quantifier,scaled = TRUE)})
  
  time[['edger-raw']] <-
    system.time({res[['edger-raw']] <- runEdgeR(targets = targets, quantifier = quantifier,scaled = FALSE)})
  
  for (meth.name in names(res)) {
    res[[meth.name]]$feature <- strsplit2(res[[meth.name]]$feature,"\\|")[,1]
  }
  
  
  time <- as.data.frame(do.call(rbind,time))
  time <- cbind('method' = rownames(time),time)
  
  return(list('res' = res, 'time' = time))
}

runMethods <- function(path,dest,quantifier){
  
  path <- normalizePath(path)
  sample.names <- basename(list.dirs(path,full.names = TRUE,recursive = FALSE))
  sample.names.split <- strsplit2(sample.names,'_')
  
  dir.create(path = dest,recursive = TRUE,showWarnings = FALSE)
  dest <- normalizePath(dest)
  if (file.exists(file.path(dest, 'targets.tsv'))) {
    message('Files are already processed. Skipping them and moving on...')
    return(invisible(NULL))
  }

  targets <- data.frame(group = gsub('group','',sample.names.split[,1]),
                        replicate = gsub('rep','',sample.names.split[,2]))
  
  targets$group <- relevel(as.factor(targets$group),ref = 'A')
  targets$names <- sample.names
  targets$sample <- targets$names
  targets$path <- file.path(path,targets$names)

  out <- callMethods(targets,quantifier)
  
  for (meth.name in names(out$res)) {
    write_tsv(x = out$res[[meth.name]],
              file = file.path(dest,paste0(meth.name,'.tsv.gz')),
              col_names = TRUE,quote = 'none')
  }
  write_tsv(x = out$time,file = file.path(dest,'time.tsv'),
            col_names = TRUE,quote = 'none')
  write_tsv(x = targets,file = file.path(dest,'targets.tsv'),
            col_names = TRUE,quote = 'none')
  
  return(invisible(NULL))
}
