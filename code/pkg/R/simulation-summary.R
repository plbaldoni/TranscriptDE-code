methodsNames <- function(){
  method <- c('edger-raw','edger-scaled','sleuth-lrt','sleuth-wt','swish')
  labels <- c('edgeR-raw','edgeR-scaled','sleuth-LRT','sleuth-Wald','Swish')
  color <- okabe_ito(length(method))

  names(method) <- names(color) <- labels
  return(list(labels = labels,method = method,color = color))
}

roundPretty <- function(x,digits = 1){
  formatC(round(x,digits),digits = digits,format = 'f')
}

#' @importFrom data.table data.table fread
loadResults <- function(path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,simulation,quantifier){

  meth <- methodsNames()
  path.time <- file.path(path,'time.tsv')
  path.method <- file.path(path,paste0(meth$method,'.tsv.gz'))
  names(path.method) <- meth$labels

  dt.scenario <- data.table('Genome' = genome,
                            'Length' = len,
                            'FC' = fc,
                            'Reads' = read,
                            'TxPerGene' = tx.per.gene,
                            'Scenario' = scenario,
                            'LibsPerGroup' = libs.per.group,
                            'Simulation' = simulation,
                            'Quantifier' = quantifier)

  dt.results <- lapply(seq_along(path.method),function(i){
    dt <- data.table('Method' = names(path.method)[i], fread(path.method[i],header = TRUE))
    setnames(x = dt,
             old = c('pvalue','pval','qvalue','qval','feature'),
             new = c('PValue','PValue','FDR','FDR','TranscriptID'),skip_absent = TRUE)

    dt <- dt[,c('Method','TranscriptID','PValue','FDR')]
    return(dt)
  })

  dt.results <- do.call(rbind,dt.results)
  dt.results <- cbind(dt.scenario,dt.results)

  dt.time <- fread(file = path.time,header = TRUE)
  dt.time <- cbind(dt.scenario,dt.time[,c('method','elapsed')])
  setnames(x = dt.time,old = c('method','elapsed'),new = c('Method','Time'))
  dt.time$Method <- names(meth$method)[match(dt.time$Method,meth$method)]

  out <- list('results' = dt.results,'time' = dt.time)

  return(out)
}

loadMetadata <- function(path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,simulation){
  path.counts <- file.path(path,'counts.tsv.gz')

  dt.scenario <- data.table('Genome' = genome,
                            'Length' = len,
                            'FC' = fc,
                            'Reads' = read,
                            'TxPerGene' = tx.per.gene,
                            'Scenario' = scenario,
                            'LibsPerGroup' = libs.per.group,
                            'Simulation' = simulation)

  dt.metadata <- fread(path.counts,header = TRUE)

  dt.transcript <- dt.metadata[,c('TranscriptID')]

  dt.metadata <- cbind(dt.scenario,dt.metadata[!dt.metadata$status == 0, c('TranscriptID', 'status')])

  # If fold-change = 1, status should be 0
  if (fc == 1) dt.metadata$status <- 0

  out <- list('simulation' = dt.metadata,'transcript' = dt.transcript)

  return(out)
}

aggregateScenario <- function(path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,quantifier,nsim){

  subpath <- paste0('simulation-',seq_len(nsim))

  ls.results <- lapply(seq_len(nsim),function(x){
    res.path <- file.path(path,subpath[x],paste0('dte-',quantifier))
    loadResults(res.path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,x,quantifier)
  })

  ls.metadata <- lapply(seq_len(nsim),function(x){
    meta.path <- file.path(path,subpath[x],'meta')
    loadMetadata(meta.path,genome,len,fc,read,tx.per.gene,scenario,libs.per.group,x)
  })

  dt.results <- do.call(rbind,lapply(ls.results,function(x){x[['results']]}))
  dt.time <- as.data.table(do.call(rbind,lapply(ls.results,function(x){x[['time']]})))
  dt.simulation <- do.call(rbind,lapply(ls.metadata,function(x){x[['simulation']]}))
  dt.transcript <- ls.metadata[[1]]$transcript

  out <- list('results' = dt.results,'time' = dt.time,'simulation' = dt.simulation,'features' = dt.transcript)
  return(out)
}

computeMetrics <- function(x,simulation,features,fdr,alpha){

  TranscriptID.DE <- x$TranscriptID[x$FDR < fdr]
  n <- length(x$TranscriptID)
  n.lt.alpha <- sum(x$PValue < alpha)

  if (length(TranscriptID.DE) > 0) {
    call.DE <- data.table(TranscriptID = TranscriptID.DE,call = 1)

    truth.DE <- simulation[Genome == x$Genome & Length == x$Length & FC == x$FC &
                             Reads == x$Reads & TxPerGene == x$TxPerGene &
                             Scenario == x$Scenario & LibsPerGroup == x$LibsPerGroup &
                             Simulation == x$Simulation, c('TranscriptID','status')]

    tb.DE <- merge(features,truth.DE,by = 'TranscriptID',all.x = TRUE)
    tb.DE <- merge(tb.DE,call.DE,by = 'TranscriptID',all.x = TRUE)

    tb.DE[is.na(status),status := 0]
    tb.DE[, status := abs(status)]
    tb.DE[is.na(call),call := 0]

    tb.DE$status <- factor(tb.DE$status,levels = c(0,1))
    tb.DE$call <- factor(tb.DE$call,levels = c(0,1))

    tb.results <- tb.DE[,table(status,call)]

    total.de <- sum(tb.results["1",])

    out <- list('N' = n,
                'N.ALPHA' = n.lt.alpha,
                'TP' = tb.results["1","1"],
                'FP' = tb.results["0","1"],
                'FPR' = tb.results["0","1"]/sum(tb.results["0",]),
                'FDR' = tb.results["0","1"]/sum(tb.results[,"1"]),
                'TPR' = ifelse(total.de == 0,NA,tb.results["1","1"]/total.de))
  } else{
    out <- list('N' = n,'N.ALPHA' = n.lt.alpha,'TP' = 0,'FP' = 0,'FPR' = 0,'FDR' = 0,'TPR' = 0)
  }
  return(lapply(out,as.double))
}

computeFDRCurve <- function(x,simulation,features,fdr,seq.n){

  truth.DE <- simulation[Genome == x$Genome & Length == x$Length & FC == x$FC &
                           Reads == x$Reads & TxPerGene == x$TxPerGene &
                           Scenario == x$Scenario & LibsPerGroup == x$LibsPerGroup &
                           Simulation == x$Simulation, c('TranscriptID','status')]

  tb.DE <- merge(features,truth.DE,by = 'TranscriptID',all.x = TRUE)
  tb.DE[is.na(status),status := 0]
  tb.DE[, status := abs(status)]

  feature.DE <- data.table(TranscriptID = x$TranscriptID,FDR = x$FDR,call = 1)

  tb.DE <- merge(tb.DE,feature.DE,by = 'TranscriptID',all.x = TRUE)
  tb.DE[is.na(call),call := 0]
  tb.DE$status <- factor(tb.DE$status,levels = c(0,1))
  tb.DE$call <- factor(tb.DE$call,levels = c(0,1))
  tb.DE <- tb.DE[order(FDR),]

  out <- lapply(seq.n,function(w){
    tb.results <- tb.DE[seq(1,w),][,table(status,call)]
    return(tb.results["0","1"])
  })

  names(out) <- paste0('n.',seq.n)

  return(out)
}

#' @importFrom thematic okabe_ito
#' @importFrom ggplot2 ggplot geom_line theme_bw scale_color_manual
#' @importFrom ggplot2 scale_x_continuous theme element_blank labs aes
#' @importFrom ggplot2 scale_y_continuous geom_abline facet_grid vars
plotFDRCurve <- function(x,max.n,base_size = 8){

  meth <- methodsNames()

  plot <- ggplot(x,aes(x = N,y = FDR,color = Method,group = Method)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(TxPerGene),scales = 'free_y') +
    geom_line(size = 0.75) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_color_manual(values = meth$color) +
    # scale_y_continuous(limits = c(0,1300)) +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"line"),
          legend.position = 'top',
          panel.grid = element_blank(),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size)) +
    labs(y = 'False discoveries',x = 'Transcripts chosen')

  return(plot)
}

#' @importFrom ggplot2 geom_col geom_text scale_fill_brewer
plotPowerBars <- function(x,fdr,max.n,base_size = 8){

  sub.byvar <- colnames(x)[-which(colnames(x) %in% c('P.SIG','TP','FP'))]

  gap <- 0.05*max(x$TP + x$FP)

  x.melt <- melt(x,id.vars = sub.byvar,
                 measure.vars = c('TP','FP'),
                 variable.name = 'Type',
                 value.name = 'Value')
  x.melt$Type <- factor(x.melt$Type,levels = c('FP','TP'),labels = c('False Positive','True Positive'))

  plot <- ggplot(x.melt,aes(x = Method,y = Value,fill = Type)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(TxPerGene)) +
    geom_col() +
    geom_text(aes(x = Method,y = (TP + FP) + gap,
                  label = roundPretty(ifelse((FP+TP) == 0,NA,100*FP/(FP+TP)),1)),
              inherit.aes = FALSE,data = x,vjust = 0,size = base_size/.pt) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_fill_brewer(palette = 'Set1') +
    labs(x = NULL,y = paste('# DE Transcripts at FDR <',roundPretty(fdr,2))) +
    scale_y_continuous(limits = c(0,max.n)) +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = 'top',
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 90,colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.key.size = unit(0.75,"line"))

  return(plot)
}

plotType1Error <- function(x,alpha,base_size = 8){

  sub.byvar <- colnames(x)[-which(colnames(x) %in% c('P.SIG','TP','FP'))]

  x.melt <- melt(x,id.vars = sub.byvar,
                 measure.vars = c('P.SIG'),variable.name = 'Type',value.name = 'Value')

  plot <- ggplot(x.melt,aes(x = Method,y = Value)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(TxPerGene)) +
    geom_col() +
    theme_bw(base_size = base_size,base_family = 'sans') +
    geom_hline(yintercept = alpha,color = 'red',linetype = 'dashed') +
    labs(x = NULL,y = paste('Proportion of transcripts with p-value <',roundPretty(alpha,2))) +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"line"),
          legend.position = 'top',
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 90),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size))

  return(plot)
}

summarizePValue <- function(x,byvar,step = 0.05){
  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  cut.sq <- seq(step,1 - step,by = step)
  cut.match <- c(0,cut.sq) + step/2
  names(cut.match) <- paste0('(',roundPretty(c(0,cut.sq),2),'-',roundPretty(c(cut.sq,1),2),']')
  n.groups <- length(cut.match)

  x.sub <- copy(x)
  x.sub[,PValue := cut(x = PValue,breaks = c(-Inf,cut.sq,Inf),labels = names(cut.match))]

  x.sub.method <- x.sub[,list(N = .N),by = byvar]

  table <- x.sub[,list(N.cat = .N),by = c(byvar,'PValue')]
  table <- merge(table,x.sub.method,by = byvar,all.x = TRUE)
  table$PValue <- factor(table$PValue,levels = names(cut.match))

  table <- table[,list(Density.Avg = mean(n.groups*N.cat/N)),by = c(sub.byvar,'PValue')]

  table$PValue.Midpoint <- cut.match[match(table$PValue,names(cut.match))]

  return(table)
}

#' @importFrom ggplot2 facet_wrap geom_histogram geom_hline scale_x_discrete rel
plotPValues <- function(x,base_size = 8){
  plot <- ggplot(data = x,aes(x = PValue,y = Density.Avg)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(Method)) +
    geom_col(col = 'black',fill = 'grey',position = position_dodge(0.9),width = 0.8) +
    geom_hline(yintercept = 1,col = 'red',linetype = 'dashed') +
    theme_bw(base_size = base_size,base_family = 'sans') +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          legend.key.size = unit(0.75,"line"),
          legend.position = 'top',
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size)) +
    labs(x = 'P-values',y = 'Density')
  return(plot)
}

#' @importFrom ggplot2 geom_bar position_dodge element_text
plotTime <- function(x,base_size = 8){
  plot <- ggplot(data = x,aes(x = Method,y = Time)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(TxPerGene)) +
    geom_bar(stat = 'identity',position = position_dodge()) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_y_continuous(limits = c(0,max(5,max(x$Time)))) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          legend.position = 'top',
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 90,colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size),
          strip.text = element_text(colour = 'black',size = base_size)) +
    labs(y = 'Time (min)',x = NULL)
  return(plot)
}

#' @importFrom data.table melt
summarizeFDRCurve <- function(x,byvar){

  cnames <- colnames(x)

  sub.cnames <- cnames[grepl('n\\.',cnames)]
  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  x.mean <- x[,lapply(.SD,mean),by = sub.byvar,.SDcols = sub.cnames]

  table <- melt(data = x.mean,id.vars = sub.byvar,variable.name = 'N',value.name = 'FDR')

  table[,N := as.numeric(gsub('n\\.','',N))]

  return(table)
}

summarizeMetrics <- function(x,byvar){

  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  table <- x[,.(P.SIG = mean(N.ALPHA/N),TP = mean(TP),FP = mean(FP)),sub.byvar]

  return(table)
}

summarizeTime <- function(x,byvar){

  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  table <- x[,.(Time = mean(Time/60)),sub.byvar]

  return(table)
}

summarizeQQ <- function(x,byvar,step = 0.001){

  sub.byvar <- byvar[!grepl('Simulation',byvar)]

  x.quant <- x[, list(q.sample = quantile(PValue,probs = seq(0,1,length.out = .N)),
                      q.theory = seq(0,1,length.out = .N)),by = byvar]

  quant.sq <- seq(step,1 - step,by = step)
  quant.match <- c(0,quant.sq) + step/2
  names(quant.match) <- paste0('(',c(0,quant.sq),'-',c(quant.sq,1),']')

  x.quant[,Q.Theory.Cat := cut(x = q.theory,breaks = c(-Inf,quant.sq,Inf),labels = names(quant.match))]

  table <- x.quant[,.(Q.Sample.Avg = mean(q.sample)),by = c(sub.byvar,'Q.Theory.Cat')]

  table$Q.Theory.Midpoint <- quant.match[match(table$Q.Theory.Cat,names(quant.match))]

  return(table)
}

plotQQPlot <- function(x,base_size = 8){
  meth <- methodsNames()

  plot <- ggplot(x,
                 aes(x = Q.Theory.Midpoint,y = Q.Sample.Avg,color = Method,group = Method)) +
    facet_grid(rows = vars(LibsPerGroup),cols = vars(TxPerGene)) +
    # geom_abline(intercept = 0,slope = 1,colour = 'black',linetype = 'dashed') +
    geom_line(alpha = 0,size = 0) +
    geom_point(pch = '.',size = 2) +
    theme_bw(base_size = base_size,base_family = 'sans') +
    scale_color_manual(values = meth$color) +
    scale_x_continuous(breaks = c(0,0.5,1),labels = c(0,0.5,1)) +
    theme(strip.text = element_text(colour = 'black',size = base_size),
          legend.background = element_rect(fill = alpha('white', 0)),
          legend.text = element_text(size = base_size),
          legend.title = element_blank(),
          legend.key.size = unit(2,"line"),
          legend.position = 'top',
          panel.grid = element_blank(),
          axis.text = element_text(colour = 'black',size = base_size),
          axis.title = element_text(colour = 'black',size = base_size)) +
    labs(y = 'Sample Quantiles',x = 'Theoretical Quantiles')

  return(plot)
}

summarizeOverdispersion <- function(path, genome, len, fc, read, tx.per.gene,
                                    scenario, libs.per.group, quantifier, nsim){

  subpath <- paste0('simulation-',seq_len(nsim))
  catchFunction <- get(ifelse(quantifier == 'salmon','catchSalmon','catchKallisto'))

  ls.results <- lapply(seq_len(nsim),function(x){
    res.path <- file.path(path,subpath[x],paste0('quant-',quantifier))
    meta.path <- file.path(path,subpath[x],'meta/counts.tsv.gz')

    meta <- fread(meta.path)

    catch <- catchFunction(list.dirs(res.path,recursive = FALSE))
    rownames(catch$annotation) <- strsplit2(rownames(catch$annotation),"\\|")[,1]

    keep <- rownames(catch$annotation) %in% meta$TranscriptID[!is.na(meta$status)]

    out <- catch$annotation$Overdispersion[keep]
    return(out)
  })

  res <- log10(unlist(ls.results))

  out <- data.table('Genome' = genome,
                    'Length' = len,
                    'FC' = fc,
                    'Reads' = read,
                    'TxPerGene' = tx.per.gene,
                    'Scenario' = scenario,
                    'LibsPerGroup' = libs.per.group,
                    'Quantifier' = quantifier,
                    'Mean' = mean(res),
                    'SD' = sd(res),
                    '2.5Pct' = quantile(res,0.025),
                    '25Pct' = quantile(res,0.25),
                    '50Pct' = quantile(res,0.5),
                    '75Pct' = quantile(res,0.75),
                    '97.5Pct' = quantile(res,0.975))

  return(out)
}

#' @importFrom data.table fwrite
summarizeQuantification <- function(path,dest,genome,fc,read,len,
                                    tx.per.gene,scenario,libs.per.group,quantifier,
                                    nsim = 20, fdr = 0.05, seq.n = seq(100,3000,100),alpha = fdr){

  byvar <- c('Genome','Length','FC','Reads','TxPerGene','Scenario','LibsPerGroup','Quantifier','Method','Simulation')

  res <- aggregateScenario(path, genome, len, fc, read, tx.per.gene, scenario, libs.per.group, quantifier, nsim)

  table.metrics <- res$results[,computeMetrics(c(.BY,.SD),simulation = res$simulation,features = res$features,fdr = fdr,alpha = alpha),by = byvar]
  table.fdr <- res$results[,computeFDRCurve(c(.BY,.SD),simulation = res$simulation,features = res$features,fdr = fdr,seq.n = seq.n),by = byvar]
  table.overdispersion <- summarizeOverdispersion(path, genome, len, fc, read, tx.per.gene, scenario, libs.per.group, quantifier, nsim)

  out <- list('fdr' = summarizeFDRCurve(table.fdr,byvar),
              'metrics' = summarizeMetrics(table.metrics,byvar),
              'time' = summarizeTime(res$time,byvar),
              'quantile' = summarizeQQ(res$results,byvar),
              'pvalue' = summarizePValue(res$results,byvar),
              'overdispersion' = table.overdispersion)

  # plotFDRCurve(out$fdr,max.n = max(seq.n))
  # plotPowerBars(out$metrics,fdr,max.n)
  # plotType1Error(out$metrics,alpha)
  # plotTime(out$time)
  # plotQQPlot(out$quantile)
  # plotPValues(out$pvalue)

  dir.create(dest,showWarnings = FALSE,recursive = TRUE)
  lapply(names(out),function(x){fwrite(x = out[[x]], file = file.path(dest,paste0(x,'.tsv.gz')),quote = FALSE,sep = '\t')})

  return(invisible())
}

summarizeScenario <- function(x,table,path,dest){
  dt <- as.character(table[x,])
  names(dt) <- colnames(table)
  table.names = c('fdr','metrics','time','quantile','pvalue','overdispersion')

  in.path <- file.path(path,do.call(file.path,as.list(dt)))
  out.path <- file.path(dest,do.call(file.path,as.list(dt)))

  # Check is simulation directory exists
  if (!dir.exists(in.path)) return(invisible())

  # Summarizing results with Salmon
  if (!all(file.exists(file.path(out.path,'dte-salmon',paste0(table.names,'.tsv.gz'))))){
    # Verbose
    summarizeQuantification(path = in.path,quantifier = 'salmon',
                            dest = file.path(out.path,'dte-salmon'),
                            genome = dt['genome'],fc = dt['fc'],read = dt['read'],
                            tx.per.gene = dt['tx.per.gene'],scenario = dt['scenario'],
                            libs.per.group = dt['libs.per.group'],len = dt['len'])
  }

  # Summarizing results with kallisto
  if (!all(file.exists(file.path(out.path,'dte-kallisto',paste0(table.names,'.tsv.gz'))))){
    # Verbose
    summarizeQuantification(path = in.path,quantifier = 'kallisto',
                            dest = file.path(out.path,'dte-kallisto'),
                            genome = dt['genome'],fc = dt['fc'],read = dt['read'],
                            tx.per.gene = dt['tx.per.gene'],scenario = dt['scenario'],
                            libs.per.group = dt['libs.per.group'],len = dt['len'])
  }
}

summarizeSimulation <- function(path,
                                dest,
                                genome = c('mm39'),
                                len = c(50,75,100,125,150),
                                fc = c(1,2),
                                read = c('single-end','paired-end'),
                                tx.per.gene = c(2,3,4,5,9999),
                                scenario = c('balanced','unbalanced'),
                                libs.per.group = c(3,5,10), workers = 1){

  path <- normalizePath(path)
  dir.create(dest,showWarnings = FALSE,recursive = TRUE)
  dest <- normalizePath(dest)

  BPPARAM <- MulticoreParam(workers = workers,progressbar = TRUE)
  register(BPPARAM)

  dt.scenario <- expand.grid('genome' = genome,
                             'len' = paste0('readlen-',len),
                             'fc' = paste0('fc',fc),
                             'read' = read,
                             'tx.per.gene' = paste0(tx.per.gene,'TxPerGene'),
                             'scenario' = scenario,
                             'libs.per.group' = paste0(libs.per.group,'libsPerGroup'),
                             stringsAsFactors = FALSE)

  bplapply(seq_len(nrow(dt.scenario)),summarizeScenario,
           table = dt.scenario,dest = dest,path = path,BPPARAM = BPPARAM)

  return(invisible())
}

#' @importFrom kableExtra kbl add_header_above kable_styling pack_rows cell_spec
#' @importFrom kableExtra landscape collapse_rows
tabulateMetrics <- function(x,cap,
                            seq.len = seq(50,150,25),
                            lib.group = c(3,5,10),
                            lib.size = c('50M','25/100M'),
                            color = TRUE,
                            font_size = NULL,
                            color.fdr = 0.05,
                            format = 'latex',...){

  dt <- copy(x)

  methods <- methodsNames()$labels

  dt$Length %<>% factor(levels = paste0(seq.len,'bp'))
  dt$LibsPerGroup %<>% mapvalues(from = paste0("#Lib/Group = ",lib.group),to = lib.group)
  dt$Scenario %<>% mapvalues(from = c('balanced','unbalanced'),to = lib.size) %<>% factor(levels = lib.size)

  dt[,A.Power := TP/3000]
  dt[,B.FDR := ifelse((FP+TP) == 0,NA,FP/(FP+TP))]

  dt.dcast <- dcast(dt, Reads + LibsPerGroup + Scenario + Length ~ Method,value.var = c('A.Power','B.FDR'))

  dt.dcast <- dt.dcast[order(Reads,LibsPerGroup,Scenario,Length),]
  dt.dcast <- dt.dcast[,c('Reads','LibsPerGroup','Scenario','Length',paste0('A.Power_',methods),paste0('B.FDR_',methods)),with = FALSE]

  dt.dcast[,c(paste0('A.Power_',methods),paste0('B.FDR_',methods)) := lapply(.SD,roundPretty,digits = 3),.SDcols = c(paste0('A.Power_',methods),paste0('B.FDR_',methods))]

  dt.dcast.color <- copy(dt.dcast)

  if(color == TRUE){
    mat.power <- as.matrix(dt.dcast[,paste0('A.Power_',methods),with = FALSE])
    class(mat.power) <- 'numeric'
    mat.fdr <- as.matrix(dt.dcast[,paste0('B.FDR_',methods),with = FALSE])
    class(mat.fdr) <- 'numeric'

    col.power <- t(sapply(1:nrow(mat.power),FUN = function(i){
      ifelse(mat.power[i,] == max(mat.power[i,which(mat.fdr[i,] < color.fdr)]) &
               mat.fdr[i,] < color.fdr,'blue','black')
    }))
    col.power[is.na(col.power)] <- 'black'
    col.fdr <- t(apply(mat.fdr,1,function(x){ifelse(x > color.fdr,'red','black')}))
    col.fdr[is.na(col.fdr)] <- 'black'

    for(imethod in methods){
      dt.dcast.color[[paste0('A.Power_',imethod)]] <-
        cell_spec(x = dt.dcast.color[[paste0('A.Power_',imethod)]],color = col.power[,paste0('A.Power_',imethod)],format = format)
      dt.dcast.color[[paste0('A.Power_',imethod)]] <- gsub("NA","-",dt.dcast.color[[paste0('A.Power_',imethod)]])

      dt.dcast.color[[paste0('B.FDR_',imethod)]] <-
        cell_spec(x = dt.dcast.color[[paste0('B.FDR_',imethod)]],color = col.fdr[,paste0('B.FDR_',imethod)],format = format)
      dt.dcast.color[[paste0('B.FDR_',imethod)]] <- gsub("NA","-",dt.dcast.color[[paste0('B.FDR_',imethod)]])
    }
  }

  kb <- kbl(dt.dcast.color,
            escape = FALSE,
            format = format,
            booktabs = TRUE,
            align = c('l',rep('r',13)),
            caption = cap,
            col.names = c('Read','Samples/Group','Library Size','Read Length',rep(methods,2)),...) %>%
    add_header_above(c(" " = 4, "Power" = length(methods), "False Discovery Rate" = length(methods))) %>%
    { if(format == 'latex'){
      dotDotDot <- match.call(expand.dots = FALSE)
      if('longtable' %in% names(dotDotDot$...)){
        getLongTable <- get('longtable')
      } else{
        getLongTable <- FALSE
      }
      if(isTRUE(getLongTable)){
        opts <- c("repeat_header")
      } else{
        opts <- c("scale_down")
      }
      kable_styling(kable_input = .,latex_options = opts,font_size = font_size) %>%
        collapse_rows(kable_input = .,columns = 1, latex_hline = "major", valign = "top")
    } else{
      kable_styling(kable_input = .,bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
        collapse_rows(kable_input = .,columns = 1, valign = "top") %>%
        landscape()
    }
    }

  return(kb)
}
