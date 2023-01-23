library(edgeR)
library(data.table)
devtools::load_all('../pkg')

path <- '../../output/quasi_poisson'
out.path <- '../../misc/quasi_poisson'

# Annotation

fa <- scanFasta('../../data/annotation/mm39/gencode.vM27.transcripts.fa.gz')

anno <- strsplit2(fa$TranscriptID,"\\|")

fa$TranscriptID <- anno[,1]
fa$GeneID <- anno[,2]

df.gene <- data.table(fa)[,.N,by = 'GeneID']
df.gene[, color := ifelse(N == 1, 'red', 'black')]

# Data setup

catch <- catchKallisto(list.dirs(file.path(path,"tech/quant-kallisto"),recursive = FALSE))
anno.catch <- strsplit2(rownames(catch$annotation),"\\|")
catch$annotation$TranscriptID <- anno[,1]
catch$annotation$GeneID <- anno[,2]
catch$annotation$N <- df.gene$N[match(catch$annotation$GeneID,df.gene$GeneID)]
catch$annotation$color <- df.gene$color[match(catch$annotation$GeneID,df.gene$GeneID)]

meta <- fread(file.path(path,'tech/meta/counts.tsv.gz'))
meta$N <- df.gene$N[match(meta$GeneID,df.gene$GeneID)]
meta$color <- df.gene$color[match(meta$GeneID,df.gene$GeneID)]

targets <- data.table(strsplit2(basename(colnames(catch$counts)),"_"))
colnames(targets) <- c('group','replicate','r')
targets[,sample := paste0(group,'-',replicate)]

# True counts

dge.true <- DGEList(counts = as.matrix(meta[,grepl('group',colnames(meta)),with = FALSE]),
                    samples = targets,
                    genes = meta[,c('TranscriptID','GeneID','Length','N','color')])

keep.true <- filterByExpr(dge.true)
dge.true <- dge.true[keep.true, , keep.lib.sizes = FALSE]

rmeans.true <- rowMeans(dge.true$counts)
rvars.true <- matrixStats::rowVars(dge.true$counts)

is.singletx.true <- dge.true$genes$N == 1

# kallisto raw counts

dge <- DGEList(counts = catch$counts,
               samples = targets,
               genes = catch$annotation)

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)

rmeans <- rowMeans(dge$counts)
rvars <- matrixStats::rowVars(dge$counts)
sqvars <- squeezeVar(rvars,df = ncol(dge$counts) - 1,robust = TRUE)$var.post

is.singletx <- dge$genes$N == 1

# Plots

agg_png(filename = file.path(out.path,"suppfigure_quasipoisson.png"),width = 7.5/2,height = 7.5,units = 'in',res = 300)

par(mfrow = c(2,1))

plot(log10(rvars.true[!is.singletx.true])~log10(rmeans.true[!is.singletx.true]),xlim = c(0,7),ylim = c(0,7),
     col = 'black',
     pch = 20,
     cex = 0.01,
     ylab = "Transcript-level variance (log10 scale)",
     xlab = "Mean transcript expression (log10 scale)",
     main = 'True (simulated) counts',
     cex.main = 8/12,
     cex.lab = 8/12,
     cex.axis = 8/12)
points(log10(rvars.true[is.singletx.true])~log10(rmeans.true[is.singletx.true]),xlim = c(0,7),ylim = c(0,7),
       col = 'red',
       pch = 20,
       cex = 0.5)

legend('topleft',
       legend = c('Transcripts from multi-transcript genes',
                  'Transcripts from single-transcript-genes'),
       pch = 16,col = c('black','red'),cex = 6/12)

abline(a = 0,b = 1,col = 'red',lty = 2)

plot(log10(rvars[!is.singletx])~log10(rmeans[!is.singletx]),xlim = c(0,7),ylim = c(0,7),
     col = 'black',
     pch = 20,
     cex = 0.01,
     ylab = "Transcript-level variance (log10 scale)",
     xlab = "Mean transcript expression (log10 scale)",
     main = 'Estimated (quantified) counts',
     cex.main = 8/12,
     cex.lab = 8/12,
     cex.axis = 8/12)
points(log10(rvars[is.singletx])~log10(rmeans[is.singletx]),xlim = c(0,7),ylim = c(0,7),
       col = 'red',
       pch = 20,
       cex = 0.5)

abline(a = 0,b = 1,col = 'red',lty = 2)

dev.off()
