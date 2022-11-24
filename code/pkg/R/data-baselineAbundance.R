#' Generator of baseline abundances based on GSE60450 data
#'
#' This object is a function returned from `stats::approxfun` which is used
#' in the simulation study to generate baseline abundances of gene-level expression.
#' The gene list used to create `baselineAbundance_function` is exported as a character
#' vector in the `baselineAbundance_genes` object.
#'
#' @docType data
#'
#' @usage 
#' data(baselineAbundance_function) 
#' data(baselineAbundance_genes)
#' 
#' @details 
#' Genes included in `baselineAbundance_genes` are protein-coding and lncRNA 
#' genes from the reference chromosomes 1, ..., 19, X, and Y from the mouse 
#' Ensembl version 106 annotation (GENCODE - Mus musculus - release M27). Only
#' genes with expected CPM>1 in at least 6 out of the 12 libraries are used.
#'
#' @format 
#' baselineAbundance_function is a function and baselineAbundance_genes is a
#' character vector
#'
#' @keywords datasets
#'
#' @references Fu NY, Rios AC, Pal B, Soetanto R et al. EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival. Nat Cell Biol 2015 Apr;17(4):365-75. PMID: 25730472
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450}{NCBI}
#'
#' @examples
#' \dontrun{
#' library(AnnotationHub) #. v3.4.0
#' library(tximeta) #. v1.14.1
#' library(SummarizedExperiment) #. v1.26.1
#' library(edgeR) #. v3.38.2
#' 
#' #. Loading annotation ("AH100674" corresponds to Ensembl version 106)
#' ah <- AnnotationHub()
#' edb <- ah[['AH100674']]
#' ensid <- keys(edb)
#' gene_ids <- select(edb,ensid,c("TXIDVERSION","GENEBIOTYPE","GENEIDVERSION","SEQNAME"))
#' 
#' #. Loading quantifications
#' files.quant <- list.files('../data/GSE60450/salmon/','quant.sf',recursive = TRUE,full.names = TRUE)
#' names(files.quant) <- basename(dirname(files.quant))
#' 
#' df.coldata <- data.frame(files = files.quant,names = names(files.quant))
#' se <- tximeta(coldata = df.coldata)
#' se.gene <- summarizeToGene(se)
#' 
#' #. Filtering protein coding and lncRNA genes from standard chromosomes
#' is.tx.relev <- gene_ids$GENEBIOTYPE %in% c('protein_coding','lncRNA') &
#'   gene_ids$SEQNAME %in% c(paste0(1:19),'X','Y')
#' gene.relev <- unique(gene_ids[is.tx.relev,'GENEIDVERSION'])
#' is.gene.relev <- rownames(se.gene) %in% gene.relev
#' se.gene.subset <- se.gene[is.gene.relev,]
#' 
#' #. Estimating baseline proportions (goodTuringProportions needs integers)
#' cts.gene <- round(assay(se.gene.subset))
#' prop <- goodTuringProportions(cts.gene)
#' 
#' #. Choose genes with expected CPM>1 in at least 6 libraries
#' is.expr <- rowSums(prop > 1e-6) >= 6
#' prop <- prop[is.expr,]
#' baselineAbundance_genes <- rownames(prop)
#' 
#' #. Creating interpolation function and saving object
#' i <- seq(from=0,to=1,length=length(prop))
#' baselineAbundance_function <- approxfun(i,sort(prop),rule=2)
#' 
#' save(baselineAbundance_genes,file = 'data/baselineAbundance_genes.rda',compress = 'xz')
#' save(baselineAbundance_function,file = 'data/baselineAbundance_function.rda',compress = 'xz')
#'  }
"baselineAbundance_function"
"baselineAbundance_genes"