#' Set of transcripts of interest from GSE60450
#'
#' Set of filtered transcripts of interest from GSE60450 from selected 
#' protein-coding genes of reference chromosomes from basic Gencode M27 mouse
#' annotation
#'
#' @docType data
#'
#' @usage data(GSE60450)
#'
#' @format data.frame
#'
#' @keywords datasets
#'
#' @references Fu NY, Rios AC, Pal B, Soetanto R et al. EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival. Nat Cell Biol 2015 Apr;17(4):365-75. PMID: 25730472
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450}{NCBI}
#'
#' @examples
#' \dontrun{
#' library(Rsubread) #. v2.10.4
#' library(limma) #. v3.52.2
#' library(data.table) #. v1.14.2
#' library(stringi) #. v1.7.8
#' library(devtools) #. v2.4.4
#' load_all()
#' 
#' #. Reading GSE60450 quantification
#' #. (see rfun::catchSalmon2, which brings in the TPM values from Salmon)
#' edger.catch <- catchSalmon2(list.dirs('../data/GSE60450/salmon/',recursive = FALSE))
#' colnames(edger.catch$counts) <- basename(colnames(edger.catch$counts))
#' colnames(edger.catch$tpm) <- basename(colnames(edger.catch$tpm))
#' 
#' #. Reading transcriptome and getting transcript duplicate information
#' contigs <- scanFasta('../../../genome/Gencode/M27/gencode.vM27.transcripts.fa.gz')
#' tx.info <- strsplit2(contigs$TranscriptID,"\\|")
#' contigs[,c('TranscriptID','GeneID','TranscriptType')] <- tx.info[,c(1,2,8)]
#' 
#' #. Selecting genes of interest (see ?baselineAbundance_function) and 
#' #. filtering their associated (unique) transcripts of interest. Genes of 
#' #. interest are protein-coding and lncRNA genes from reference chromosomes with 
#' #. expected CPM>1 in at least 6 out of the 12 libraries in GSE60450. Transcripts
#' #. of interest are those protein-coding and lncRNA transcripts from the selected 
#' #. genes.
#' data("baselineAbundance_genes")
#' is.relev.transcript <- 
#'   (contigs$Unique == TRUE & 
#'      contigs$TranscriptType %in% c('protein_coding','lncRNA') &
#'      contigs$GeneID %in% baselineAbundance_genes)
#' 
#' contigs.sub <- contigs[is.relev.transcript, ]
#' 
#' #. Ranking filtered transcripts based on Salmon's TPM
#' mat.tpm <- edger.catch$tpm[match(contigs.sub$TranscriptID,rownames(edger.catch$tpm)),]
#' 
#' dt.tpm <- data.table(TranscriptID = contigs.sub$TranscriptID,
#'                      GeneID = contigs.sub$GeneID,
#'                      TPM = mat.tpm)
#' 
#' dt.tpm$aveTPM <- exp(rowMeans(log1p(dt.tpm[,paste0('TPM.',colnames(edger.catch$tpm)),with = FALSE]))) - 1
#' dt.tpm$aveTPM <- 1e6*dt.tpm$aveTPM/sum(dt.tpm$aveTPM)
#' 
#' dt.tpm[order(-aveTPM),Rank := seq_len(.N),by = 'GeneID']
#' 
#' #. Bringing rankings back to the annotation and saving data
#' contigs.sub$Rank <- dt.tpm$Rank[match(contigs.sub$TranscriptID,dt.tpm$TranscriptID)]
#' GSE60450 <- contigs.sub[order(contigs.sub$GeneID,contigs.sub$Rank),]
#' rownames(GSE60450) <- NULL
#' save(GSE60450,file = 'data/GSE60450.rda',compress = 'xz')
#'  }
"GSE60450"