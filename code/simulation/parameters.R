OK <- requireNamespace("data.table", quietly = TRUE)
if (!OK) {
  stop("data.table package required but is not installed (or can't be loaded)")
}
library(data.table)

dest.data <- '../../output/simulation/data'

for (genome in c('mm39')) {
  for (len in c(50, 75, 100, 125, 150)) {
    for (fc in c(1, 2)) {
      for (scenario in c('balanced', 'unbalanced')) {
        for (paired.end in c(FALSE, TRUE)) {
          for (max.tx in c(2, 3, 4, 5, 9999)) {
            for (libs.per.group in c(3, 5, 10)) {
              for (simulation in 1:20) {
                dest <- file.path(dest.data,
                                  genome,
                                  paste0('readlen-',len),
                                  paste0('fc',gsub("\\.","_",as.character(fc))),
                                  paste0(ifelse(paired.end,'paired','single'),'-end'),
                                  paste0(max.tx,'TxPerGene'),
                                  scenario,
                                  paste0(libs.per.group,'libsPerGroup'),
                                  paste0('simulation-',simulation))
                dir.create(dest,recursive = TRUE,showWarnings = FALSE)

                if (!length(list.files(dest, 'time.tsv', recursive = TRUE)) == 2) {
                  dt <- data.table(dest = dest,
                                   genome = genome,
                                   rlen = len,
                                   fc = fc,
                                   pe = paired.end,
                                   mtx = max.tx,
                                   scenario = scenario,
                                   libs = libs.per.group,
                                   simulation = simulation)

                  fwrite(x = dt,
                         file = 'parameters.txt',
                         sep = "\t",
                         quote = FALSE,
                         append = TRUE,
                         row.names = FALSE)
                }
              }
            }
          }
        }
      }
    }
  }
}
