#' @importFrom wasabi prepare_fish_for_sleuth
#' @importFrom jsonlite fromJSON
runSalmon <- function(bin,index,targets,dest,
                      options = '-l A --numBootstraps 100 --validateMappings -p 6'){
  
  dir.create(dest,showWarnings = FALSE,recursive = TRUE)
  dest <- normalizePath(dest)
  
  paired.end <- ifelse(ncol(targets) == 2,TRUE,FALSE)
  
  for (lib.num in seq_len(nrow(targets))) {
    
    targets.file <- as.character(targets[lib.num,])
    sample.name <- gsub('.fastq.gz','',basename(targets.file[1]))
    dest.sample.name <- file.path(dest,sample.name)
    
    # Check if Salmon needs to be run
    if (file.exists(file.path(dest.sample.name, 'aux_info/meta_info.json')) &
        file.exists(file.path(dest.sample.name, 'aux_info/bootstrap/bootstraps.gz'))) {
      js <- fromJSON(file.path(dest.sample.name, 'aux_info/meta_info.json'))
      boot <- gzcon(file(file.path(dest.sample.name, 'aux_info/bootstrap/bootstraps.gz'), open = "rb"))
      
      nboot.tx <- js$num_valid_targets
      nboot.json <- js$num_bootstraps
      nboot.actual <- readBin(boot, what = "double", n = nboot.tx * nboot.json)
      nboot.tx*nboot.json == length(nboot.actual)
      
      run <- nboot.tx*nboot.json != length(nboot.actual)
    } else{
      run <- TRUE
    }
    
    if (isTRUE(run)) {
      if (isFALSE(paired.end)) {
        cmd.sample <- paste('-r',targets.file)
      } else {
        cmd.sample <- paste(c('-1','-2'),targets.file,collapse = ' ')
      }
      cmd <- paste('quant',
                   '-i',index,
                   options,
                   cmd.sample,
                   '-o',dest.sample.name)
      system2(command = bin,args = cmd)
    }
    
    # Preparing salmon for sleuth
    if (!file.exists(file.path(dest.sample.name, 'abundance.h5'))) {
      prepare_fish_for_sleuth(dest.sample.name)
    }
  }
}

#' @importFrom rhdf5 h5ls
runKallisto <- function(bin,index,targets,dest,
                        options = '--bootstrap-samples=100 --threads=6'){
  
  dir.create(dest,showWarnings = FALSE,recursive = TRUE)
  dest <- normalizePath(dest)
  
  paired.end <- ifelse(ncol(targets) == 2,TRUE,FALSE)
  
  for (lib.num in seq_len(nrow(targets))) {
    
    targets.file <- as.character(targets[lib.num,])
    sample.name <- gsub('.fastq.gz','',basename(targets.file[1]))
    dest.sample.name <- file.path(dest,sample.name)
    
    # It is not possible to know if kallisto has been run to completion without
    # manual inspection of the abundance.h5 file
    if (file.exists(file.path(dest.sample.name, 'run_info.json'))) {
      nboot <- as.numeric(strsplit2(options,"--bootstrap-samples=| ")[1,2])
      run <- tryCatch(!sum(grepl('\\/bootstrap',h5ls(file.path(dest.sample.name,'abundance.h5'))[,'group'])) == nboot,error = function(e) TRUE)
    } else{
      run <- TRUE
    }
    
    if (isTRUE(run)) {
      if (isFALSE(paired.end)) {
        cmd.sample <- paste('-l 180 -s 40 --single',targets.file)
      } else {
        cmd.sample <- paste(targets.file,collapse = ' ')
      }
      cmd <- paste('quant',
                   '-i',index,
                   options,
                   '-o',dest.sample.name,
                   cmd.sample)
      system2(command = bin,args = cmd)
    }
  }
}