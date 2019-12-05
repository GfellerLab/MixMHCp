#!/usr/bin/env Rscript

####
## Calling ggseqlogoMOD to build sequence logos.
## Written by Marthe Solleder.
####


## packages
library(ggplot2)

## input arguments
args <- commandArgs(trailingOnly = TRUE)
pathLib <- args[1]
input <- args[2]
output <- args[3]
alphabet <- args[4]
clusters_min <- args[5]
clusters_max <- args[6]
inputType <- args[7]




## load functions
source(paste(pathLib, 'R/ggseqlogo.r', sep = ''))
source(paste(pathLib, 'R/col_schemes.r', sep = ''))
source(paste(pathLib, 'R/heights.r', sep = ''))

pdf(NULL)

## check additional residues in alphabet
addAlph <- c()
alph <- c('A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','V','W','Y')
A <- as.list(strsplit(alphabet, '')[[1]])
for(a in A){
  if(!(a %in% alph)){
    addAlph <- append(addAlph, a)
  }
}
addAlph <- paste(addAlph, collapse = '')

## define ylimit for plots by ggseqlogoMOD depending on alphabet size
if(length(A) <= 20){
  ylimit <- c(0, log2(20))
} else{
  ylimit <- c(0, log2( length(A) ))
}

## function to build PWM
buildPWM <- function(peptides, alphabet, lengthP){
  
  L <- length(peptides)
  ## build empty PWM
  pwm <- matrix(0, nrow=lengthP, ncol=length(alphabet))
  colnames(pwm) <- alphabet
  ## fill PWM
  for(s in 1:L){
    u <- strsplit(peptides[s], split = '')
    for(p in 1:lengthP){
      pwm[p,u[[1]][p]] <- pwm[p,u[[1]][p]] + 1
    }
  }

  ## return PWM
  return(pwm)
}

if (inputType == 'pwm'){
  
  print('Input type read from MixMHCp: PWM')
  
  ## read all pwm of directory into dataframe
  pathInput <- paste0(input, 'ggseqlogo')
  files <- list.files(path=pathInput, pattern="*.txt", full.names=T, recursive=FALSE)
  df <- lapply(files, read.table, row.names = 1)
  N <- names(df) <- gsub('.*pwm_', '', files)
  
} else {
  
  print('Input type read from MixMHCp: responsibility file')
  
  ## go over all clusters for which motif deconvolution was performed
  pathInput <- paste0(input, 'responsibility/')
  listPeptides <- list()
  df <- list()
  N <- c()
  lengthsOfData <- c()
  
  ## go over all clusters defined at input
  for(i in clusters_min:clusters_max){
    
    ## open responsibilty file
    d <- read.table(paste0(pathInput, 'resp_', i, '.txt'),
                    header = TRUE,
                    sep = "\t")
    
    ## depending on amount of clusters (1, ... m), for each peptide, select sub DF
    respValues <- d[,c(2:(i+2))]
    
    ## check for all rows which resp. value is the highest
    cols <- colnames(respValues)[apply(respValues,1,which.max)]
    
    ## iterate over all peptides (length/rows of df)
    for(r in 1:length(cols)){
      p <- as.character(d[r,1])
      lengthsOfData <- append(lengthsOfData, nchar(p))
      
      ## decide on name (number of lcuster or trash)
      if(cols[r] == 'Trash'){
        k <- paste0('L', nchar(p), '_', i, '_Trash' )
      }else{
        n <- substr( cols[r], 2, nchar(cols[r]) )
        k <- paste0('L', nchar(p), '_', i, '_', n )
      }
      N <- append(N, k)
      listPeptides[[k]] <- append(listPeptides[[k]], p)
    }
    
    N <- N[!duplicated(N)]
    lengthsOfData <- lengthsOfData[!duplicated(lengthsOfData)]
  }
    
  ## check if for all lengths, every cluster c has 1:c subclusters
  for(peptL in lengthsOfData){
    #for(c1 in 1:args$clusters){
     for(c1 in clusters_min:clusters_max){
        if( !(exists(paste0('L', peptL, '_', c1, '_Trash'), listPeptides)) ){
        N <- append( N, paste0('L', peptL, '_', c1, '_Trash') )
      }
      for(c2 in 1:c1){
        if( !(exists(paste0('L', peptL, '_', c1, '_', c2), listPeptides)) ){
          N <- append( N, paste0('L', peptL, '_', c1, '_', c2) )
        }
      }
    }
  }

  ## build pwm
  for(c in N){
    ## length of peptides
    l <- as.numeric(substr(strsplit(c, '_')[[1]][1], 2, nchar(strsplit(c, '_')[[1]][1])))
    if(exists(c, where = listPeptides)){
      ## build pwm of peptide list
      pwm <- t(buildPWM(listPeptides[[c]], A, l))
      df[[c]] <- pwm
    } else{
      ## build pwm of zeros
      pwm <- matrix(0L, nrow = length(A), ncol = l)
      df[[c]] <- pwm
    }
  }

}


## go through each PWM in dataframe and plot sequence logo
for(w in 1:length(df)){

  pwm <- df[[w]]
  
  ## check if matrix is 0
  if( all(pwm == 0) == FALSE ){

    ## column sum for amount of peptides
    len <- colSums(pwm)[1]
    
    ## normalize pwm
    pwm <- pwm/colSums(pwm)
    pwm <- as.matrix(pwm)
    
    ## build sequence logo
    l <- ggseqlogoMOD(data = pwm,
                      libPath = pathLib,
                      additionalAA = addAlph,
                      ylim = ylimit,
                      axisTextSizeX = 6,
                      axisTextSizeY = 6,
                      axisTitleSize = 6) +
      theme(plot.margin = margin(0.05, 0.01, -0.45, -0.29, "cm"),
            axis.line = element_line(colour = "black", size = 0.2),
            axis.ticks = element_line(colour = 'black', size = 0.2))
    
    ## save plot
    f <- paste0('logos/logo_', gsub('.txt', '',N[w]), '-', as.integer(len))
    print(paste0('logos/logo_', gsub('.txt', '',N[w]), '-', as.integer(len)))
    ggsave(paste0(f, '.png'), 
           device = 'png',
           path = output, 
           width = 0.6 * (1 + dim(pwm)[2]), 
           height = 4, 
           units = 'cm')

  } else{
    l <- ggplot()
    f <- paste0('logos/logo_', gsub('.txt', '', N[w]), '-', 0)
    ggsave(paste0(f, '.png'), 
           device = 'png',
           path = output, 
           width = 0.6 * (1 + dim(pwm)[2]), 
           height = 4, 
           units = 'cm')
  }

}



####
## END
####
