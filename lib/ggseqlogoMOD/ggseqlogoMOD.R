####
## Sequence logo generator for modified sequences
## Written by Marthe Solleder. Build upon ggseqlogo by Omar Wagih.
###


#### function to plot sequence logo
ggseqlogoMOD <- function( data,
                          libPath = NULL,
                          smallSampleCorr = FALSE,
                          col_scheme = 'phosphorylated',
                          additionalAA = 'sty',
                          seq_type = 'aa',
                          font = 'helvetica_phosphorylated',
                          legendText = FALSE,
                          ylim = c(0,log2(23)),
                          title = NULL,
                          titleSize = 24,
                          titlePos = 0.5,
                          axisTextSizeX = 18,
                          axisTextSizeY = 18,
                          axisTitleSize = 18,
                          ...){
  ## find length of peptides
  if(typeof(data) == 'character'){
    lengthP = nchar(data[[1]])
  } else if(typeof(data) == 'double'){
    lengthP = ncol(data)
  }

  #### plot sequence logo
  p = ggplot() +

    #### do the plotting of the sequence logo (by Omar Wagih's ggseqlogo)
    geom_logo(data=data,
              smallSampleCorr=smallSampleCorr,
              col_scheme=col_scheme,
              additionalAA=additionalAA,
              font=font,
              legendText=legendText, 
              libPath=libPath,
              ...) +

    #### size of title, size of x and y tick marks, size of y axis description
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.5),
          axis.ticks = element_line(colour = 'black', size = 0.5),
          plot.title = element_text(size=titleSize, hjust=titlePos, vjust = 0.1, family="sans"),
          axis.text.x = element_text(size=axisTextSizeX, color = 'black', family="sans"),
          axis.text.y = element_text(size=axisTextSizeY, color = 'black', family="sans"),
          axis.title.y = element_text(size=axisTitleSize, family="sans")) +

    #### remove space between axis and plot space, define x ticks
    scale_y_continuous(expand = c(0, 0)) +

    #### costum y axis range
    coord_cartesian(ylim = ylim) +
    scale_x_continuous(breaks = seq(1, lengthP, 1)) +

    #### add title if given
    ggtitle(title)

  return(p)
}
