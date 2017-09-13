### Imports
.libPaths(c("/n/home11/adarnell/R/x86_64-unknown-linux-gnu-library/3.3/", .libPaths()))
library(tidyverse)  # tables, plotting, read-write
library(stringr)  # string manipulation
library(BSgenome)  # genome manipulation
library(magrittr)  # nice piping operators
library(Cairo)  # plot printing library

# color palette for codons
cbPalette <- c("#CC79A7", "#E69F00", "#D55E00",
               "#0072B2", "#56B4E9", "#009E73")

# function to read codon density and convert to single column
getdata <- function (file) {
  # read data
  data  <- read_tsv(file)
  # gather all columns except codon
  longdata  <- (
    gather(data, 'position', 'density', -codon)
    %>% mutate(position=as.integer(position))
  )
  longdata
}

# function for pairing data for later ribo density subtraction
pairsamples <- function(data, samplepair) {
  # get all points for sample1
  df1 <- data %>% filter(sample == samplepair[['bottom']])
  # get all points for sample2
  df2 <- data %>% filter(sample == samplepair[['top']])
  # join both data frames
  df <- (
    df1
    %>% inner_join(df2, by=c("codon", 'position', 'aa'))
    # choose only aa that will be plotted for this sample pair
    %>% filter(aa==samplepair[['aa']])
  )
  df
}

# get codon-aa pairings
codon.table  <- (
  GENETIC_CODE  
  %>% tbl_df()  
  %>% rownames_to_column(var='codon') 
  %>% dplyr::rename(., aa=value)
)

# get all codon ribo density files
codonfiles  <- list.files('../processeddata/',
                          pattern='codon.ribo.density.tsv',
                          recursive=TRUE,
                          full.names=TRUE
                         )
# parse sample names
samples  <- str_match(codonfiles, 'processeddata//(.*)/')[,2]
# get sample pairs for ribosome density diff calc
# these pairs are written in the '.org' file and exported to '.tsv'
samplepairs = read_tsv('../tables/samplepairs_codon_density.tsv')

# read all data
data <- lapply(codonfiles, getdata)
names(data) <- samples
# join all data and add 'aa' for each codon
(data
  %<>% bind_rows(.id='sample')
  %>% left_join(codon.table, by='codon')
)
# parse new columns containing cell line, starv, modification, drug
# for creating plot titles later
(data
  %<>% mutate(
         cellline=str_extract(sample, '293t|hela|hct116'),
         starv=str_extract(sample, 'arg|leu|rich'),
         modification=str_extract(sample, 'gcn2ko|hrgfp|ragbq99l'),
         drug=str_extract(sample, 'torin')
         )
)

################### cell type plot ##############################

# calculate density difference and rename sample pair for plot title
subtractsamples <- function(data) {
  (data
    %<>% mutate(
           density.diff=density.y-density.x,
           samplepairs=str_c(starv.y, ' | ', cellline.y)
         )
    # discard columns having .x and .y at the end
    %>% select(-matches(".*x$|.*y$"))
  )
  data
}

# select only pairs for plot 1
plotpairs <- samplepairs %>% filter(plotnumber == 1)
# get difference for all samplepairs
diffdata <- bind_rows(lapply(
  1:nrow(plotpairs),
  function(n) subtractsamples(pairsamples(data, plotpairs[n,]))))

#### plot 1
cairo_pdf('../figures/codondensitydifference_celltypes.pdf', width=6, height=4)
p <- ggplot(diffdata, aes(x=position, y=density.diff, color=codon))
(p
  + geom_line()
  + theme_classic()
  + theme(
      strip.background=element_blank(),
      axis.text=element_text(size=6),
      axis.title=element_text(size=8),
      legend.text=element_text(size=6),
      legend.title=element_text(size=8),
      legend.key.size = unit(0.6, 'lines')
    )
  + labs(x='Distance from codon (nt)', y='Ribosome Density\n(Starved \u2013 Rich)')
  + scale_x_continuous(limits=c(-100, 100))
  + scale_y_continuous(limits=c(-0.5, 1.5))
  + scale_colour_manual(values=rep(cbPalette, 2))
  # change 'fixed' to 'free' if we want separate axis for each panel
  + facet_wrap(~samplepairs, ncol=3, scales='fixed')
 + guides(color=guide_legend(title="Codon"))
)
dev.off()

###################### hrgfp, ragbq99l, gcn2ko plots #####################
# calculate density difference and rename sample pair for plot title
subtractsamples <- function(data) {
  (data
    %<>% mutate(
           density.diff=density.y-density.x,
           samplepairs=str_c(starv.y, ' | ', modification.y)
         )
    # discard columns having .x and .y at the end
    %>% select(-matches(".*x$|.*y$"))
  )
  data
}

# select only pairs for plot 2
plotpairs <- samplepairs %>% filter(plotnumber == 2)
# get difference for all samplepairs
diffdata <- bind_rows(lapply(
  1:nrow(plotpairs),
  function(n) subtractsamples(pairsamples(data, plotpairs[n,]))))

# change plot order
diffdata$samplepairs = factor(diffdata$samplepairs,
                              levels=c(
                                "arg | hrgfp",
                                "arg | ragbq99l",
                                "arg | gcn2ko",
                                "leu | hrgfp",
                                "leu | ragbq99l",
                                "leu | gcn2ko"))
#### plot 2
cairo_pdf('../figures/codondensitydifference_gcn2_tor.pdf', width=6, height=4)
p <- ggplot(diffdata, aes(x=position, y=density.diff, color=codon))
(p
  + geom_line()
  + theme_classic()
  + theme(
      strip.background=element_blank(),
      axis.text=element_text(size=6),
      axis.title=element_text(size=8),
      legend.text=element_text(size=6),
      legend.title=element_text(size=8),
      legend.key.size = unit(0.6, 'lines')
    )
  + labs(x='Distance from codon (nt)', y='Ribosome Density\n(Starved \u2013 Rich)')
  + scale_x_continuous(limits=c(-100, 100))
  + scale_y_continuous(limits=c(-0.25, 3.0))
  + scale_colour_manual(values=rep(cbPalette, 2))
  # change 'fixed' to 'free' if we want separate axis for each panel
  + facet_wrap(~samplepairs, ncol=3, scales='fixed')
 + guides(color=guide_legend(title="Codon"))
)
dev.off()

###################### torin plots #####################
# calculate density difference and rename sample pair for plot title
subtractsamples <- function(data) {
  (data
    %<>% mutate(
           density.diff=density.y-density.x,
           samplepairs=str_c(starv.y, ' | ', drug.y)
         )
    # discard columns having .x and .y at the end
    %>% select(-matches(".*x$|.*y$"))
  )
  data
}

# select only pairs for plot 3
plotpairs <- samplepairs %>% filter(plotnumber == 3)
# get difference for all samplepairs
diffdata <- bind_rows(lapply(
  1:nrow(plotpairs),
  function(n) subtractsamples(pairsamples(data, plotpairs[n,]))))

#### plot 3
cairo_pdf('../figures/codondensitydifference_torin.pdf', width=4, height=2)
p <- ggplot(diffdata, aes(x=position, y=density.diff, color=codon))
(p
  + geom_line()
  + theme_classic()
  + theme(
      strip.background=element_blank(),
      axis.text=element_text(size=6),
      axis.title=element_text(size=8),
      legend.text=element_text(size=6),
      legend.title=element_text(size=8),
      legend.key.size = unit(0.6, 'lines')
    )
  + labs(x='Distance from codon (nt)', y='Ribosome Density\n(Starved \u2013 Rich)')
  + scale_x_continuous(limits=c(-100, 100))
  + scale_y_continuous(limits=c(-0.5, 1.5))
  + scale_colour_manual(values=rep(cbPalette, 2))
  # change 'fixed' to 'free' if we want separate axis for each panel
  + facet_wrap(~samplepairs, ncol=3, scales='fixed')
 + guides(color=guide_legend(title="Codon"))
)
dev.off()
