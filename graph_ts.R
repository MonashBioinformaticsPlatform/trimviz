#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(ape)
library(reshape2)
library(gridExtra)


## Collect arguments
args <- commandArgs(TRUE)
if (length(args) != 4){
  print (args)
  stop('incorrect number of arguments given to graph_ts.R')
}
#setwd(args[3])
aggFN = args[3]
maxAggN = as.integer(args[4])
df <- read.table(file = args[1], comment.char = '', header=TRUE)
#df <- read.table(file = './trimVisTmpFiles/trimviz_readData.tsv', comment.char = '', header=TRUE)
graphchunk = 5
rnames = unique(df$read)
maxLen =  max( df$position )
for (col in (7:ncol(df))){
  df[df[,col]==0,col] = NA  # set 0 to NA to colour black
}

# choose first adapter for graph (col 7)
df['consec_adapt_residues'] = as.numeric(df[,7]) 
df$seq=as.character(df$seq)
genSeq=F
if ('genomic_seq' %in% colnames(df)){
  df$genomic_seq = as.character(df$genomic_seq)
  genSeq=T
  
}
print('current directory:')
print(getwd())
#print(args[2])
pdf(args[2], width = 14, height = 8)
maxcol=max(df$consec_adapt_residues, na.rm = T)
for (i in 1:as.integer(length(rnames)/graphchunk)){
  tograph = rnames[(i*graphchunk-1):(i*graphchunk+graphchunk-2)]
  temp = df[df$read %in% tograph,]
  temp=temp[order(match(temp$read,df$df), temp$position),]
  temp$read = substr(temp$read, 21, nchar(as.character(temp$read)))
  temp$read=factor(temp$read, levels = unique(temp$read))
  gr <- ggplot(data=temp, aes(x=position, y=qual, label=seq)) +
    geom_rect(data= temp, aes(xmax = tp_cutoff-0.5, xmin = fp_cutoff-0.5, ymin = -12, ymax = max(qual)+10), size=0.01, colour = 'white', fill = 'white') +
    geom_hline(data = data.frame(yint=c(0:5)*10), aes( yintercept = yint), colour = "#DDDDDD", size=0.5) +
    geom_line() + geom_point() + 
    geom_text(data=temp, mapping=aes(x=position, y=-4, label=seq, colour=consec_adapt_residues), size=2.7, fontface="bold") +
    scale_colour_gradient(low="blue", high="red", na.value = "black") +
    geom_vline(aes(xintercept = fp_cutoff-0.5), col="red") +
    geom_vline(aes(xintercept = tp_cutoff-0.5), col="blue") +
    facet_grid(read ~ .) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.background = element_rect(fill = '#F0F0F0', colour = '#F0F0F0'),
          legend.title = element_text(size = 10))+
    ggtitle("grey zone = trimmed")+
    labs(title = "Random sample of reads in trimmed read file", x = "Position on read. Grey zones indicate trimmed bases.", y = 'Quality of base call')
  if(genSeq){
  gr <- gr +  
    geom_text(data=temp, mapping=aes(x=position, y=-10.5, label=genomic_seq), col='black', size=2.7) +
    geom_text(data=temp, mapping=aes(x=position, y=-8.3, label=ifelse(genomic_seq==seq,' ', '*')), size=3.5, fontface="bold", col='red') + 
    coord_cartesian(ylim=c(-10.8, max(temp$qual)+5), xlim=c(0,maxLen))
  } else {
    gr <- gr + coord_cartesian(ylim=c(-8, max(temp$qual)+5), xlim=c(0,maxLen))
  }
  print(gr)
}
dev.off()


df2 <- read.table('./trimVisTmpFiles/seq3psites.txt', header = T)
if (nrow(df2) == 0){
  stop("No reads in trimmed-class dataset for making aggregate trim plots. Perhaps there are no trimmed reads in the data?")
}
subSampN=min(maxAggN, nrow(df2))
df2 = df2[sample(1:subSampN, replace = F),]

blocks = list()
lblocks = list()
reord = list()
first=T
absentData=list()
for (prefix in c('s','q','g')){
  scols=grepl(pattern= paste0('^',prefix,'[\\d]+$'),colnames(df2), perl=T)
  if(sum(scols) > 0){
    smat = as.matrix(df2[,scols])
    rownames(smat) = df2$readID
    aggFlank=(sum(scols)-1)/2
    if (prefix %in% c('s','g')){ # seq
      #smat[smat=='X']='N'
      d = ape::dist.gene(smat)
    }
    else { # numerical val
      d = dist(smat)
    }
    hc = hclust(d)
    slong = melt(data = smat)
    colnames(slong)=c('readId', 'pos', prefix)
    slong$pos=as.numeric(substr(slong$pos, 2, 9999)) - 1 - aggFlank
    sreord = data.frame(readId = hc$labels[hc$order], order = 1:length(hc$labels))
    colnames(sreord)[2] = paste0('ord.',prefix)
    blocks[[prefix]] = smat
    lblocks[[prefix]] = slong
    reord[[prefix]] = sreord
    sresults = merge(slong, sreord, by = 'readId')
    if (first){
      allres = sresults
      first=F
    } else {
      allres=merge(allres, sresults, by=c('readId','pos'))
      allres=allres[order(allres$pos),]
    }
  } else {
    absentData = c(absentData, prefix)
  }
}


anchorPlot <- function(x, datatype='s', clustby='s'){
  ordCol=paste0('ord.',clustby)
  gr = ggplot(x) + geom_raster(aes_string(x = 'pos', y = ordCol, fill = datatype)) + theme_bw() + 
  #gr = ggplot(x, aes(x = pos,y = ord.s)) + geom_raster(aes(fill = s)) + theme_bw() + 
    geom_vline(xintercept = 0, linetype='dashed', col='black')
  if (is.numeric(x[[datatype]])){
    gr = gr + scale_fill_gradientn(colours = heat.colors(80))
  } else {
    gr = gr + scale_fill_manual(values = c("orange","green","blue","888888","red", "black"))  # A C G N T X
  }
  return(gr)
}

#print(aggFN)
pdf(aggFN, width=15, height=50)
if (length(absentData) > 0 && absentData == 'g'){
  grid.arrange(anchorPlot(allres, 's','s'), anchorPlot(allres, 'q','s'), 
               anchorPlot(allres, 's','q'), anchorPlot(allres, 'q','q'), 
               ncol=2)
} else {
grid.arrange(anchorPlot(allres, 's','s'), anchorPlot(allres, 'q','s'), anchorPlot(allres, 'g','s'), 
             anchorPlot(allres, 's','q'), anchorPlot(allres, 'q','q'), anchorPlot(allres, 'g','q'),
             anchorPlot(allres, 's','g'), anchorPlot(allres, 'q','g'), anchorPlot(allres, 'g','g'),
             ncol=3)
}
dev.off()
