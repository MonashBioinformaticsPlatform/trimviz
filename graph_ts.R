#!/usr/bin/env Rscript
#library(dplyr)
library(ggplot2)
library(ape)
library(reshape2)
library(gridExtra)

anchorPlot <- function(x, datatype='s', clustby='s', EOI = '3p', gdiff=F){
  ordCol=paste0('ord.',clustby)
  ttls1 = list(s='Read sequence; ', g='Genomic sequence; ', q='Base-quality; ', t='Trim location; ')
  ttls2 = list(s='Clustered by read sequence', g='Clustered by genomic sequence', q='Clustered by base-quality patterns', t='Clustered by trim location')
  x$xalpha=F
  if (gdiff && (datatype=='g')){
    x$xalpha=(x$g==x$s)
    ttls2[[clustby]] = '(diff, c.f. read) '
  }
  gr = ggplot(x, aes_string(x = 'pos', y = ordCol, fill = datatype)) + 
    geom_raster(aes(alpha = as.numeric(pos == 0 | xalpha))) + # | s == 'N'))) + 
    theme_bw() + 
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"), legend.margin=margin(t = 0, unit='cm')) +
    scale_alpha(range = c(1,0)) + guides(alpha=FALSE)+
    geom_vline(xintercept = 0, linetype='dashed', col='black') + ggtitle(label = paste0(ttls1[[datatype]],ttls2[[clustby]])) +
    xlab(paste0('nucleotide position relative to ',EOI,'-trim site')) + ylab('')
  if (is.numeric(x[[datatype]])){
    minq = min (x[[datatype]][x$s %in% c('A','C','G','T')])
    maxq = max (x[[datatype]][x$s %in% c('A','C','G','T')]) 
    gr = gr + scale_fill_gradientn(colours = heat.colors(80), limits = c(minq, maxq)) +
    guides(fill = guide_colourbar(barwidth = 0.5, title = 'base-call quality', 
                                      title.position = "left", title.theme = element_text(angle = 90, size = 11)))
  } else {
    gr = gr + scale_fill_manual(values = c("orange","green","#4444FF","#707070","red", "black")) + 
    guides(fill = guide_legend( keywidth = 0.5, title = 'Sequence (X = insertion)',
                                    title.position = "left", title.theme = element_text(angle = 90, size = 11)))
    # A C G N T X
  }
  return(gr)
}



## Collect arguments
args <- commandArgs(TRUE)
if (length(args) != 3){
  print (args)
  stop('incorrect number of arguments given to graph_ts.R')
}
out_DN = args[1]
maxAggN = as.integer(args[2])
gdiff = (as.character(args[3])=='True')
reads_FN = paste0(out_DN, '/trimVisTmpFiles/trimviz_readData.tsv')
#seq3p_FN = paste0(out_DN, '/trimVisTmpFiles/seq3psites.txt')
#seq5p_FN = paste0(out_DN, '/trimVisTmpFiles/seq5psites.txt')
#aggFN = paste0(out_DN, '/TVheatmap.pdf')

############################
# part 1: indiv read plots #
############################

df <- read.table(file = reads_FN, comment.char = '', header=TRUE)
#df <- read.table(file = './trimVisTmpFiles/trimviz_readData.tsv', comment.char = '', header=TRUE)
graphchunk = 5
rnames = unique(df$read)
maxLen =  max( df$position )
# if all colnames are standard, it means there is no adapter column
ad_cols = ! (colnames(df) %in% c('read', 'position', 'seq', 'qual', 'fp_cutoff', 'tp_cutoff','trim_class','genomic_seq'))
if (sum(ad_cols) == 0){
  df$dummy_adapter=0
  ad_cols = c(ad_cols, T)
}

# choose first adapter for graph (col 7)
df['consec_adapt_residues'] = as.numeric(df[,which(ad_cols)[1]]) 
maxcol=max(df$consec_adapt_residues, na.rm = T)
df$consec_adapt_residues[df$consec_adapt_residues==0] = NA  # set 0 to NA to colour black

maxQual = max(df$qual)
df$seq=as.character(df$seq)
genSeq=F
if ('genomic_seq' %in% colnames(df)){
  df$genomic_seq = as.character(df$genomic_seq)
  genSeq=T
}
pdf(paste0(out_DN, '/indiv_reads.pdf'), width = 14, height = 8)
for (i in 0:floor(length(rnames)/graphchunk)){
  tograph = rnames[((i*graphchunk)+1):(i*graphchunk+graphchunk)]
  temp = df[df$read %in% tograph,]
  if (nrow(temp) > 0){
    temp=temp[order(match(temp$read,df$df), temp$position),]
    temp$read = substr(temp$read, 21, nchar(as.character(temp$read)))
    temp$read=factor(temp$read, levels = unique(temp$read))
    gr <- ggplot(data=temp, aes(x=position, y=qual, label=seq)) +
      geom_rect(data= temp, aes(xmax = tp_cutoff-0.5, xmin = fp_cutoff-0.5, ymin = -12, ymax = max(qual)+10), size=0.01, colour = 'white', fill = 'white') +
      geom_hline(data = data.frame(yint=c(0:5)*10), aes( yintercept = yint), colour = "#DDDDDD", size=0.5) +
      geom_line() + geom_point() + 
      scale_y_continuous(breaks = seq(0,40, by=10)) +
      geom_text(data=temp, mapping=aes(x=position, y=-4, label=seq, colour=consec_adapt_residues), size=2.7, fontface="bold") +
      scale_colour_gradient(low="blue", high="red", na.value = "black", limits=c(5, maxcol)) +
      geom_vline(aes(xintercept = fp_cutoff-0.5), col="red") +
      geom_vline(aes(xintercept = tp_cutoff-0.5), col="blue") +
      guides(colour = guide_colourbar(barwidth = 0.5, title = 'conseq. adapter bases', 
                                      title.position = "left", title.theme = element_text(angle = 90))) +
      facet_grid(read + trim_class ~ .) +
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
      coord_cartesian(ylim=c(-10.8, max(temp$qual)+5), xlim=c(0,maxLen)) +
      geom_text(data = data.frame(lab = c('read:','ref:'), y=c(-3.5, -10)), aes(y=y, label=lab), x=-0.2, size=3.5, hjust=1, fontface="bold")
    } else {
      gr <- gr + coord_cartesian(ylim=c(-8, max(temp$qual)+5), xlim=c(0,maxLen))
    }
    print(gr)
  }
}
dev.off()

####################
# part 2: heatmaps #
####################
for (EOI in c('3p','5p')){
  seq_FN = paste0(out_DN, '/trimVisTmpFiles/seq',EOI,'sites.txt')
  df2 <- read.table(seq_FN, header = T, colClasses = 'character') # columns of Ts can be interpreted as logical
  colConvert=data.frame(pattern=c('^readID$', 'CutPos', '^[gs][\\d]+$','^q[\\d]+$'), fun=1:4)
  colConvert$fun=c(as.character, as.integer, function(x){factor(x, levels=c('A', 'C','G','T','X','N'))}, as.integer)
  for (rw in 1:nrow(colConvert)){
  	  pattern = colConvert$pattern[rw]
    scols = which(grepl(pattern=pattern, colnames(df2), perl=T))
      for (scol in scols){
  	        df2[,scol] = colConvert$fun[rw][[1]](df2[,scol])
      }
  }
  if (nrow(df2) == 0){
    print(paste0("No reads in ",EOI,"-trimmed class dataset for making aggregate trim plots. Perhaps there are no ",EOI,"-trimmed reads in the data?"))
    next
  }
  if (EOI=='5p'){
    df2$anchor = df2$fpCutPos
  }else if (EOI == '3p'){
    df2$anchor = df2$tpCutPos
  }
  print(paste0('quantiles of ', EOI, '-clipping positions in reads:'))
  print(quantile(df2$anchor, (0:10)/10))
  
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
        # while we're here, reset the zero of quals
        minq = min(smat, na.rm = T)
        maxq = max(smat[smat != 78 & smat != 88], na.rm = T) # get rid of N's and X's inserted into qvals
        zeroQual=0
        if (minq > 32 && maxq < 75){
          print ("Guessing phred encoding as either Illumina 1.8+ Phred+33 (0-41 raw) or Sanger Phred+33 (0-40 raw)")
          zeroQual = 33
        } else {
          print ("Cannot guess phred encoding. Leaving q-vals as-is.")
          print (paste0("Max qval:", maxq, "  Min qval:", minq))
        }
        smat = smat - zeroQual
        for (cl in which(scols)){
          df2[,cl] = df2[,cl] - zeroQual
        }
      }
      hc = hclust(d)
      slong = melt(data = smat)
      colnames(slong)=c('readID', 'pos', prefix)
      slong$pos=as.numeric(substr(slong$pos, 2, 9999)) - 1 - aggFlank
      sreord = data.frame(readID = hc$labels[hc$order], order = 1:length(hc$labels))
      colnames(sreord)[2] = paste0('ord.',prefix)
      blocks[[prefix]] = smat
      lblocks[[prefix]] = slong
      reord[[prefix]] = sreord
      sresults = merge(slong, sreord, by = 'readID')
      if (first){
        allres = sresults
        first=F
      } else {
        allres=merge(allres, sresults, by=c('readID','pos'))
        allres=allres[order(allres$pos),]
      }
    } else {
      absentData = c(absentData, prefix)
    }
  }
  
  sreord = data.frame(df2[!duplicated(df2$readID),c('readID', 'anchor')])
  sreord = merge(sreord, allres[!duplicated(allres$readID),])
  x = order(sreord$anchor, sreord$ord.s)
  sreord=sreord[x,]
  sreord$ord.t=1:nrow(sreord)
  #sreord$ord.t = rank(sreord$tpCutPos, sreord$ord.s, sreord$ord.q, ties.method = 'random')
  
  allres2 = merge(allres, sreord[,c('readID', 'ord.t')], all.x=T, by='readID')
  allres2$q[allres2$s %in% c('N','X')]=NA
  
  
  if (length(absentData) > 0 && absentData == 'g'){
    pdf( paste0(out_DN, '/TVheatmap_S_',EOI,'.pdf'), width=15, height=20)
      grid.arrange(anchorPlot(allres2, 's','s', EOI), anchorPlot(allres2, 'q','s', EOI), ncol=2)
    dev.off()
    pdf( paste0(out_DN, '/TVheatmap_Q_',EOI,'.pdf'), width=15, height=20)
      grid.arrange(anchorPlot(allres2, 's','q', EOI), anchorPlot(allres2, 'q','q', EOI), ncol=2)
    dev.off()
  } else {
    pdf( paste0(out_DN, '/TVheatmap_S_',EOI,'.pdf'), width=15, height=20)
      grid.arrange(anchorPlot(allres2, 's','s', EOI), anchorPlot(allres2, 'q','s', EOI, gdiff), anchorPlot(allres2, 'g','s', EOI, gdiff), ncol=3)
    dev.off()
    pdf( paste0(out_DN, '/TVheatmap_Q_',EOI,'.pdf'), width=15, height=20)
      grid.arrange(anchorPlot(allres2, 's','q', EOI), anchorPlot(allres2, 'q','q', EOI), anchorPlot(allres2, 'g','q', EOI, gdiff), ncol=3)
    dev.off()
    pdf( paste0(out_DN, '/TVheatmap_G_',EOI,'.pdf'), width=15, height=20)
      grid.arrange(anchorPlot(allres2, 's','g', EOI), anchorPlot(allres2, 'q','g', EOI), anchorPlot(allres2, 'g','g', EOI, gdiff), ncol=3)
    dev.off()
  }
  
  ###############################
  # part 3: trm-length profiles #
  ###############################
  
  df3=df2[,c('readID','fpCutPos','tpCutPos','anchor')] 
  df3=df3[order(df3$anchor),]
  df3$vert=1:nrow(df3)
  pdf()
  grL <- ggplot(df3) + geom_segment(aes(x=fpCutPos, y=vert, xend=tpCutPos, yend=vert), size=1) + 
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank())+
                       #axis.text.y=element_blank(),
                       #axis.ticks.y=element_blank(),
                       #axis.title.y=element_blank()) + 
    ylab(label = 'number of reads') + xlab(label = 'position on raw read')
  
  pdf(paste0(out_DN, '/profile_',EOI,'cut.pdf'), width = 6, height = 6)
  print (grL)
  dev.off()

}


