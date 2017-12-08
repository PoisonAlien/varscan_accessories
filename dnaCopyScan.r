#Usage: From command line		
#Rscript dnaCopyScan.r <foo.copynumber>
#Author: Anand Mayakonda (Cancer Science Institute of Singapore; NUS)

dnaCopy=function(copynumber_file){
  suppressPackageStartupMessages(library(DNAcopy))
  suppressPackageStartupMessages(library(data.table))
  tn = data.table::fread(copynumber_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
  colnames(tn)[1] = 'contig'
  tn$contig = gsub(pattern = 'chr', replacement = '', x = tn$contig, fixed = T)
  tnxy = tn[contig %in% c('X', 'Y')]
  tn = tn[!contig %in% c('X', 'Y')]
  #tn = tn[!contig == 'Y']
  tn = tn[order(as.numeric(tn$contig)),]
  tn = rbind(tn, tnxy)
  tn = as.data.frame(tn)
  #samp.name = colnames(tn)[5] 
  
  samp.name = gsub(pattern = '.copynumber', replacement = '', x = copynumber_file)
  
  
  cn = DNAcopy::CNA(genomdat = as.numeric(tn[,7]),chrom = as.character(tn[,1]),maploc = as.numeric(tn[,2]), data.type="logratio", sampleid= basename(samp.name), presorted = T)
  cn = DNAcopy::smooth.CNA(cn)
  cn = DNAcopy::segment(cn, alpha = 0.01, nperm = 10000, p.method = 'hybrid', min.width = 2, kmax = 25, nmin = 210, 
             eta = 0.05, trim = 0.025, undo.SD = 3, undo.prune = 0.05, undo.splits = 'none')
  
  tiff(filename = paste(samp.name,"tiff",sep="."), width = 12, height = 4, pointsize = 9, units = 'in', bg = 'white', res = 100)
  plot(cn, ylim = c(-5,5), pt.col = colors()[c(16,17)], segcol = 'midnightblue')
  abline(h=0,lty=1,lwd=0.3)
  abline(h=c(-1,1),lty=2,lwd=0.3)
  dev.off()
  
  segs = cn$output
  colnames(segs) = c("Sample",'Chromosome','Start','End','Num_Probes','Segment_Mean')
  write.table(segs, paste(samp.name, '_cbs.seg', sep=''), quote = F, row.names = F, sep='\t')
  
  save(cn, file = paste(samp.name, '_cbs.RData', sep=''))
}

args = commandArgs(trailingOnly = TRUE)
dnaCopy(copynumber_file = args[1])
