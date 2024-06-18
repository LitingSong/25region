
##################### ldsc
## peak width= 2000
## 5' 
## non 5'
## non-regulating

Dist <- c('1k','2k')[1]

.libPaths(c(.libPaths(),'/sc/arion/projects/roussp01a/jaro/programs/R_libs_4_0_2/'))
primarySourceDir="/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/"
#primarySourceDir="/sc/arion/projects/roussp01a/liting/NYGC_AD_R01/"

source(file.path(primarySourceDir,"STEP4.R"))

root='/sc/arion/projects/roussp01a/liting'

load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata')

#half of length
if(Dist=='1k'){dista <- 1000}
if(Dist=='2k'){dista <- 2000}

enhancer_nonreg <- list()
for (brain_Rg in names(GRs)){
  
  raw_peak_file <- paste0('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/peaks/',brain_Rg,'/',brain_Rg,'.bed')
  raw_peak <- read.table(raw_peak_file)
    colnames(raw_peak)[1:3] <- c('chr','start','end')
  
  br_gr <- makeGRangesFromDataFrame(raw_peak)
  hits = findOverlaps( br_gr, merged_peak_500)
  active_enh <- merged_peak_500[subjectHits(hits),] # active enhancers
  hits = findOverlaps(  GRs$ARC, active_enh)
  non_regul <- active_enh[unique(-subjectHits(hits)), ]
  non_regul <- non_regul[!duplicated(names(non_regul)),]
  non_regul <- as.data.frame(non_regul)
  non_regul$start <-  non_regul$start - ((dista/2)-250)
  non_regul$end <- non_regul$start + dista
  
  non_regul <- makeGRangesFromDataFrame(non_regul)
  
  enhancer_nonreg[[brain_Rg]] <- non_regul
  
}


epLinks_5 <- list()
for (brain_Rg in names(GRs)){
  grs <- as.data.frame(GRs[[brain_Rg]])
  #grs$chr <- grs$seqnames
  grs$start <- merged_peak_500[grs$eID,]@ranges@start - ((dista/2)-250) 
  #grs$end <- merged_peak_500[grs$eID,]@ranges@start + 1250
  grs$end <- grs$start + dista
  grs <- grs[order(grs$order),]
  grs <- grs[!duplicated(grs$ID),]
  grs <- subset(grs, order==1)
  #grs$TargetGene_symbol <- ensg_2genesymb[grs$TargetGene,]
  epLinks_5[[brain_Rg]] <- makeGRangesFromDataFrame(grs[,-4:-6], keep.extra.columns = T)
}


epLinks_n5 <- list()
for (brain_Rg in names(GRs)){
  grs <- as.data.frame(GRs[[brain_Rg]])
  #grs$chr <- grs$seqnames
  grs$start <- merged_peak_500[grs$eID,]@ranges@start - ((dista/2)-250)
  #grs$end <- merged_peak_500[grs$eID,]@ranges@start + 1250
  grs$end <- grs$start + dista
  grs <- grs[order(grs$order),]
  grs <- grs[!duplicated(grs$ID),]
  grs_5 <- subset(grs, order==1)
  grs <- subset(grs, !eID%in%grs_5$eID)
  
  #grs$TargetGene_symbol <- ensg_2genesymb[grs$TargetGene,]
  epLinks_n5[[brain_Rg]] <- makeGRangesFromDataFrame(grs[,-4:-6], keep.extra.columns = T)
}



  
outDir=paste0(root,"/",'ldsc_',Dist,'_isoform')
call_step4 <- function(peaks,name){
  myStep4CellPeaks=peaks
  myStep4CellPeaks=myStep4CellPeaks[sapply(myStep4CellPeaks,length)>2000]
  if(length(myStep4CellPeaks)>0){
    peaksMetadata=data.frame(
      peakSets=names(myStep4CellPeaks),
      peaksFullName=names(myStep4CellPeaks),
      stringsAsFactors=F
    )
    outDir=paste0(root,"/",name)
    step4jobGen( 
      peakSets=myStep4CellPeaks,
      peaksMetadata=peaksMetadata,
      outDir=outDir,
      GENOME_VERSION = "hg38",
      step4cores=12,
      doLdscAnalysis=T, #FIXME:
      doConsAnalysis=T,
      ldscPadding=0,
      STEP4_IS_DRY_RUN=T
    )
  }else{
    return('no peak left')
  }
}

call_step4(epLinks_5,paste0('ldsc_',Dist,'_5'))
call_step4(epLinks_n5,paste0('ldsc_',Dist,'_n5'))
call_step4(enhancer_nonreg,paste0('ldsc_',Dist,'_nonreg'))

## number of elements

enhancer_n <- as.data.frame(cbind(unlist(lapply(enhancer_nonreg,length)),unlist(lapply(epLinks_5,length)),unlist(lapply(epLinks_n5,length))))
colnames(enhancer_n) <- c('nonABC',"5'","non5'")
enhancer_n$brain <- rownames(enhancer_n)
enhancer_n <- melt(enhancer_n)
enhancer_sum <- aggregate(data=enhancer_n,value~brain,sum)

enhancer_stat <- merge(enhancer_n,enhancer_sum, by="brain")%>%mutate(prop = value.x/value.y)


round(mean(subset(enhancer_stat, variable=="non5'")$prop)*100,2)
round((sd(subset(enhancer_stat, variable=="non5'")$prop)/3)*100,2)

round(mean(subset(enhancer_stat, variable=="5'")$prop)*100,2)
round((sd(subset(enhancer_stat, variable=="5'")$prop)/3)*100,2)


# number
pdf(file = "./figures/supp/sp3_n_enhancers_ldsc.pdf", width = 4.5, height = 4.5)

ggplot(enhancer_n, aes(x=brain,y=value, fill=variable))+
 geom_bar(stat = 'identity',position = 'dodge')+
 #facet_wrap(~variable)+
 geom_text(aes(label=value),size=2.5,
           position=position_dodge(width = 1) )+
 scale_fill_manual("Enhancer",values = c("5'"="#FEE090",
                                         "non5'"="#4575B4",
                                         "nonABC"="#FB9A99"))+
 theme_classic()+theme(legend.position = 'top',
                       axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
 )+xlab('') +ylab('Number of elements')

dev.off()

# proportion
pdf(file = "./figures/supp/sp3_p_enhancers_ldsc.pdf", width = 4.5, height = 4.5)

ggplot(enhancer_stat, aes(x=brain,y=prop, fill=variable))+
 geom_bar(stat = 'identity',position = 'stack')+
 scale_fill_manual("Enhancer",values = c("5'"="#FEE090",
                                         "non5'"="#4575B4",
                                         "nonABC"="#FB9A99"))+
 theme_classic()+theme(legend.position = 'top',
                       axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
 )+xlab('') +ylab('Proportion of elements')

dev.off()




