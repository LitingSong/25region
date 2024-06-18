##################### ldsc for ABC promoter
## peak width= 1000
## 5' 
## non 5'


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

dista <- 1000
active_promoter_nonreg <- list()
active_promoter_5 <- list()
active_promoter_n5 <- list()



for (brain_Rg in names(GRs)){
  
  grs_exp <- as.data.frame(exp_promoters[[brain_Rg]])
  grs_exp$start_tss <- ifelse(grs_exp$strand=='+',grs_exp$start,grs_exp$end )
  
  grs_exp$start <- grs_exp$start_tss - (dista/2)
  grs_exp$end <- grs_exp$start + dista
  grs_exp$strand <- '*'
  
  grs <- as.data.frame(GRs[[brain_Rg]])
  grs <- grs[!duplicated(grs$TargetTranscript),]
  grs$start <- grs$TargetTranscriptTSS - (dista/2)
  grs$end <- grs$start + dista
  
  grs_5 <- subset(grs,order==1)
  grs_n5 <- subset(grs,order!=1)
  grs_nonreg <- subset(grs_exp,!symbol %in% grs$TargetTranscript )
    
  active_promoter_5[[brain_Rg]] <- makeGRangesFromDataFrame(grs_5, keep.extra.columns = T)
  active_promoter_n5[[brain_Rg]] <- makeGRangesFromDataFrame(grs_n5, keep.extra.columns = T)
  active_promoter_nonreg[[brain_Rg]] <- makeGRangesFromDataFrame(grs_nonreg, keep.extra.columns = T)
  
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

call_step4(active_promoter_5,paste0('ldsc_',Dist,'_Promoter_5'))
call_step4(active_promoter_n5,paste0('ldsc_',Dist,'_Promoter_n5'))
call_step4(active_promoter_nonreg,paste0('ldsc_',Dist,"_Promoter_nonreg"))

# number of elements
promoter_n <- as.data.frame(cbind(unlist(lapply(active_promoter_nonreg,length)),
                                 unlist(lapply(active_promoter_5,length)),
                                 unlist(lapply(active_promoter_n5,length))))
colnames(promoter_n) <- c('nonABC',"5'","non5'")
promoter_n$brain <- rownames(promoter_n)
promoter_n <- melt(promoter_n)

promoter_sum <- aggregate(data=promoter_n,value~brain,sum)
promoter_stat <- merge(promoter_n,promoter_sum, by="brain")%>%mutate(prop = value.x/value.y)

round(mean(subset(promoter_stat, variable=="non5'")$prop)*100,2)
round((sd(subset(promoter_stat, variable=="non5'")$prop)/3)*100,2)

round(mean(subset(promoter_stat, variable=="5'")$prop)*100,2)
round((sd(subset(promoter_stat, variable=="5'")$prop)/3)*100,2)

pdf(file = "./figures/sp3_n_promoter_ldsc.pdf", width = 4.5, height = 4.5)

ggplot(promoter_n, aes(x=brain,y=value, fill=variable))+
 geom_bar(stat = 'identity',position = 'dodge')+
 #facet_wrap(~variable)+
 geom_text(aes(label=value),size=2,
           position=position_dodge(width = 1) )+
 scale_fill_manual("Promoter",values = c("5'"="#FEE090",
                                         "non5'"="#4575B4",
                                         "nonABC"="#FB9A99"))+
 theme_classic()+theme(legend.position = 'top',
                       axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
 )+xlab('') +ylab('Number of elements')

dev.off()


pdf(file = "./figures/sp3_p_Promoter_ldsc.pdf", width = 4.5, height = 4.5)

ggplot(promoter_stat, aes(x=brain,y=prop, fill=variable))+
 geom_bar(stat = 'identity',position = 'stack')+
 scale_fill_manual("Promoter",values = c("5'"="#FEE090",
                                         "non5'"="#4575B4",
                                         "nonABC"="#FB9A99"))+
 theme_classic()+theme(legend.position = 'top',
                       axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)
 )+xlab('') +ylab('Proportion of elements')

dev.off()

