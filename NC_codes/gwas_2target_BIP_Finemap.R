# fine mapping SNP To target gene using ABC (BIP)
# ABC-MAX

library(dplyr)
library(stringr)
library(reshape2)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)
library(data.table)

setwd('~/PF_25BR/data/')

load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata')
ensg_2genesymb <- read.csv('/sc/arion/projects/roussp01a/liting/Pf_25/ensg_2genesymb.csv', row.names = 1)
#load('ABC_summary.Rdata')
#ensg_2genesymb <- read.csv('ensg_2genesymb.csv', row.names = 1)

dista <- 1000

fmp_bip <- read.delim(gzfile('daner_bip_pgc3_only_finemap_only_finemap_all.txt.gz'))%>%subset(PIP !='PIP' )
fmp_bip$PIP <- as.numeric(fmp_bip$PIP)
fmp_bip <- subset(fmp_bip,PIP > 0.01)
fmp_bip <- plyr::rename(fmp_bip,c("CHR"="chr", "BP"='start','PIP'='finemap_posterior_probability'))
fmp_bip$chr <- paste0('chr',fmp_bip$chr)
fmp_bip$start <- as.numeric(fmp_bip$start)
fmp_bip$end <- fmp_bip$start
fmp_bip <- makeGRangesFromDataFrame(fmp_bip, keep.extra.columns = T)

ch = import.chain("/hpc/users/songl05/PF_25BR/data/hg19ToHg38.over.chain")

fmp_bip_38 <- liftOver(fmp_bip,ch)

# to perform closest gene,

write.table(unique(as.data.frame(fmp_bip_38)[, c("seqnames", "start", "end" ,"A1",'A2','SNP' )]),
            sep='\t', file="gwas_finemap_bip.avinput", quote = F,row.names = F,col.names = F)


system('sh /hpc/users/songl05/PF_25BR/codes/annovar_bip.sh')
system('sh /hpc/users/songl05/PF_25BR/codes/find.closest_bip.sh')


anno_hg38 <- read.table("finemap_bip.hg38_multianno.txt", sep='\t',header = T )
closest_gene <- read.table("gwas_finemap_bip.answer.bed", sep='\t' )


closest_gene$closest_G <- str_split(closest_gene$V15,'\\s|;',simplify = T)[,8]
closest_gene$closest_ensG <- str_split(closest_gene$V15,'\\s|;|\\.',simplify = T)[,2]
closest_gene$closest_distance <- str_split(closest_gene$V15,'\\|',simplify = T)[,2]


epLinks <- list()
for (brain_Rg in names(GRs)){
  grs <- as.data.frame(GRs[[brain_Rg]])
  #grs$chr <- grs$seqnames
  grs$start <- merged_peak_500[grs$eID,]@ranges@start - ((dista/2)-250)
  grs$end <- grs$start + dista
  grs$TargetGene_symbol <- ensg_2genesymb[grs$TargetGene,]
  epLinks[[brain_Rg]] <- makeGRangesFromDataFrame(grs[,-4:-6], keep.extra.columns = T)
}

fm_link <- list()
for(br in names(epLinks)){
  hits = findOverlaps(epLinks[[br]], fmp_bip_38)
  results = cbind.data.frame(epLinks[[br]][hits@from,], fmp_bip_38[hits@to,])
  results$BR_Region <- br
  fm_link[[br]] <- results
}

bip_final_link = data.frame(do.call("rbind.data.frame", fm_link))
bip_final_link <- left_join(bip_final_link,closest_gene[,c('V1','V2','closest_ensG','closest_distance')], c("seqnames.1" = "V1","start.1"="V2"))
bip_final_link <- left_join(bip_final_link,anno_hg38[,c('Chr','Start','Func.refGene','Gene.refGene')], c("seqnames.1" = "Chr","start.1"="Start"))

bip_final_link <- bip_final_link%>%group_by(SNP,BR_Region)%>%top_n(1,ABC.Score)
bip_final_link$closest_distance <- abs(as.numeric(bip_final_link$closest_distance))

bip_final_link <- as.data.frame(rename(bip_final_link, c("gene"="TargetGene",'rsid'='SNP')))

save(bip_final_link, file='/sc/arion/projects/roussp01a/liting/Pf_25/output/bip_finemap_final_link.RData')


finemap_mull <- c("SCN2A","TRANK1","DCLK3","INSYN2B","SYNE1","THSD7A","CACNA1B","TUBBP5","PLCB3","PRDX5","KCNK4","AP001453.3","TRPT1","FKBP2","DNAJC4","RASGRP1","FURIN","FES","YWHAE","DPH1","GSDMB","MED24","THRA","EEF1A2","KCNQ2")
