# fine mapping SNP To target gene using ABC (SCZ)
# ABC-MAX

library(dplyr)
library(stringr)
library(reshape2)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)
library(data.table)

load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata')
ensg_2genesymb <- read.csv('/sc/arion/projects/roussp01a/liting/Pf_25/ensg_2genesymb.csv', row.names = 1)
#load('ABC_summary.Rdata')
#ensg_2genesymb <- read.csv('ensg_2genesymb.csv', row.names = 1)

dista <- 1000

fmp_scz <- read.delim("/sc/arion/projects/roussp01a/jaro/data/dual_assay/data/scz_finemap.tsv")%>%subset(finemap_posterior_probability> 0.01)

fmp_scz <- fmp_scz[,c(1,3:5,13,19)]
fmp_scz <- plyr::rename(fmp_scz,c("chromosome"="chr", "position"='start'))
fmp_scz$chr <- paste0('chr',fmp_scz$chr)
fmp_scz$end <- fmp_scz$start
fmp_scz <- makeGRangesFromDataFrame(fmp_scz, keep.extra.columns = T)

ch = import.chain("/hpc/users/songl05/PF_25BR/data/hg19ToHg38.over.chain")

fmp_scz_38 <- liftOver(fmp_scz,ch)

# to perform closest gene,
#write.table( subset(as.data.frame(fmp_scz_38)[,c("seqnames","start",'start')]), row.names = F, col.names = F, sep='\t',quote = F,
#             file='/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/BI_25/data/finemap_scz_38.bed' )


closest_gene <- read.table('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/BI_25/data/answer.sort.bed',sep='\t')
#closest_gene <- read.table('answer.sort.bed',sep='\t')

closest_gene$closest_G <- str_split(closest_gene$V12,'\\s|;',simplify = T)[,8]
closest_gene$closest_ensG <- str_split(closest_gene$V12,'\\s|;|\\.',simplify = T)[,2]
closest_gene$closest_distance <- str_split(closest_gene$V12,'\\|',simplify = T)[,2]



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
  hits = findOverlaps(epLinks[[br]], fmp_scz_38)
  results = cbind.data.frame(epLinks[[br]][hits@from,], fmp_scz_38[hits@to,])
  results$BR_Region <- br
  fm_link[[br]] <- results
}

SCZ_final_link = data.frame(do.call("rbind.data.frame", fm_link))
SCZ_final_link <- left_join(SCZ_final_link,closest_gene[,c('V1','V2','closest_ensG','closest_distance')], c("seqnames.1" = "V1","start.1"="V2"))
SCZ_final_link <- SCZ_final_link%>%group_by(rsid,BR_Region)%>%top_n(1,ABC.Score)
SCZ_final_link$closest_distance <- abs(as.numeric(SCZ_final_link$closest_distance))

SCZ_final_link <- as.data.frame(rename(SCZ_final_link, c("gene"="TargetGene")))

save(SCZ_final_link, file='/sc/arion/projects/roussp01a/liting/Pf_25/output/SCZ_finemap_final_link.RData')

