# ucsc/scz INPUT
# SCZ finemaped SNP

library(stringr)
library(reshape2)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(cowplot)
library(data.table)


load('/sc/arion/projects/roussp01a/liting/Pf_25/output/SCZ_finemap_final_link.RData') # get from gwas_2target_szc_finemap.R
load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata')


# SCZ
SCZ_final_link$internalPromoter <- promoter_anno[SCZ_final_link$TargetTranscript,]$internalPromoter
SCZ_final_link$tss_start <- promoter_anno[SCZ_final_link$TargetTranscript,]@ranges@start
SCZ_final_link$tss_chr <- SCZ_final_link$seqnames
SCZ_final_link$tss_end <- promoter_anno[SCZ_final_link$TargetTranscript,]@ranges@start +1 
SCZ_final_link$tss_strand <- as.data.frame(promoter_anno)[SCZ_final_link$TargetTranscript,'strand']

# finemap
fmp_scz <- read.delim("/sc/arion/projects/roussp01a/jaro/data/dual_assay/data/scz_finemap.tsv")%>%subset(finemap_posterior_probability> 0.01)
rownames(fmp_scz) <- fmp_scz$rsid


# pops score
scz_pops <- read.table('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/POPS/files/final_pops/added_ridge_sz3.preds', sep='\t', header = T, row.names =1 )
scz_pops$quantile <- rank(scz_pops$PoPS_Score)/length(scz_pops$PoPS_Score) 
ensg_2genesymb <- read.csv('/sc/arion/projects/roussp01a/liting/Pf_25/ensg_2genesymb.csv', row.names = 1)

SCZ_final_link$quantile <- scz_pops[SCZ_final_link$gene,'quantile']
SCZ_final_link$pops <- scz_pops[SCZ_final_link$gene,'quantile']

# prioritized gene/isoform
SCZ_final_link$finemap_posterior_probability <- fmp_scz[SCZ_final_link$rsid,'finemap_posterior_probability']

scz_target_supple <- unique(SCZ_final_link[,c('rsid','finemap_posterior_probability','gene','TargetGene_symbol','TargetTranscript','order','quantile','BR_Region')])
scz_target_supple <- rename(scz_target_supple,c("PoPS_quantile"="quantile", "TargetGene"='gene'))
scz_target_supple$PoPS_quantile <- round(scz_target_supple$PoPS_quantile, 3)

write.table( unique(scz_target_supple), file='/sc/arion/projects/roussp01a/liting/Pf_25/data/supple_tables/prioritized_genes_scz.txt',quote = F,sep='\t',row.names = F)

 
# Prepare data for ucsc/scz genome brower 
#Candi <- subset(SCZ_final_link, TargetGene_symbol%in%c('SLC66A2','CALN1','MYT1L','PLCH2','TBL1XR1')) 
Candi <- subset(SCZ_final_link, rsid%in%c("rs11660941","rs72980082",'rs2944829', "rs6687012","rs6688934","rs2494638")|TargetGene_symbol=='ABHD2')
Candi <- SCZ_final_link

## candidate snp bed file
candi_snp <- unique(Candi[,c('seqnames.1','start.1','end.1','rsid')])
candi_snp$start.1 <- candi_snp$start.1 - 1
candi_snp <- cbind(candi_snp, score=1000,strand='.',candi_snp[,c('start.1','end.1')],itemRgb='241,66,60')
write.table(candi_snp, file='/hpc/users/songl05/PF_25BR/data/ucsc/scz/candi_snp.txt',sep='\t',quote = F,row.names = F, col.names = F)

candi_snp$PIP <- fmp_scz[candi_snp$rsid,'finemap_posterior_probability']
write.table( candi_snp[,c('seqnames.1','start.1','end.1','PIP')],  
             file='/hpc/users/songl05/PF_25BR/data/ucsc/scz/candi_snp_pip.txt',sep='\t',quote = F,row.names = F, col.names = F )


## candidate ABC links
## track type=interact name="interact Example One" description="An interact file" interactDirectional=true maxHeightPixels=200:100:50 visibility=full
Br_color <- read.table('/hpc/users/songl05/PF_25BR/data/color_br.txt',sep='\t',comment.char = '', header = T, row.names = 1)

Candi$enhancer_start1 <- Candi$start
Candi$enhancer_end1 <- Candi$end -1 

link_start <- apply(Candi[,c('tss_start','tss_end','enhancer_start1','enhancer_end1')],1,min)
link_end <- apply(Candi[,c('tss_start','tss_end','enhancer_start1','enhancer_end1')],1,max)

link_inter <- cbind(chr=Candi$seqnames.1,link_start, link_end, ID=Candi$ID, uint=500,
                    value= Candi$ABC.Score*10, exp=Candi$ID,color=Candi$BR_Region,
                    Candi[,c('seqnames.1','enhancer_start1','enhancer_end1')],
                    sourceName=Candi$eID,sourceStrand=Candi$strand, targetChrom= Candi$seqnames.1, 
                    Candi[,c('tss_start','tss_end')],targetName=paste0('I',Candi$order),targetStrand=Candi$tss_strand)


link_inter$color <- Br_color[link_inter$color,'rgb']
write.table(unique(link_inter), file='/hpc/users/songl05/PF_25BR/data/ucsc/scz/link_inter.txt',sep='\t',quote = F,row.names = F, col.names = F)


## promoter isoforms bed
###track name="candidate snp" visibility=2 color="255,0,0"
promoter_isoform_tss <- as.data.frame(subset(promoter_anno, gene_id%in%unique(Candi$gene)))[,c('seqnames','start','end','order','score','strand','start','end')]
promoter_isoform_tss$order <- paste0('I',promoter_isoform_tss$order)
promoter_isoform_tss$score <- 500
promoter_isoform_tss$itemRgb <- ifelse(rownames(promoter_isoform_tss)%in%Candi$TargetTranscript,'227,33,29','0,0,0')
promoter_isoform_tss$itemRgb[promoter_isoform_tss$order==1] <- '56,126,184'
write.table(promoter_isoform_tss, file='/hpc/users/songl05/PF_25BR/data/ucsc/scz/promoter_isoform_tss.txt',sep='\t',quote = F,row.names = F, col.names = F)


# GWAS 
library(data.table)
gwas_scz <- fread('/hpc/users/songl05/PF_25BR/data/PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv',select = c(1,2,3,11))
gwas_scz$CHROM <- paste0('chr',gwas_scz$CHROM)
lead_hg19 <-  unique(subset(gwas_scz,ID%in%Candi$index.snp)[,1:4])
gwas_scz_candi <- c()
for (i in 1:nrow(lead_hg19)){
  gwas_scz_candi <- rbind(gwas_scz_candi, subset(gwas_scz, CHROM==as.character(lead_hg19[i,1]) & 
                                                   POS > as.numeric(lead_hg19[i,3])-1000000  & 
                                                   POS < as.numeric(lead_hg19[i,3]) + 1000000))
  print( paste0(lead_hg19[i,1],':',as.numeric(lead_hg19[i,3])-1000000,'-',as.numeric(lead_hg19[i,3]) + 1000000))
}
ch = import.chain("/hpc/users/songl05/PF_25BR/data/hg19ToHg38.over.chain")

gwas_19 = with(gwas_scz_candi, GRanges(CHROM, IRanges(start=as.numeric(POS), end=as.numeric(POS) ), PVAL=PVAL,ID=ID ))
gwas_pos_38 <- as.data.frame(liftOver(gwas_19,ch))[,c(3:5,8)]
gwas_pos_38$PVAL <-  -log10(gwas_pos_38$PVAL)
write.table(gwas_pos_38, file='/hpc/users/songl05/PF_25BR/data/ucsc/scz/gwas_pos_38_scz.txt',sep='\t',quote = F,row.names = F, col.names = F)

## lead_snp bed file
#track name="Lead snp" visibility=2 color="255,0,0"
index_snp <- unique(Candi$index.snp)
gwas_pos_38  <- as.data.frame(liftOver(gwas_19,ch))
lead_snp <- subset(gwas_pos_38, ID %in%index_snp)[,c('seqnames','start','end','ID')]
lead_snp <- cbind(lead_snp, score=1000,strand='.',lead_snp[,c('start','end')],itemRgb='76,117,156')
write.table(lead_snp, file='/hpc/users/songl05/PF_25BR/data/ucsc/scz/lead_snp.txt',sep='\t',quote = F,row.names = F, col.names = F)


# TF disrupt

# to identify genes with disrupted footprint and TF motif


# MOTIF break (SNP)
load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/test/neuron_norm_none_BIC_4_in_0.05_CPM_1_in_0.1/files/SCZ_finemap_.Rdata' )
motifBrk <- result

# atac-footprint
atac_footPrt <- read.table(header = T,'/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/tobias/bindetect/bindetect_results.txt' )

# identify snp both disrupted tf motif and corespnsing tf footprint
mtf_footp <- subset(motifBrk,geneSymbol %in% str_split(atac_footPrt$output_prefix,'_',simplify = T)[,1] )

xxx <- subset(SCZ_final_link,rsid%in%names(mtf_footp)  #& order!=1 
              &gene %in%rownames(scz_pops)[(scz_pops$quantile> 0.9)])

unique(xxx$TargetGene_symbol)

# find example: rs and gene

ft_path <- '/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/tobias/bindetect/'
for(RS in unique(xxx$trait_buddy_rs)){
  # RS <- 'rs35709455'
  gene <- subset(xxx,rsid==RS )
  TF <- subset(mtf_footp,SNP_id==RS)$geneSymbol
  footPrt_n <- dir(ft_path,pattern = TF)
  footPrt <- read.table(paste0(ft_path,footPrt_n,'/',footPrt_n, "_overview.txt"),header = T)
  overlap11 <- findOverlaps(makeGRangesFromDataFrame(footPrt[,1:3]), makeGRangesFromDataFrame(gene[1,4:6]))
  if( length(overlap11) >= 1){print(c(gene[,'TargetGene_symbol'],RS))}
}


# snp to atac-seq
ML_LI <- read.table('/sc/arion/projects/roussp01a/pengfei/ref/BICCN/Li_science_2023_multibrainregion_atacseq/Schizophrenia.99credset.PPA.hg19ToHg38.bdg')
ML_LI$id <- paste0(ML_LI$V1,":",ML_LI$V3)

subset(ML_LI, id%in%"chr3:180898828")

atac_footPrt <- read.table(header = T,'/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/tobias/bindetect//NR0B1_M6381_1.02/NR0B1_M6381_1.02_overview.txt' )


















target_n5 <- c()
for (Br in unique(SCZ_final_link$BR_Region)){
  grs <- subset(SCZ_final_link, BR_Region==Br)
  grs <- grs[order(grs$order),]
  grs_uniqueid <- grs[!duplicated(grs$ID),]
  grs_5 <- subset(grs_uniqueid, order==1)
  grs <- subset(grs, !ID%in%grs_5$ID)
  
  target_n5 <- rbind(target_n5,grs)
}

library(reshape2)

gene_brain_region <- dcast(data =target_n5, formula = gene~ BR_Region,length)
rownames(gene_brain_region) <- gene_brain_region$gene
gene_brain_region <- gene_brain_region[,-1]
gene_brain_region <- gene_brain_region[!is.na(ensg_2genesymb[rownames(gene_brain_region),]),]

gene_brain_region[gene_brain_region > 1] <- 1
gene_brain_region <- gene_brain_region[rowSums(gene_brain_region)<=3,]
#gene_brain_region <- gene_brain_region[rowSums(gene_brain_region)==9,]


