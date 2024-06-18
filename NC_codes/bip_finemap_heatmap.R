# fine mapping SNP To target gene using ABC (bip)
# ABC-MAX

library(RColorBrewer)
library(dplyr)
library(stringr)
library(reshape2)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)
library(data.table)

load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata')
ensg_2genesymb <- read.csv('/sc/arion/projects/roussp01a/liting/Pf_25/ensg_2genesymb.csv', row.names = 1)
#source('/hpc/users/songl05/PF_25BR/codes/gwas_2target_BIP_Finemap.R')
load('/sc/arion/projects/roussp01a/liting/Pf_25/output/bip_finemap_final_link.RData')

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

# we prioritized xx genes across brain regions for xx out of xxx fine-mapped SNP overlapped with ABC enhancers
length(unique(bip_final_link$rsid))
length(unique(bip_final_link$gene))
length(fmp_bip)

# pops score
bip_pops <- read.table('/sc/arion/projects/roussp01a/liting/Pf_25/PF_25BR/data/added_ridge_bip2.preds', sep='\t', header = T, row.names =1 )
bip_pops$quantile <- rank(bip_pops$PoPS_Score)/length(bip_pops$PoPS_Score) 
bip_pops$gene <- ensg_2genesymb[rownames(bip_pops),'gene_name']

# closest_gene
closest_gene <- read.table("/sc/arion/projects/roussp01a/liting/Pf_25/PF_25BR/data/gwas_finemap_bip.answer.bed", sep='\t' )
closest_gene$closest_G <- str_split(closest_gene$V15,'\\s|;',simplify = T)[,8]
closest_gene$closest_ensG <- str_split(closest_gene$V15,'\\s|;|\\.',simplify = T)[,2]
closest_gene$closest_distance <- str_split(closest_gene$V15,'\\|',simplify = T)[,2]


dista <- 1000

target_5 <- subset(bip_final_link, order==1)
target_n5 <-subset(bip_final_link, order!=1)

target_n5_gene <- dcast(data =target_n5, formula = gene~ BR_Region,length) 
target_n5_gene[,-1] <- ifelse(target_n5_gene[,-1]==0,0,1) # non5'
target_5_gene <- dcast(data =target_5, formula = gene~ BR_Region,length)
target_5_gene[,-1] <- ifelse(target_5_gene[,-1]==0,0,2) # 5'


#target_gene <- left_join(target_n5_gene, target_5_gene,by='gene')
target_gene <- merge(target_n5_gene, target_5_gene,by='gene',all.x=T,all.y=T)

rownames(target_gene) <- target_gene$gene
target_gene <- target_gene[,-1]
target_gene[is.na(target_gene)] <- 0
target_isoform <- target_gene[,c(1:9)] + target_gene[,c(10:18)] # 0: none ;1: non5'; 2:5'; 3: both;
colnames(target_isoform) <- gsub('\\.x','',colnames(target_isoform))


col_fun <- colorRamp2(c(0:3), c("#FFFFFF","#4393C3FF","#FEE090FF","#92C5DEFF"))
fill_col <- col_fun(c(1:3))
Br_color <- read.table('/sc/arion/projects/roussp01a/liting/Pf_25/PF_25BR/data/color_br.txt',sep='\t',comment.char = '', header = T, row.names = 1)

#target_isoform <- target_isoform[target_isoform$BasalGanglia>0,]
{ 
  p1 <- Heatmap(target_isoform[,rownames(Br_color)], col = col_fun, 
                row_labels = ensg_2genesymb[rownames(target_isoform),],
                column_labels = Br_color$subRG,
                row_names_gp=gpar(fontsize = 9),
                column_names_gp=gpar(fontsize = 9),
                cluster_columns = F,
                show_heatmap_legend=F,
                column_title = 'ABC target promoter',
                column_title_gp = gpar(fill = "#E8F0F4",fontsize = 9, col = "black", border = "#E8F0F4"),
                rect_gp = gpar(col = "gray", lwd = 0.2),show_column_dend = F, show_row_dend = F )
  
  lgd1  <- Legend(labels = c("Non 5'","5'",'Both'), 
                  title = "ABC target", direction = "horizontal",
                  legend_gp = gpar(fill = fill_col,col='black'))#,#nr=1)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PANEL 2: snp- isoform(5', non5') # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PANEL 3:   differential expression  # # #  
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
{
  #load("DE_isoform.RData")
  load('/hpc/users/songl05/PF_25BR/data/DE_isoform.RData')
  select_vs <- c('Forebrain','MidDen','BasalGanglia','Limbic__Neocortex','MDT__Dien',
                 'RMTG__Midbrain','DRN__Midbrain','HAB__Dien','ARC__Dien','Forebrain__MidDen')
  
  for (vs in select_vs){
    if(nrow(Nup[[vs]])>0) {Nup[[vs]]$comp <- vs}
    if(nrow(Ndown[[vs]])>0) {Ndown[[vs]]$comp <- vs}
    
  }
  
  
  up_trans <- do.call(rbind, Nup[select_vs])
  down_trans <- do.call(rbind, Ndown[select_vs])
  up_trans <- subset(up_trans, PeakID %in%bip_final_link$TargetTranscript)
  up_trans <- data.table::dcast(up_trans,gene_id~comp)
  rownames(up_trans) <- up_trans$gene_id
  up_trans <- up_trans[,-1]
  up_trans <- up_trans[rownames(target_isoform),]
  rownames(up_trans) <- rownames(target_isoform)
  #up_trans[is.na(up_trans)] <- 0
  up_trans[up_trans>=1] <- "Y"
  up_trans[up_trans==0] <- NA
  
  #up_trans <- up_trans[,intersect(colnames(up_trans), c('BasalGanglia','Forebrain','MidDen','Limbic__Neocortex',
  #                        "ARC__Dien","HAB__Dien","MDT__Dien","DRN__Midbrain"))]
  up_trans <- up_trans[, c('BasalGanglia','Forebrain','MidDen')]
  p3 <- Heatmap(up_trans, col = structure("#4393C3FF",names="Y"), cluster_columns = F,cluster_rows = F,
                row_labels = ensg_2genesymb[rownames(up_trans),],
                na_col = 'white',
                column_labels = c('BasGan','Forebrain','MidDen'),
                column_title = 'DE isoform',
                row_names_gp=gpar(fontsize = 9),
                column_names_gp=gpar(fontsize = 9),
                column_title_gp = gpar(fill = "#E8F0F4", fontsize = 9, col = "black", border = "#E8F0F4"),
                show_heatmap_legend=F,rect_gp = gpar(col = "gray", lwd = 0.2),
                show_column_dend = F, show_row_dend = F)
  
  lgd3  <- Legend(labels = c('Y'), title = "DE promoter isoform", direction = "horizontal",
                  legend_gp = gpar(fill = fill_col[1]))#,#nr=1)
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PANEL 4:   eQTL # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #
{ 
  fm_38 <- as.data.frame(fmp_bip_38)
  fm_38$index <- paste(fm_38$seqnames, fm_38$start, sep=":")
  
  #eqtl_dir <- './'
  eqtl_dir <- '/hpc/users/songl05/PF_25BR/data/GTEx_Analysis_v8_eQTL/'
  
  eqtl_brain <- dir(path=eqtl_dir, pattern = 'Brain.*signif_variant_gene_pairs.txt.gz',full.names = F)
  #eqtl_brain <- eqtl_brain[ !grepl('Cerebell|Spinal' ,eqtl_brain) ]
  
  eqtl_target  <- c()
  for (i in 1: length(eqtl_brain)){
    eqtl_sig <- read.delim(gzfile(paste0( eqtl_dir, eqtl_brain[i])), header = T)[,c('gene_id','variant_id')]
    eqtl_sig[,c('chr','pos')] <- str_split(eqtl_sig$variant_id, '_',simplify = T)[,1:2]
    
    eqtl_sig$index <- paste(eqtl_sig$chr, eqtl_sig$pos, sep=":")
    
    e_gene <- unique(subset(eqtl_sig, index%in%fm_38$index)[,'gene_id'])
    eqtl_target <- rbind(eqtl_target,cbind(e_gene, eqtl_brain[i]))
  }
  
  eqtl_target <- as.data.frame(eqtl_target)
  e_target <- dcast(eqtl_target,e_gene~V2)
  rownames(e_target) <- str_split(e_target$e_gene,'\\.',simplify = T)[,1]
  e_target <- e_target[,-1]
  colnames(e_target ) <- str_split_fixed(str_split(colnames(e_target ) ,'\\.',simplify = T)[,1],'_',2)[,2]
  
  
  e_target <- as.data.frame(e_target)[rownames(target_isoform),]
  rownames(e_target) <- rownames(target_isoform)
  e_target <- ifelse(!is.na(e_target),1,0)
  
  #col_fun <- colorRamp2(c(0:3), c("#FFFFFF","#91bfdb","#ffffbf","#fc8d59"))
  
  p4 <- Heatmap(e_target[,c("Frontal_Cortex_BA9",'Cortex','Hippocampus','Amygdala','Caudate_basal_ganglia','Nucleus_accumbens_basal_ganglia','Putamen_basal_ganglia','Hypothalamus')], 
                col = col_fun, 
                column_title = 'eQTL',
                row_names_gp=gpar(fontsize = 9),
                column_names_gp=gpar(fontsize = 9),
                column_title_gp = gpar(fill = "#E8F0F4", fontsize = 9, col = "black", border = "#E8F0F4"),
                row_labels = ensg_2genesymb[rownames(e_target),],cluster_columns = F,cluster_rows = F,
                show_heatmap_legend=F,rect_gp = gpar(col = "gray", lwd = 0.2), show_column_dend = F, show_row_dend = F)
  
  e_target_combined <- as.data.frame(e_target)
  e_target_combined$eqtl <- rowSums(e_target_combined==1)
  e_target_combined$eqtl <- ifelse(e_target_combined$eqtl>=1,1,0)
  
  p4_combined <- Heatmap(e_target_combined[,c("eqtl")], 
                         col = col_fun, 
                         column_title = 'eQTL',
                         column_labels = 'Brain',
                         row_names_gp=gpar(fontsize = 9),
                         column_names_gp=gpar(fontsize = 9),
                         column_title_gp = gpar(fill = "#E8F0F4", fontsize = 9, col = "black", border = "#E8F0F4"),
                         row_labels = ensg_2genesymb[rownames(e_target_combined),],cluster_columns = F,cluster_rows = F,
                         show_heatmap_legend=F,rect_gp = gpar(col = "gray", lwd = 0.2), show_column_dend = F, show_row_dend = F)
  
  lgd4  <- Legend(labels = c('Y'), title = "eQTL", direction = "horizontal",
                  legend_gp = gpar(fill = fill_col[1]))#,#nr=1)
  
  
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PANEL 5:   pathway   # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

#write.table(row.names(target_isoform), file='target.genes.txt', row.names = F, col.names = F, quote = F)
# run metascape, and read meta result
{ 
  meta_enrich <- (read.delim('meta_enrich_all.txt') %>% subset(grepl('Summary',GroupID) & Category=='GO Biological Processes'))[1:10, ]
  rownames(meta_enrich)<- meta_enrich$Term
  meta_anno <- read.delim('meta_anno_all.txt',check.names = F, row.names = 1)
  
  colnames(meta_anno) <- str_split(colnames(meta_anno),' ',simplify = T)[,1]
  meta_anno <- meta_anno[,rownames(meta_enrich)]
  meta_anno[is.na(meta_anno)] <- 0
  p5 <- Heatmap(meta_anno[row.names(target_isoform),], 
                rect_gp = gpar(col = "gray", lwd = 0.2),
                column_title = 'GO: BP',
                show_heatmap_legend=F,
                row_names_gp=gpar(fontsize = 9),
                column_names_gp=gpar(fontsize = 9),
                column_title_gp = gpar(fill = "#E8F0F4", fontsize = 9, col = "black", border = "#E8F0F4"),
                column_labels = meta_enrich[colnames(meta_anno),'Description'],
                col = col_fun, cluster_columns = F,cluster_rows = F)
  
  lgd5  <- Legend(labels = c('Y'), title = "GO: BP", direction = "horizontal",
                  legend_gp = gpar(fill = fill_col[1]))#,#nr=1)
  
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PANEL 6: distance:  the distance between candi ld snp (enhancer) and target isoform (TargetTranscript) # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
{ 
  bip_link_ep <- bip_final_link[,c('gene','TargetTranscript','eID','order','BR_Region','start.1')]
  bip_link_ep$prom_pos <- promoter_anno[bip_link_ep$TargetTranscript,]@ranges@start
  bip_link_ep$enhancer_pos <- merged_peak_500[bip_link_ep$eID,]@ranges@start
  bip_link_ep$ep_dist <- abs(bip_link_ep$enhancer_pos - bip_link_ep$prom_pos)  
  bip_link_ep$vp_dist <- abs(bip_link_ep$start.1 - bip_link_ep$prom_pos)  
  
  vp_dist_max <- aggregate(vp_dist~ gene + BR_Region ,bip_link_ep,mean)
  vp_dist_max <- dcast(vp_dist_max, gene~BR_Region, value=vp_dist)
  rownames(vp_dist_max) <- vp_dist_max$gene
  vp_dist_max <- vp_dist_max[,-1]
  vp_dist_max <- vp_dist_max[rownames(target_isoform),]
  max_d <- log10(max(vp_dist_max, na.rm = T))
  min_d <- log10(min(vp_dist_max, na.rm = T))
  #vp_dist_max[is.na(vp_dist_max)] <- 1
  
  p6 <- Heatmap(log10(vp_dist_max[,rownames(Br_color)]),na_col = "white",
                col=colorRamp2(rev(seq(min_d,max_d,length=5)), brewer.pal(n = 5, name = "RdYlBu")),
                cluster_rows = F, cluster_columns = F,
                row_labels = ensg_2genesymb[rownames(target_isoform),],
                rect_gp = gpar(col = "gray", lwd = 0.2),
                row_names_gp=gpar(fontsize = 9),
                column_names_gp=gpar(fontsize = 9),
                column_labels = Br_color$subRG,
                column_title = 'E-P distance',
                column_title_gp = gpar(fill = "#E8F0F4", fontsize=9,col = "black", border = "#E8F0F4"),
                show_heatmap_legend=F)
  
  lgd6  <- Legend(col_fun = colorRamp2(rev(seq(min_d,max_d,length=5)), 
                                       brewer.pal(n = 5, name = "RdYlBu")), 
                  title = "log10(distance)",
                  direction = "horizontal")
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PANEL 6b: abc score
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
{ 
  bip_link_abc <- bip_final_link[,c('gene','order','ABC.Score','BR_Region')]
  
  ep_abc_mean <- aggregate(ABC.Score~ gene + BR_Region ,bip_link_abc,mean)
  ep_abc_mean <- dcast(ep_abc_mean, gene~BR_Region, value=ABC.Score)
  rownames(ep_abc_mean) <- ep_abc_mean$gene
  ep_abc_mean <- ep_abc_mean[,-1]
  ep_abc_mean <- ep_abc_mean[rownames(target_isoform),]
  
  p6b <- Heatmap(ep_abc_mean[,rownames(Br_color)],na_col = "white",
                 col=colorRamp2(rev(seq(0,0.1,length=5)), 
                                brewer.pal(n = 5, name = "RdYlBu")),
                 cluster_rows = F, cluster_columns = F,
                 row_labels = ensg_2genesymb[rownames(target_isoform),],
                 rect_gp = gpar(col = "gray", lwd = 0.2),
                 row_names_gp=gpar(fontsize = 9),
                 column_names_gp=gpar(fontsize = 9),
                 column_labels = Br_color$subRG,
                 column_title = 'ABC Score',
                 column_title_gp = gpar(fill = "#E8F0F4", fontsize=9,col = "black", border = "#E8F0F4"),
                 show_heatmap_legend=F)
  
  lgd6b  <- Legend(col_fun =colorRamp2(rev(seq(0,0.1,length=5)), 
                                       brewer.pal(n = 5, name = "RdYlBu")),
                   at=c(0,0.05,0.1),
                   title = "ABC Score",
                   direction = "horizontal")
  
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# PANEL 7:   Pops score           # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #


{
  p7 <- Heatmap(bip_pops[rownames(target_isoform),'quantile'],na_col = 'white',
                cluster_rows = F, cluster_columns = F,
                column_labels = c('PoPS'),
                rect_gp = gpar(col = "gray", lwd = 0.2),
                column_title = 'PoPS',
                show_heatmap_legend=F,
                row_names_gp=gpar(fontsize = 9),
                column_names_gp=gpar(fontsize = 9),
                column_title_gp = gpar(fill = "#E8F0F4", fontsize = 9, col = "black", border = "#E8F0F4"),
                row_labels = ensg_2genesymb[rownames(target_isoform),],
                colorRamp2(rev(seq(0,1,length=5)), brewer.pal(5, "RdYlBu"))) 
  
  lgd7  <- Legend(col_fun = colorRamp2(rev(seq(0,1,length=5)), 
                                       brewer.pal(5, "RdYlBu")) ,
                  at=c(0,0.5,1),
                  title = "PoPS",
                  direction = "horizontal")
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # #
## PANEL 8: closest gene
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

clos <- unique(bip_final_link[,c("TargetGene_symbol","gene","rsid",
                                 "Func.refGene",'closest_ensG')])
clos$same_G <- ifelse(clos$gene%in%closest_gene$closest_ensG,1,2)
clos <- clos%>%dplyr::group_by(gene) %>%dplyr::mutate(impact_m=paste0(sort(unique(Func.refGene)), collapse = ','))
#clos <- clos%>%group_by(gene) %>% mutate(same_G_m=paste0(sort(unique(same_G)), collapse = ','))
clos <- as.data.frame(unique(clos[,c('gene','impact_m','same_G')]))
clos$impact_m <- gsub('_variant','',clos$impact_m)
clos$impact_m[!clos$impact_m%in%c('intronic','downstream',
                                  'upstream','intergenic','exonic',
                                  'UTR5','UTR3')] <- 'Others'
rownames(clos) <- clos$gene


p8_1 <- Heatmap(clos[rownames(target_isoform),"same_G"],
                col = col_fun, 
                column_labels = c('Closest gene'),
                column_names_gp=gpar(fontsize = 9),
                rect_gp = gpar(col = "gray", lwd = 0.2),show_column_dend = F, show_row_dend = F,
                column_title_gp = gpar(fill = "#E8F0F4", fontsize = 9, col = "black", border = "#E8F0F4"),
                show_heatmap_legend=F, column_title = 'Prox.')#'Closest G')

lgd8_1  <- Legend(labels = c("Y","N"), 
                  title = "Closest gene", direction = "horizontal",
                  legend_gp = gpar(fill = fill_col[1:2]))#,#nr=1)



p8_2 <- Heatmap(clos[rownames(target_isoform),"impact_m"], 
                #col=structure(brewer.pal(8,'Set3'), names=unique(clos$impact_m)),
                col=structure(brewer.pal(9,'Set3')[-4], names=c('intronic','UTR5','Others','UTR3','upstream','intergenic','exonic')),
                column_labels = c('SNP Anno'),
                column_names_gp=gpar(fontsize = 9),
                column_title_gp = gpar(fill = "#E8F0F4", fontsize = 9, col = "black", border = "#E8F0F4"),
                show_heatmap_legend=F, column_title = 'Func')#Variant')

lgd8_2  <- Legend(legend_gp = gpar(fill = brewer.pal(9,'Set3')[-4]),nrow = 1,
                  title = "SNP Annotation", direction = "horizontal",
                  labels = c('intron',"5' UTR",'Others',"3' UTR",'upstream','intergenic','exonic'))#,#nr=1)



# combine heatmap

anno = anno_mark(at = which(bip_pops[rownames(target_isoform),'quantile']> 0.95), 
                 labels_gp = gpar(fontface='italic',fontsize=10),
                 labels = ensg_2genesymb[rownames(target_isoform),][which(bip_pops[rownames(target_isoform),'quantile']> 0.95)], 
                 which = "row")

pdf(file = "../figures/bip/bip_heatmap_finemap.pdf", width = 7, height = 8)

#draw(p1+p6+p6b+p3+p4_combined+p8_1+p8_2+p7+rowAnnotation(mark = anno),padding=unit(c(2.5,0,0,0),'cm'))

#draw(x = unit(0.4, "npc"), y = unit(0.07, "npc"),packLegend(lgd1,lgd6,lgd6b,lgd3,lgd4,lgd8_1,lgd7,lgd8_2,
                                                            #direction = "horizontal",max_width = unit(25, "cm")))
draw(p1+p6b+p8_1+p8_2+p7+rowAnnotation(mark = anno),padding=unit(c(2.5,0,0,0),'cm'))
draw(x = unit(0.4, "npc"), y = unit(0.07, "npc"),packLegend(lgd1,lgd6b,lgd8_1,lgd7,lgd8_2,
                                                            direction = "horizontal",max_width = unit(12, "cm")))

dev.off()



# # # # # # # # # # # # # # # # # # # # 
# # statistics
# # # # # # # # # # # # # # # # # # # # 

brain_target_5 <- list()
brain_target_n5 <- list()
bip_final_link$Region <- Br_color[bip_final_link$BR_Region,"Region"]

bip_final_link1 <- bip_final_link
bip_final_link1$Region <- ifelse(bip_final_link1$Region%in%c('NEX','Limbic'),'ForeBr',bip_final_link1$Region)

brain_target <- list()
for (region in unique(bip_final_link1$Region)){
  brain_target[[region]] <- unique(subset(bip_final_link1,Region==region)[,'TargetGene_symbol'])
}


p1 <- venn.diagram(brain_target, filename = NULL,na="remove",
                   height = 1000 ,  width = 1000 , resolution = 300,
                   col=c("#E31A1C","#33A02C", '#FDBF6F', '#FF7F00'),
                   fill = c(alpha("#E31A1C",0.3), alpha('#33A02C',0.3),
                            alpha('#FDBF6F',0.3), alpha('#FF7F00',0.3)),
)
pdf(file = "../figures/bip/venn.pdf", width = 3.5, height = 3.5)
plot_grid(p1)
dev.off()


## prop of  5', non 5' and both level for each region
isof <- c('non5'="Non 5'",'5'="5'",'both'='Both')

n_isoforms <- rbind(colSums(target_isoform==1),
                    colSums(target_isoform==2),
                    colSums(target_isoform==3))


colnames(n_isoforms) <- Br_color[colnames(n_isoforms) ,"subRG"]
rownames(n_isoforms) <- c('non5','5','both')
n_isoforms <- apply(n_isoforms, 2, function(x){x/sum(x)}) # for prop
n_isoforms <- melt(n_isoforms)
n_isoforms$Var1 <- factor(isof[n_isoforms$Var1], levels = c('Both',"5'","Non 5'"))

mean(subset(n_isoforms, Var1=="Non 5'")$value)
sd(subset(n_isoforms, Var1=="Non 5'")$value)

pdf(file = "../figures/bip/prop_iso.pdf", width = 4, height = 4)
ggplot(n_isoforms, aes(x=factor(Var2, levels = Br_color$subRG),y=value, fill=Var1))+
  geom_bar(stat = 'identity',width = 0.6)+
  scale_fill_manual('',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF"))+
  theme_classic()+xlab('')+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), legend.position = 'top',
        axis.line = element_line(colour = "black",linewidth = 0.2))+
  xlab('')+ylab('Proportion')# ('Number of target genes')
dev.off()



## pops score

target_pops <- c()

for (br in colnames(target_isoform)){
  targ_n5 <- rownames(target_isoform)[target_isoform[,br]==1]
  targ_5 <- rownames(target_isoform)[target_isoform[,br]==2]
  #targ_both <- rownames(target_isoform)[target_isoform[,br]==3]
  targ_backg <-   setdiff(closest_gene$closest_ensG, rownames(target_isoform)[target_isoform[,br]!=0])
  
  target_pops <- rbind(rbind(cbind(targ_n5,"Non 5'",br),
                             cbind(targ_5,"5'",br),
                             #cbind(targ_both,'Both',br),
                             cbind(targ_backg,'Closest gene',br)),target_pops)
}

target_pops <- as.data.frame(target_pops)
colnames(target_pops) <- c('gene','promoter','brain')
target_pops$brain <- Br_color[target_pops$brain,'subRG']
target_pops$PoPS <- bip_pops[target_pops$gene,'quantile']


#library(ggsignif)

pop_mean <- aggregate(data=target_pops,PoPS~brain+promoter,mean)
pop_sd <- aggregate(data=target_pops,PoPS~brain+promoter,sd)
pop_n <- aggregate(data=target_pops,PoPS~brain+promoter,length)
pop_mean$sd <- (pop_sd$PoPS)/sqrt(pop_n$PoPS)

pdf(file = "/sc/arion/projects/roussp01a/liting/Pf_25/PF_25BR/figures/bip/pops_updatedse.pdf.pdf", width = 4, height = 4)
ggplot(pop_mean, aes(x=factor(brain,levels = Br_color$subRG),color='black',
                     fill=factor(promoter, levels = c("5'","Non 5'",'Closest gene')),
                     y=PoPS))+
  geom_bar(stat = 'identity',position = 'dodge',width = 0.7)+
  geom_errorbar(aes(ymin=PoPS-sd, ymax=PoPS+sd,color="black"), width=.2,
                position=position_dodge(0.7)) + 
  scale_color_manual('',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF",'Closest gene'='orange'))+
  scale_fill_manual('',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF",'Closest gene'='orange'))+
  guides(color=guide_legend(override.aes = aes(label = "")))+
  theme_classic()+xlab('')+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), legend.position = 'top',
        axis.line = element_line(colour = "black",linewidth = 0.2))
dev.off()

write.table(pop_n, file='/sc/arion/projects/roussp01a/liting/Pf_25/PF_25BR/figures/bip/pop_n.txt',sep = '\t',quote = F)

wilcox.test(target_pops[target_pops$promoter=="Non 5'",'PoPS'], target_pops[target_pops$promoter=="5'",'PoPS'], alternative = 'greater')
wilcox.test(target_pops[target_pops$promoter=="Non 5'",'PoPS'], target_pops[target_pops$promoter=='Closest gene','PoPS'], alternative = 'greater')


# 
ggplot(target_pops, aes(x=brain,color=factor(promoter, levels = c("5'","Non 5'",'Closest gene')),
                        fill=factor(promoter, levels = c("5'","Non 5'",'Closest gene')),
                        y=PoPS))+
  geom_boxplot(outlier.size = 0.01,alpha=0.3) +
  scale_color_manual('Promoter isoform',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF",'Closest gene'='orange'))+
  scale_fill_manual('Promoter isoform',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF",'Closest gene'='orange'))+
  guides(color=guide_legend(override.aes = aes(label = "")))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), legend.position = '',
        axis.line = element_line(colour = "black",linewidth = 0.2))

#
ggplot(target_pops, aes(x=factor(promoter, levels = c("5'","Non 5'",'Closest gene')),
                        y=PoPS,fill=promoter))+
  facet_grid(~brain) +
  geom_boxplot(outlier.size = 0.1) +
  scale_color_manual('Promoter isoform',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF",'Closest gene'='orange'))+
  scale_fill_manual('Promoter isoform',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF",'Closest gene'='orange'))+
  guides(color=guide_legend(override.aes = aes(label = "")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), legend.position = '')+
  #stat_compare_means(label ="p.signif"  ) +
  geom_signif(
    comparisons = list(
      c("5'", "Non 5'"),
      #c("Non 5'", "Both"),
      c("Non 5'", "Closest gene")
    ),#,#map_signif_level = function(p) sprintf("%.2g", p)
    map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "ns."=2))+
  xlab('')+ylab('PoPS score')+ylim(c(0,1.2))

dev.off()





























