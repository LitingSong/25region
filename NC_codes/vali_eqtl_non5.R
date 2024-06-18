# validation the ABC links 
# 5' vs non 5'
# enhancer peak length: 1000

library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(reshape2)


setwd('/hpc/users/songl05/PF_25BR/')

# validate the ABC links using biccn brain data
options(stringsAsFactors = F)
load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata')

epLinks <- list()
for (brain_Rg in names(GRs)){
  grs <- as.data.frame(GRs[[brain_Rg]])
  #grs$chr <- grs$seqnames
  grs$start <- merged_peak_500[grs$eID,]@ranges@start - 250
  grs$end <- grs$start + 1000
  #grs$TargetGene_symbol <- ensg_2genesymb[grs$TargetGene,]
  epLinks[[brain_Rg]] <- makeGRangesFromDataFrame(grs[,-4:-6], keep.extra.columns = T)
}

GRs <- epLinks

Br_color <- read.table('./data/color_br.txt',sep='\t',comment.char = '', header = T, row.names = 1)

# get region-specific IDS (specific_link)

ABC <- c()

for( br in names(GRs)){
  raw_score <- as.data.frame(GRs[[br]][,c('ID','eID','ABC.Score','CellType','order','TargetGene')])#[,-1:-5]
  raw_score <- raw_score[order(raw_score$order),]
  raw_score <- raw_score[!duplicated(raw_score$ID),]
 
  ABC <- rbind(ABC,raw_score)
  
}

ABC_5 <- subset(ABC, order==1)
ABC_n5 <- subset(ABC, order!=1)

get_specific_link <- function(type){
  if(type=='5'){ABC <- ABC_5}
  if(type=='n5'){ABC <- ABC_n5}
  
  ABC_m <- dcast(ABC,ID ~ CellType, value.var = "ABC.Score")
  rownames(ABC_m) <- ABC_m$ID
  ABC_m <- ABC_m[,-1]
  ABC_m <- ABC_m[ , names(GRs)]
  
  specific_link <- vector(mode='list', length=ncol(ABC_m))
  for (i in 1:ncol(ABC_m)){
    specific_link[[i]] = rownames(ABC_m)[!is.na(ABC_m[,i]) & rowSums(is.na(ABC_m[,-i]))==(ncol(ABC_m)-1)]
  }
  
  return(specific_link)
  
}

fishe_test_spec <- function(prop_validated,i,j){
  
  prop_df_vali <- matrix(prop_validated$V1,i,j, byrow = T)
  prop_df_over <- matrix(prop_validated$V2,i,j, byrow = T)
  
  n_vali = prop_df_vali 
  n_other_vali=rowSums(prop_df_vali)-prop_df_vali
  n_over=prop_df_over
  n_other_over=rowSums(prop_df_over)-prop_df_over
  
  
  
  specif_matrix <- data.frame(n_vali=as.vector(n_vali),
                              n_other_vali=as.vector(n_other_vali),
                              n_over=as.vector(n_over),
                              n_other_over=as.vector(n_other_over))
  
  p_val <- apply(specif_matrix, 1, function(x){fisher.test( matrix(x,nrow=2) , alternative='greater')$p.value})
  or_val <- apply(specif_matrix, 1, function(x){fisher.test( matrix(x,nrow=2) , alternative='greater')$estimate})
  
  pv_matrix <- matrix(p_val, i,j)
  or_matrix <- matrix(or_val, i,j)
  
  return(list(pv_matrix=pv_matrix,or_matrix=or_matrix))
}

# using BICCN data
BICCN <- function(type ){
  
  if(type=='5'){ABC_d <- ABC_5}
  if(type=='n5'){ABC_d <- ABC_n5}
  

  biccn_brain <- dir(path='/sc/arion/projects/roussp01a/liting/Pf_25/biccn',pattern = '*pos.bedpe',full.names = T)
  biccn_brain <- biccn_brain[ !grepl('ACBGM|BFEXA|MGC|SIGA|THMGA|merged|OPC|OGC|ASCT|CBGRC|PV_Ch|NP' ,biccn_brain) ]
  biccn_brs <- str_split(biccn_brain, '\\/|\\.',simplify = T)[,9]
  
  Br_color <- read.table('./data/color_br.txt',sep='\t',comment.char = '', header = T, row.names = 1)
  
  ensg_2genesymb <- read.csv('/sc/arion/projects/roussp01a/liting/Pf_25/ensg_2genesymb.csv', row.names = 1)
  
  prop_validated <- c()
  for (i in 1:length(biccn_brs)){
    
    biccn_link <- read.table(biccn_brain[i],sep='\t')
    colnames(biccn_link)[1:7] <- c("chr", "start_p", "end_p", "chr_d",	"start",	"end","gene" )
    
    gr0 <- as(biccn_link, "GRanges")
    
    specific_link <- get_specific_link(type=type)
    for (j in 1:length(GRs)){
      
      ABC <- subset(ABC_d, CellType==names(GRs)[j] & ID%in%specific_link[[j]])
      #gr1 <- GRs[[j]][!duplicated(GRs[[j]]$ID)]
      #gr1 <-  subset(gr1, ID%in%specific_link[[j]])
      #ABC <- as.data.frame(gr1)
      gr1 <- as(ABC, "GRanges")
      
      hits = findOverlaps(gr0, gr1)
      
      overlap_abc_biccn <- cbind(as.data.frame(gr0)[queryHits(hits),], 
                                 ABC[subjectHits(hits), c('TargetGene','ID') ])
      overlap_abc_biccn$gene <- str_split(overlap_abc_biccn$gene,'\\|',simplify = T)[,1]
      overlap_abc_biccn$gene_sym <- ensg_2genesymb[overlap_abc_biccn$TargetGene,]
      
      n_vali <- length(unique(overlap_abc_biccn$ID[overlap_abc_biccn$gene==overlap_abc_biccn$gene_sym]))
      n_over <- length(unique(overlap_abc_biccn$ID)) 
      
      prop_validated <- rbind(prop_validated,c(n_vali,n_over ))
      
    }
  }
  
  #prop_validated <- as.data.frame(prop_validated)%>%mutate(prop=V1/V2)
  prop_validated <- as.data.frame(prop_validated)%>%mutate(prop=ifelse(V2>=50,V1/V2,NA))
  
  rownames(prop_validated) <- (as.data.frame(expand.grid(i=names(GRs),j=biccn_brs))%>%mutate(id=paste(i,j,sep='-')))[,'id']
  #prop_validated$prop[prop_validated$V2< 100] <- NA
  
  
  pv_df <- fishe_test_spec(prop_validated,i,j)$pv_matrix
  rownames(pv_df) <- biccn_brs
  colnames(pv_df) <- names(GRs)
  
  or_df <- fishe_test_spec(prop_validated,i,j)$or_matrix
  rownames(or_df) <- biccn_brs
  colnames(or_df) <- names(GRs)
  
  prop_df <- matrix(prop_validated$prop,length(biccn_brs),length(GRs), byrow = T)
  
  rownames(prop_df) <- biccn_brs
  colnames(prop_df) <- names(GRs)
  prop_df <- prop_df*100
  
  
  
  return(list(prop_df = prop_df, pv_df = pv_df, or_df=or_df))
}

# using eqtl data
Eqtl <- function(type){
  
  if(type=='5'){ABC_d <- ABC_5}
  if(type=='n5'){ABC_d <- ABC_n5}
  
  
  eqtl_dir <- c('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc/fine/GTex/GTEx_vcf/') 
  eqtl_brain <- dir(path=eqtl_dir, pattern = '.*_eqtls.txt.gz',full.names = F)
  
  eqtl_brain <- eqtl_brain[ !grepl('Cerebell|Spinal' ,eqtl_brain) ]
  
  plink_ld.8 <- read.table('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc/fine/GTex/GTEx_vcf/plink.ld.0.8',head=T)
  
  prop_validated <- c()
  for (i in 1:length(eqtl_brain)){
    print(i)
    eqtl_sig <- read.table(gzfile(paste0( eqtl_dir, eqtl_brain[i])), header = T)[,c('gene_id','variant_id')]
    eqtl_sig$CHR_A <- gsub('chr','',str_split(eqtl_sig$variant_id,'_', simplify = T)[,1])
    eqtl_sig$BP_A <- as.numeric(str_split(eqtl_sig$variant_id,'_', simplify = T)[,2])
    
    eqtl_sig_ld <- merge(eqtl_sig, plink_ld.8, by=c('CHR_A','BP_A'))
    
    ld_target <- eqtl_sig_ld[,c('CHR_B','BP_B','gene_id')]
    colnames(ld_target) <- c('CHR_A','BP_A','gene_id')
    
    eqtl_ld <- unique(rbind(ld_target, eqtl_sig[,c('CHR_A','BP_A','gene_id')]))
    eqtl_ld <-within(eqtl_ld, {
      chr=paste('chr',CHR_A,sep='')
      start=BP_A
      stop=BP_A
      gene_id = str_split(gene_id,'\\.',2,simplify = T)[,1] })
    
    eqtl_ld <- eqtl_ld[,c("chr","start","stop", "gene_id")]
    
    
    gr0 = with(eqtl_ld, GRanges(chr, IRanges(start=start, end=stop)))
    specific_link <- get_specific_link(type=type)
    
    for (j in 1:length(GRs)){
      print(j)
      ABC <- subset(ABC_d, CellType==names(GRs)[j] & ID%in%specific_link[[j]])
      #gr1 <- GRs[[j]][!duplicated(GRs[[j]]$ID)]
      #gr1 <-  subset(gr1, ID%in%specific_link[[j]])
      #ABC <- as.data.frame(gr1)
      gr1 <- as(ABC, "GRanges")
      
      # hits = findOverlaps(gr0, gr1)
      # gr1 <- GRs[[j]][!duplicated(GRs[[j]]$ID)]
      # gr1 <-  subset(gr1, ID%in%specific_link[[j]])
      # ABC <- as.data.frame(gr1)
      hits = findOverlaps(gr0, gr1)
      
      overlap_abc_eqtl <- cbind(eqtl_ld[queryHits(hits),], ABC[subjectHits(hits), c('TargetGene','ID') ])
      
      #
      n_vali <- length(unique(overlap_abc_eqtl$ID[overlap_abc_eqtl$gene_id==overlap_abc_eqtl$TargetGene]))
      n_over <- length(unique(overlap_abc_eqtl$ID)) 
      
      
      prop_validated <- rbind(prop_validated,c(n_vali,n_over ))
    }
  }
  
  #prop_validated <- as.data.frame(prop_validated)%>%mutate(prop=V1/V2)
  prop_validated <- as.data.frame(prop_validated)%>%mutate(prop=ifelse(V2>=50,V1/V2,NA))
  
  eqtl_region <- str_split(gsub('Brain_','',eqtl_brain),'\\.',simplify = T)[,1]
  
  pv_df <- fishe_test_spec(prop_validated,i,j)$pv_matrix
  rownames(pv_df) <- eqtl_region
  colnames(pv_df) <- names(GRs)
  
  or_df <- fishe_test_spec(prop_validated,i,j)$or_matrix
  rownames(or_df) <- eqtl_region
  colnames(or_df) <- names(GRs)
  
  #prop_df <- matrix(prop_validated$V2,length(eqtl_brain),length(GRs), byrow = T)
  prop_df <- matrix(prop_validated$prop,length(eqtl_brain),length(GRs), byrow = T)
  
  
  
  rownames(prop_df) <- eqtl_region
  colnames(prop_df) <- names(GRs)
  prop_df <- prop_df*100
  return(list(prop_df=prop_df,pv_df=pv_df, or_df=or_df))
}


# run validation
biccn_v_5 <- BICCN('5')
eqtl_v_5 <- Eqtl('5')

biccn_v_n5 <- BICCN('n5')
eqtl_v_n5 <- Eqtl('n5')

# 5'
prop_df_biccn_5 <- biccn_v_5[['prop_df']]
prop_df_eqtl_5 <- eqtl_v_5[['prop_df']]

or_df_biccn_5 <- biccn_v_5[['or_df']]
or_df_eqtl_5 <- biccn_v_5[['or_df']]

pv_df_biccn_5 <- biccn_v_5[['pv_df']]
pv_df_eqtl_5 <- eqtl_v_5[['pv_df']]

# non 5'

prop_df_biccn_n5 <- biccn_v_n5[['prop_df']]
prop_df_eqtl_n5 <- eqtl_v_n5[['prop_df']]

or_df_biccn_n5 <- biccn_v_n5[['or_df']]
or_df_eqtl_n5 <- biccn_v_n5[['or_df']]

pv_df_biccn_n5 <- biccn_v_n5[['pv_df']]
pv_df_eqtl_n5 <- eqtl_v_n5[['pv_df']]


save(prop_df_biccn,prop_df_eqtl,prop_df_epimap, file='./data/prop_df_biccn.RData')

# heatmap  
plot_heat <- function(or_df,pv_df, legend_name,column_title){
  
  col_fun <- colorRamp2(seq(max(or_df,na.rm = T),min(or_df,na.rm = T),length=8), (brewer.pal(8, "RdBu")))
  
  Heatmap(or_df,col_fun ,name = legend_name,
          column_labels = Br_color[colnames(or_df),'subRG'],#cluster_columns  = F,
          column_title = column_title,
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(pv_df[i, j] <0.05 )
              grid.text(sprintf("%1.1e", pv_df[i, j]), x, y, gp = gpar(fontsize = 8))},
          
          bottom_annotation = HeatmapAnnotation(
            Region= Br_color[colnames(or_df),'Region'],
            BroadRegion = Br_color[colnames(or_df),'Group'],
            col = list(Region = c("Dien" = "#FDBF6F", 
                                  "BasGan" = "#33A02C", 
                                  "MidBr" = "#FF7F00",
                                  'Limbic'='#FB9A99',
                                  "NEX"='#E31A1C'),
                       BroadRegion=c("MidDien"="#FDBF6F",
                               "BasGan"="#33A02C",
                               "ForeBr"="#FB9A99"))
            
          ) ) }

## original value
## BICCN
pdf(file = "./figures/sp2_abc_biccn_valid_heatmap_eid.pdf", width = 14, height = 7)
five <- plot_heat(or_df=prop_df_biccn_5,pv_df=pv_df_biccn_5,legend_name='Prop', column_title= "5'")
n5 <- plot_heat(or_df=prop_df_biccn_n5,pv_df=pv_df_biccn_n5,legend_name='Prop',  column_title= "non-5'")
draw(five+n5)
#plot_heat(prop_df=prop_df_biccn,legend_name='Proportion (%)')
dev.off()

## EQTL
pdf(file = "./figures/sp2_abc_eqtl_valid_heatmap_5_pair.pdf", width = 14, height = 7)
five <-plot_heat(or_df=prop_df_eqtl_5,pv_df=pv_df_eqtl_5,legend_name='Prop', column_title= "5'")
n5 <-plot_heat(or_df=prop_df_eqtl_n5,pv_df=pv_df_eqtl_n5,legend_name='Prop',  column_title= "non-5'")
draw(five+n5)
#plot_heat(prop_df=prop_df_eqtl,legend_name='Proportion (%)')
dev.off()



# barplot to compare the validated proportion for 5' non5' links
prop_eqtl <- cbind(melt(prop_df_eqtl_5 ), melt(prop_df_eqtl_n5 )[,-1:-2])
colnames(prop_eqtl) <- c("eqtl","abc","5'","non5'")
prop_eqtl$pair <- paste(prop_eqtl$abc,prop_eqtl$eqtl,sep=' -- ')

prop_eqtl <- subset(prop_eqtl,pair%in%c("Neocortex -- Frontal_Cortex_BA9",
                           "Neocortex -- Cortex"  ,
                           "Limbic -- Hippocampus",
                          
                           "BasalGanglia -- Caudate_basal_ganglia",
                           "BasalGanglia -- Nucleus_accumbens_basal_ganglia",
                           "BasalGanglia -- Putamen_basal_ganglia",
                           "ARC -- Hypothalamus"
                           )) #"Limbic -- Amygdala" ,


library(ggplot2)
library(reshape2)
prop_eqtl <- melt(prop_eqtl[,-1:-2])

prop_eqtl_pair_bar <- ggplot(prop_eqtl, aes(x=pair,y=value,fill=variable))+geom_bar(stat='identity', position='dodge', width=0.7)+
  coord_flip()+ theme_bw() + ylab('Validated Proportion') + xlab('')+
  scale_fill_manual('Links',values= c("5'"="#FEE090FF" , "non5'"="#4393C3FF"  ))
  
  

ggsave(prop_eqtl_pair_bar, file='/hpc/users/songl05/PF_25BR/figures/sp2_prop_eqtl_pair_bar.pdf', width=6, height=4)







# 
# 
# ## scaled value
# pdf(file = "./figures/abc_biccn_valid_heatmap_scale_prox.pdf", width = 9, height = 7)
# prop_df_scale <- t(apply(prop_df_biccn, 1,scale ))
# colnames(prop_df_scale) <- colnames(prop_df_biccn)
# plot_heat(prop_df=prop_df_scale,pv_df=pv_df_biccn,legend_name='Scaled proportion')
# dev.off()
# 
# pdf(file = "./figures/abc_eqtl_valid_heatmap_scale_prox.pdf", width = 9, height = 7)
# prop_df_scale <- t(apply(prop_df_eqtl, 1,scale ))
# colnames(prop_df_scale) <- colnames(prop_df_eqtl)
# plot_heat(prop_df=prop_df_scale,pv_df=pv_df_eqtl,legend_name='Scaled proportion')
# dev.off()
# 
# pdf(file = "./figures/abc_epimap_valid_heatmap_scale_prox.pdf", width = 9, height = 7)
# prop_df_scale <- t(apply(prop_df_epimap, 1,scale ))
# colnames(prop_df_scale) <- colnames(prop_df_epimap)
# plot_heat(prop_df=prop_df_scale,pv_df=pv_df_epimap,legend_name='Scaled proportion')
# dev.off()
# 
# 
# 
# 
# 
# 
# 
