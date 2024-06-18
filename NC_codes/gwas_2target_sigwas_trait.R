# significant gwas snp to target gene across neuropsychiatric disorders

library(stringr)
library(reshape2)
library(rtracklayer)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(cowplot)
library(data.table)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# part1: get targets
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

setwd('~/PF_25BR/data/')

#source('../codes/get_sig_gwas.R') # significant snp
load('gsums_sig.RData')
load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata')

ensg_2genesymb <- read.csv('/sc/arion/projects/roussp01a/liting/Pf_25/ensg_2genesymb.csv', row.names = 1)
Br_color <- read.table('color_br.txt',sep='\t',comment.char = '', header = T, row.names = 1)

gsums_sig <- plyr::rename(gsums_sig,c("CHR"="chr", "BP"='start'))
gsums_sig$chr <- paste0('chr',gsums_sig$chr)
gsums_sig$end <- gsums_sig$start

gsums_sig <- makeGRangesFromDataFrame(gsums_sig, keep.extra.columns = T)

ch = import.chain("hg19ToHg38.over.chain")

gsums_sig_38 <- liftOver(gsums_sig,ch)


# for annotation and find closest gene
write.table(unique(as.data.frame(gsums_sig_38)[, c("seqnames", "start", "end" ,"A1",'A2','SNP' )]),
            sep='\t', file="gwas_sig.avinput", quote = F,row.names = F,col.names = F)

gsums_sig_38_trait <- as.data.frame(gsums_sig_38)[, c("seqnames", "start", "end" ,"A1",'A2','SNP','trait' )]

system('sh /hpc/users/songl05/PF_25BR/codes/annovar.sh')
system('sh /hpc/users/songl05/PF_25BR/codes/find.closest.sh')

anno_hg38 <- read.table("gwas_sig.anno.hg38_multianno.txt", sep='\t',header = T )

closest_gene <- read.table("gwas_sig.answer.bed", sep='\t' )

closest_gene$closest_G <- str_split(closest_gene$V15,'\\s|;',simplify = T)[,8]
closest_gene$closest_ensG <- str_split(closest_gene$V15,'\\s|;|\\.',simplify = T)[,2]
closest_gene$closest_distance <- str_split(closest_gene$V15,'\\|',simplify = T)[,2]
closest_gene$index <- paste(closest_gene$V1,closest_gene$V2)

# target gene/isoform

epLinks <- list()
dista <- 1000
for (brain_Rg in names(GRs)){
  grs <- as.data.frame(GRs[[brain_Rg]])
  #grs$chr <- grs$seqnames
  grs$start <- merged_peak_500[grs$eID,]@ranges@start - ((dista/2)-250)
  grs$end <- grs$start + dista
  grs$TargetGene_symbol <- ensg_2genesymb[grs$TargetGene,]
  epLinks[[brain_Rg]] <- makeGRangesFromDataFrame(grs[,-4:-6], keep.extra.columns = T)
}

sig_link <- list()
for(br in names(epLinks)){
  hits = findOverlaps(epLinks[[br]], gsums_sig_38)
  results = cbind.data.frame(epLinks[[br]][hits@from,], gsums_sig_38[hits@to,])
  results$BR_Region <- br
  sig_link[[br]] <- results
}


final_link = data.frame(do.call("rbind.data.frame", sig_link))
final_link <- final_link%>%group_by(SNP,BR_Region,trait)%>%top_n(1,ABC.Score)
final_link <- as.data.frame(rename(final_link, c("gene"="TargetGene")))


final_link <- left_join(final_link,closest_gene[,c('V1','V2','closest_G','closest_ensG','closest_distance')], c("seqnames.1" = "V1","start.1"="V2"))
final_link <- left_join(final_link,anno_hg38[,c('Chr','Start','Func.refGene','Gene.refGene')], c("seqnames.1" = "Chr","start.1"="Start"))

save(final_link, file='siggwas_final_link.RData')


final_link_supple <- unique(final_link[,c('SNP','P','trait','gene','TargetGene_symbol','TargetTranscript','order','BR_Region')])
final_link_supple <- final_link_supple[order(final_link_supple$trait),]

write.table(final_link_supple, file='/sc/arion/projects/roussp01a/liting/Pf_25/data/supple_tables/prioritized_genes_traits_supple.txt',
            quote = F,row.names = F,sep='\t')
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # #    part2: statistic   functions         # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

## get 5', non 5' and both target matrix
get_isoform_matrix <- function(Trait){
  
  trait_final_link <- subset(final_link,trait==Trait  )
  
  target_5 <- subset(trait_final_link, order==1)
  target_n5 <-subset(trait_final_link, order!=1)
  
  
  target_n5_gene <- dcast(data =target_n5, formula = gene~ BR_Region,length) 
  target_n5_gene[,-1] <- ifelse(target_n5_gene[,-1]==0,0,1) # non5'
  target_5_gene <- dcast(data =target_5, formula = gene~ BR_Region,length)
  target_5_gene[,-1] <- ifelse(target_5_gene[,-1]==0,0,2) # 5'
  
  print(nrow(target_n5_gene))
  print(nrow(target_5_gene))
  
  #target_gene <- left_join(target_n5_gene, target_5_gene,by='gene')
  target_gene <- merge(target_n5_gene, target_5_gene,by='gene',all.x=T,all.y=T)
  
  rownames(target_gene) <- target_gene$gene
  target_gene <- target_gene[,-1]
  target_gene[is.na(target_gene)] <- 0
  target_isoform <- target_gene[,c(1:9)] + target_gene[,c(10:18)] # 0: none ;1: non5'; 2:5'; 3: both;
  colnames(target_isoform) <- gsub('\\.x','',colnames(target_isoform))
  
  return(target_isoform)
  
}

## abc score
get_abc <- function(tt){
  
  link_abc <- subset(final_link,trait==tt  )[,c('gene','order','ABC.Score','BR_Region')]

  ep_abc_mean <- aggregate(ABC.Score~ gene + BR_Region ,link_abc,mean)
  ep_abc_mean <- reshape2::dcast(ep_abc_mean, gene~BR_Region, value=ABC.Score)
  rownames(ep_abc_mean) <- ep_abc_mean$gene
  ep_abc_mean <- ep_abc_mean[,-1]
  
  target_isoform <- get_isoform_matrix(tt)
    
  ep_abc_mean <- ep_abc_mean[rownames(target_isoform),]
  
  
  ep_abcSc <- cbind(melt(ep_abc_mean), melt(target_isoform))[,-3]
  colnames(ep_abcSc) <- c('brain','abc','promoter')
  ep_abcSc <- subset(ep_abcSc, promoter%in%c(1,2,3)) 
  ep_abcSc$promoter <- ifelse(ep_abcSc$promoter==1,"Non 5'",ifelse(ep_abcSc$promoter==2,"5'",'Both'))
  ep_abcSc$brain <- Br_color[ep_abcSc$brain,'subRG']
  
  ep_abcSc <- cbind(melt(ep_abc_mean), melt(target_isoform))[,-3]
  colnames(ep_abcSc) <- c('brain','abc','promoter')
  ep_abcSc <- subset(ep_abcSc, promoter%in%c(1,2,3)) 
  ep_abcSc$promoter <- ifelse(ep_abcSc$promoter==1,"Non 5'",ifelse(ep_abcSc$promoter==2,"5'",'Both'))
  ep_abcSc$brain <- Br_color[ep_abcSc$brain,'subRG']
  ep_abcSc$trait <- tt
  return(ep_abcSc)
}

## e-p distance
get_distance <- function(tt){

  link_ep <- subset(final_link,trait==tt  )[,c('gene','TargetTranscript','eID','order','BR_Region','start.1')]
  link_ep$prom_pos <- promoter_anno[link_ep$TargetTranscript,]@ranges@start
  link_ep$enhancer_pos <- merged_peak_500[link_ep$eID,]@ranges@start
  link_ep$ep_dist <- abs(link_ep$prom_pos - link_ep$enhancer_pos )
  link_ep$vp_dist <- abs(link_ep$start.1 - link_ep$prom_pos )
  
  vp_dist_max <- aggregate(vp_dist~ gene + BR_Region ,link_ep,mean)
  vp_dist_max <- reshape2::dcast(vp_dist_max, gene~BR_Region, value=vp_dist)
  
  target_isoform <- get_isoform_matrix(tt)
  
  rownames(vp_dist_max) <- vp_dist_max$gene
  vp_dist_max <- vp_dist_max[,-1]
  vp_dist_max <- vp_dist_max[rownames(target_isoform),]
  
  vp_dist <- cbind(melt(vp_dist_max), melt(target_isoform))[,-3]
  colnames(vp_dist) <- c('brain','dis','promoter')
  vp_dist <- subset(vp_dist, promoter%in%c(1,2,3)) 
  vp_dist$promoter <- ifelse(vp_dist$promoter==1,"Non 5'",ifelse(vp_dist$promoter==2,"5'",'Both'))
  vp_dist$brain <- Br_color[vp_dist$brain,'subRG']
  vp_dist$trait <- tt
  
  
  return(vp_dist)
}


## 5 and non 5 proportion
get_n5_prop <- function(tt){
  
  isof <- c('non5'="Non 5'",'5'="5'",'both'='Both')
  
  n_isoforms <- rbind(colSums(target_isoform[[tt]]==1),
                      colSums(target_isoform[[tt]]==2),
                      colSums(target_isoform[[tt]]==3))
  
  #n_isoforms <- apply(target_isoform[[tt]], 2, table)[-1,]
  colnames(n_isoforms) <- Br_color[colnames(n_isoforms) ,"subRG"]
  rownames(n_isoforms) <- c('non5','5','both')
  p_isoforms <- apply(n_isoforms, 2, function(x){x/sum(x)}) # for prop
  p_isoforms <- melt(p_isoforms)
  p_isoforms$Var1 <- factor(isof[p_isoforms$Var1], levels = c('Both',"5'","Non 5'"))
  p_isoforms$trait <- tt
  p_isoforms$n <- melt(n_isoforms)$value
  return(p_isoforms)
  
}

## get pops score
get_pops <- function(tt){
  
  
  TT_pops <- read.table(paste0('added_ridge_',tt,'.preds'), sep='\t', header = T, row.names =1 )
  TT_pops$quantile <- rank(TT_pops$PoPS_Score)/length(TT_pops$PoPS_Score) 
  
  #gsums_sig_tt <- as.data.frame(gsums_sig_38)%>%subset(trait==tt)
  gsums_sig_tt <- get_lead(tt)
  

  target_pops <- c()
  target_isoform_t <- target_isoform[[tt]]
  for (br in colnames(target_isoform_t)){
    targ_n5 <- rownames(target_isoform_t)[target_isoform_t[,br]==1]
    targ_5 <- rownames(target_isoform_t)[target_isoform_t[,br]==2]
    #targ_both <- rownames(target_isoform_t)[target_isoform_t[,br]==3]
    targ_backg <-   setdiff(unique(left_join(gsums_sig_tt,closest_gene,c("seqnames" = "V1","start"="V2"))[,'closest_ensG']),
                          rownames(target_isoform_t)[target_isoform_t[,br]!=0])
    
    
    target_pops <- rbind(rbind(cbind(targ_n5,"Non 5'",br),
                               cbind(targ_5,"5'",br),
                               #cbind(targ_both,'Both',br),
                               cbind(targ_backg,'Closest gene',br)),target_pops)
  }
  
  target_pops <- as.data.frame(target_pops)
  colnames(target_pops) <- c('gene','promoter','brain')
  target_pops$brain <- Br_color[target_pops$brain,'subRG']
  target_pops$trait <- Trait
  target_pops$PoPS <- TT_pops[target_pops$gene,'quantile']
  
  return(target_pops)
  

}

## plot venn
plot_venn <- function(tt){
  
  
  tt_final_link  <- subset(final_link,trait==tt  )
  
  tt_final_link$Region <- Br_color[tt_final_link$BR_Region,"Region"]
  
  
  tt_final_link$Region <- ifelse(tt_final_link$Region%in%c('NEX','Limbic'),'ForeBr',tt_final_link$Region)
  
  brain_target <- list()
  for (region in unique(tt_final_link$Region)){
    brain_target[[region]] <- unique(subset(tt_final_link,Region==region)[,'gene'])
  }
  
  
  p1 <- venn.diagram(brain_target, filename=NULL,
                     sub=tt,
                     #filename = paste0('../figures/combined_trait/p2.',tt,'_venn_target.png'),
                     na="remove",
               height = 1000 ,  width = 1000 , resolution = 300,
               cex=0.7,cat.cex = 0.7,
               col=c("#E31A1C","#33A02C", '#FDBF6F', '#FF7F00'),
               fill = c(alpha("#E31A1C",0.3), alpha('#33A02C',0.3),alpha('#FDBF6F',0.3), alpha('#FF7F00',0.3)),
  )
  return(p1)
}  


## get lead snp
library(GenomicRanges)
get_lead <- function(tt){
  df <- subset(as.data.frame(gsums_sig_38),trait==tt)
   
  # make a genomic range object
  gr = GRanges(seqnames=df$seqnames,
               IRanges(start=df$start,end=df$start),P=df$P)
  names(gr) = df$SNP
  
  # you can change this
  FLANK = 10000 
  
  REGIONS <- reduce(flank(gr,FLANK,both=TRUE))
  # each SNP can only be matched to one merged region
  # so we just find overlap between region and snp
  # and assign the snp to the region
  gr$region = subjectHits(findOverlaps(gr,REGIONS))
  # order by pvalue
  gr = gr[order(gr$P),]
  # keep only the top snp in each region
  gr <- as.data.frame(gr[!duplicated(gr$region)])
  gr$index <- paste(gr$seqnames,gr$start)
  return(gr)
}

## eqtl
get_eqtl <- function(tt){ 
  
  sig_38 <- subset(gsums_sig_38_trait,trait==tt)
  sig_38$index <- paste(sig_38$seqnames, sig_38$start, sep=":")
  
  #eqtl_dir <- './'
  eqtl_dir <- '/hpc/users/songl05/PF_25BR/data/GTEx_Analysis_v8_eQTL/'
  
  eqtl_brain <- dir(path=eqtl_dir, pattern = 'Brain.*signif_variant_gene_pairs.txt.gz',full.names = F)
  #eqtl_brain <- eqtl_brain[ !grepl('Cerebell|Spinal' ,eqtl_brain) ]
  
  eqtl_target  <- c()
  for (i in 1: length(eqtl_brain)){
    eqtl_sig <- read.delim(gzfile(paste0( eqtl_dir, eqtl_brain[i])), header = T)[,c('gene_id','variant_id')]
    eqtl_sig[,c('chr','pos')] <- str_split(eqtl_sig$variant_id, '_',simplify = T)[,1:2]
    
    eqtl_sig$index <- paste(eqtl_sig$chr, eqtl_sig$pos, sep=":")
    
    e_gene <- unique(subset(eqtl_sig, index%in%sig_38$index)[,'gene_id'])
    eqtl_target <- c(eqtl_target, e_gene)
  }
  eqtl_target <- str_split_fixed(eqtl_target,'\\.',2)[,1]
  target_isoform_t <- target_isoform[[tt]]
  prop <- length(intersect(rownames(target_isoform_t),eqtl_target))/nrow(target_isoform_t)
  
  return(prop)
 
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # 
##    part3:  statistic main programs
# # # # # # # # # # # # # # # # # # # # # # # # # # # # 

target_isoform <- list()
for (tt in setdiff(unique(final_link$trait),c("insomn2",'adhd'))){
  print(tt)
  target_isoform[[tt]] <- get_isoform_matrix(tt)
}

# p1 number and proportion of  5 and non 5
trait_info <- read.table('./trait.txt',row.names = 1,sep='\t')
n5_prop <- list()
for (tt in setdiff(unique(final_link$trait),c("insomn2",'adhd'))){
  print(tt)
  n5_prop[[tt]] <- get_n5_prop(tt)
}
n5_prop_com <- do.call('rbind',n5_prop)

pdf(file = "../figures/combined_trait/p1.prop_iso.pdf", width = 8, height = 4.5)
n5_prop_com$trait <- trait_info[n5_prop_com$trait,]
ggplot(n5_prop_com, aes(x=factor(Var2, levels = Br_color$subRG),y=value, fill=Var1))+
  geom_bar(stat = 'identity',width = 0.6)+
  scale_fill_manual('Promoter isoform',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF"))+
  theme_bw()+ facet_grid(.~trait)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), legend.position = 'top')+
  xlab('')+ylab('Proportion')# ('Number of target genes')
dev.off()

mean(subset(n5_prop_com, Var1=="Non 5'")$value)
(sd(subset(n5_prop_com, Var1=="Non 5'")$value))/sqrt(54)

pdf(file = "../figures/combined_trait/p1.n_iso.pdf", width = 8, height = 4.5)
ggplot(n5_prop_com, aes(x=factor(Var2, levels = Br_color$subRG),y=n, fill=Var1))+
  geom_bar(stat = 'identity',width = 0.6)+
  scale_fill_manual('Promoter isoform',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF"))+
  theme_bw()+ facet_wrap(~trait,nrow = 1)+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), legend.position = 'top')+
  xlab('')+ylab( 'Number of target genes')
dev.off()

## xxxx
targets <- unique(final_link[,c('gene','CellType','trait')])
targets_table <- aggregate(data=targets, CellType~gene+trait, length)
targets_count <- aggregate(data=targets_table, .~CellType+trait, length)
targets_total <- aggregate(data=targets_count, gene~trait, sum)

target_sta <- merge(targets_count,targets_total,by="trait")%>%mutate(prop=gene.x/gene.y)%>%subset(CellType==1 &!trait%in%c("insomn2",'adhd') )
mean(target_sta$prop)
sd(target_sta$prop)/sqrt(6*9)

library(ggsci)
pdf(file = "../figures/supp/sp_traits_prop_ntargets.pdf", width = 5, height = 5)
targets_table$trait <- trait_info[targets_table$trait,]
ggplot(subset(targets_table,!trait%in%c("insomn2",'adhd')), aes(x=trait, fill=as.character(CellType)))+
  geom_bar(position = 'fill',width=0.7)+theme_bw() +
  scale_fill_jco() +
  xlab('disorder') +ylab('Proportion of the number of involved regions')+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), 
        axis.line = element_line(colour = "black",linewidth = 0.2))+
  guides(fill=guide_legend(title="Number of regions")) 
  
dev.off()

# p2 venn diagrame (overlap among brain regions)
pdf(file = "../figures/combined_trait/p2.venn.pdf", width = 7.5, height = 5)
p2_v1 <- plot_venn("sz3")
p2_v2 <- plot_venn("bip2")
p2_v3 <- plot_venn("neu2")
p2_v4 <- plot_venn("smoking")
p2_v5 <- plot_venn("drinking")
p2_v6 <- plot_venn("mdd_wray")
plot_grid(p2_v1,p2_v2,p2_v3,p2_v4,p2_v5,p2_v6)
dev.off()


## p3 abc score distribution

ep_abcSc <- list()
for (Trait in setdiff(unique(final_link$trait),c("insomn2",'adhd')) ){
  print(Trait)
  ep_abcSc[[Trait]] <-  get_abc(Trait)
}
ep_abcSc <- as.data.frame(do.call(rbind,ep_abcSc))


pdf(file = "../figures/combined_trait/p3.abcSc.pdf", width = 7.5, height = 5)
ep_abcSc$trait <- trait_info[ep_abcSc$trait,]
ggplot(subset(ep_abcSc,promoter!='Both'), aes(x=factor(brain, levels = Br_color$subRG),y=abc, fill=factor(promoter, levels = c("5'","Non 5'",'Both'))))+
  geom_boxplot(outlier.size = 0.1)+
  scale_fill_manual('Promoter isoform',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF"))+
  theme_bw()+
  facet_wrap(~trait)+ylim(c(0.02,0.25))+
  guides(color=guide_legend(override.aes = aes(label = "")))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), legend.position = 'top')+
  xlab('')+ylab('Mean E-P ABC score')+
  stat_compare_means(method = 'wilcox.test',label ="p.signif" ) 

dev.off()


# p4. e-p distance

ep_dist <- list()
for (Trait in setdiff(unique(final_link$trait),c("insomn2",'adhd')) ){
  print(Trait)
  ep_dist[[Trait]] <-  get_distance(Trait)
}
ep_dist <- as.data.frame(do.call(rbind,ep_dist))

pdf(file = "../figures/combined_trait/p4.ep_dist.pdf", width = 7.5, height = 5)
ep_dist$trait <- trait_info[ep_dist$trait,]
ggplot(subset(ep_dist,promoter!='Both'), aes(x=factor(brain, levels = Br_color$subRG),
                                             y=dis/1000, fill=factor(promoter, levels = c("5'","Non 5'",'Both'))))+
  geom_boxplot(outlier.size = 0.1)+
  scale_fill_manual('Promoter isoform',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF"))+
  theme_bw()+
  facet_wrap(~trait)+
  guides(color=guide_legend(override.aes = aes(label = ""))) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), legend.position = 'top')+
  xlab('')+ylab('Mean E-P distance (kb)')+ylim(c(0,100))+
  stat_compare_means(label ="p.signif", method = 'wilcox.test') 

dev.off()


# p5. pops score
pop_trait <- list()
for (Trait in setdiff(unique(final_link$trait),c("insomn2",'adhd','smoking')) ){
  print(Trait)
  pop_trait[[Trait]] <-  get_pops(Trait)
}

pops_trait <- as.data.frame(do.call(rbind,pop_trait))
pops_trait$PoPS <- as.numeric(pops_trait$PoPS)
library(ggsignif)

pdf(file = "../figures/combined_trait/p5.pops_lead.pdf", width = 8, height = 7)

pops_trait$trait <- trait_info[pops_trait$trait,]

pops_trait$brain <- factor(pops_trait$brain, levels = Br_color$subRG)

ggplot(pops_trait, aes(x=factor(promoter, levels = c("5'","Non 5'",'Both','Closest gene')),
                       y=PoPS,fill=promoter))+
  facet_grid(trait~brain)+
  geom_boxplot(outlier.size = 0.1) +
  scale_color_manual('Promoter isoform',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF",'Closest gene'='orange'))+
  scale_fill_manual('Promoter isoform',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF",'Closest gene'='orange'))+
  guides(color=guide_legend(override.aes = aes(label = "")))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), legend.position = '')+
  #stat_compare_means(label ="p.signif"  ) +
  geom_signif(test='wilcox.test',
    comparisons = list(
      c("5'", "Non 5'"),
      #c("Non 5'", "Both"),
      c("Non 5'", "Closest gene")
    ),#,#map_signif_level = function(p) sprintf("%.2g", p)
    map_signif_level=c("***"=0.001,"**"=0.01, "*"=0.05, "ns."=2))+
  xlab('')+ylab('PoPS score')+ylim(c(0,1.2))

dev.off()


  
## p6. closest gene of lead snp #######################
#gsums_sig_38_trait$index <- paste(gsums_sig_38_trait$seqnames,gsums_sig_38_trait$start)

props <- c()
for (tt in setdiff(unique(final_link$trait),c("insomn2",'adhd'))){
  lead_snp_t <- get_lead(tt)
  closest_gene_t <- unique((subset(closest_gene, index %in%lead_snp_t$index))[,'closest_ensG'])
  
    target_isoform_tt <- target_isoform[[tt]]
    prop <- apply(target_isoform_tt, 2, function(x){ 
      final_link_tt_gene <- rownames(target_isoform_tt)[x!=0]
      prop_new <- 100*(length(setdiff(final_link_tt_gene,closest_gene_t))/length(final_link_tt_gene))
      return(prop_new)
      })
    
    prop <- cbind(prop,tt,colnames(target_isoform_tt))
    props <- rbind(prop,props)
}

props <- as.data.frame(props)
props$prop <- as.numeric(props$prop)
props$V3 <- Br_color[props$V3,'subRG']

mean(props$prop)
sd(props$prop)/sqrt(6*9)

pdf(file = "../figures/combined_trait/p6.nonclosest.pdf", width = 7.5, height = 4)
props$tt <- trait_info[props$tt,]
ggplot(props, 
       aes(x=factor(V3,levels = Br_color$subRG),y=prop))+
  facet_wrap(~tt)+
  geom_bar(stat = 'identity',width = 0.6,fill='lightblue')+
  scale_fill_manual('Promoter isoform',values= c("5'"="#FEE090FF" , "Non 5'"="#4393C3FF"  ,'Both'="#92C5DEFF"))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), legend.position = 'top')+
  xlab('')+ylab('Proportion of non-closest targets (%)')# ('Number of target genes')
dev.off()



# p7.eqtl validation proportion
eqtl_prop <- c()
for (tt in setdiff(unique(final_link$trait),c("insomn2",'adhd'))){ 
  eqtl_prop <- c( eqtl_prop, get_eqtl(tt) )
}
eqtl_prop <- as.data.frame( cbind(eqtl_prop,setdiff(unique(final_link$trait),c("insomn2",'adhd'))))

pdf(file = "../figures/combined_trait/p7.eqtl_vali.pdf", width = 4, height = 4)

ggplot(eqtl_prop, 
       aes(x=V2,y=as.numeric(eqtl_prop)*100))+
  geom_bar(stat = 'identity',width = 0.6,fill='lightblue')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5), legend.position = '')+
  xlab('')+ylab('Proportion of eqtl supported gene')# ('Number of target genes')
dev.off()



