library("GenomicRanges")
library(BiocManager)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)

setwd('/hpc/users/songl05/PF_25BR/codes')

# Abc result

load('/sc/arion/projects/CommonMind/roussp01a/INGELHEIM/atacseq/step3_neuron_width500/files/abc_piso/fine/ABC_summary.Rdata')

Br_color <- read.table('../data/color_br.txt',sep='\t',comment.char = '', header = T, row.names = 1)
Br_color$n_abc <- ''

ABC <- c()
for( br in names(GRs)){
  raw_score <- as.data.frame(GRs[[br]][,c('ID2','ABC.Score','CellType')])[,-1:-5]
  raw_score <- raw_score[order(-raw_score$ABC.Score),]
  raw_score <- raw_score[!duplicated(raw_score$ID2),]
  Br_color[br,'n_abc'] <- nrow(raw_score)
  ABC <- rbind(ABC,raw_score)
}

N_EPI <- aggregate(data=ABC,ID2~CellType, length)
max(N_EPI$ID)
min(N_EPI$ID)



ABC_gene <- c()
for( br in names(GRs)){
  raw_score <- as.data.frame(GRs[[br]][,c('ID','ABC.Score','CellType')])[,-1:-5]
  raw_score <- raw_score[order(-raw_score$ABC.Score),]
  raw_score <- raw_score[!duplicated(raw_score$ID),]
  Br_color[br,'n_abc'] <- nrow(raw_score)
  ABC_gene <- rbind(ABC_gene,raw_score)
}

N_EG <- aggregate(data=ABC_gene,ID~CellType, length)
max(N_EG$ID2)
min(N_EG$ID2)

ABC_m <- dcast(ABC,ID2~ CellType,value.var = "ABC.Score")
rownames(ABC_m) <- ABC_m$ID2
ABC_m <- ABC_m[,-1]
#ABC_m[is.na(ABC_m)] <- 0

# 
pearson_sim <- function(ix){
  A = (X[,ix[1]])
  B = (X[,ix[2]])
  return( cor(A,B) )
}  

spearman_sim <- function(ix){
  A = (X[,ix[1]])
  B = (X[,ix[2]])
  return( cor(A,B,use ='complete',method='spearman' ) )
}  

# pairwise intersection
intersect_n <- function(ix){
  A = all_id[X[,ix[1]]>0]
  B = all_id[X[,ix[2]]>0]
  return(length(intersect(A,B)))
}

all_id <- rownames(ABC_m)
X <- as.matrix(ABC_m)
n <- ncol(X) 
cmb <- expand.grid(i=1:n, j=1:n) 

# sim matrix
sim_matrix <- matrix(apply(cmb,1,spearman_sim),n,n)
rownames(sim_matrix) <- colnames(ABC_m)
colnames(sim_matrix) <- colnames(ABC_m)

# intersect_matrix
intersect_matrix <- matrix(apply(cmb,1,intersect_n),n,n)
rownames(intersect_matrix) <- colnames(ABC_m)
colnames(intersect_matrix) <- colnames(ABC_m)


# 1. heatmap for similirity
#col_fun = colorRamp2(c(0.15,0.35,0.55,0.75,0.95), rev(brewer.pal(5, "RdBu")))

# 1. heatmap for similirity

col_fun = colorRamp2(seq(max(sim_matrix),min(sim_matrix),length=5), (brewer.pal(5, "RdBu")))

pdf(file = "../figures/f4_abc_heatmap_iso.pdf", width = 9, height = 7)

Heatmap(sim_matrix,col=col_fun,show_column_names = TRUE,
        column_labels = Br_color[colnames(sim_matrix),'subRG'],
        row_labels = Br_color[rownames(sim_matrix),'subRG'],
        cluster_rows = TRUE,cluster_columns = T,show_row_names = TRUE,column_names_side = c("bottom"),
        row_names_side = c("left"),
        show_row_dend=TRUE,show_column_dend=F,
        rect_gp = gpar(type = "none"),
        row_title = NULL,column_title = NULL,
        width = ncol(sim_matrix)*unit(15, "mm"), 
        height = nrow(sim_matrix)*unit(15, "mm"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        } ,
        left_annotation = rowAnnotation(show_legend=F,
                                        Region= Br_color[colnames(ABC_m),'Region'],
                                        BroadRegion = Br_color[colnames(ABC_m),'Group'],
                                        col = list(Region = c("Dien" = "#FDBF6F", 
                                                              "BasGan" = "#33A02C", 
                                                              "MidBr" = "#FF7F00",
                                                              'Limbic'='#FB9A99',
                                                              "NEC"='#E31A1C'),
                                                   BroadRegion=c("MidDien"="#FDBF6F",
                                                           "BasGan"="#33A02C",
                                                           "ForeBr"="#FB9A99"))
        ),show_heatmap_legend=F)


draw(x = unit(0.6, "npc"), y = unit(0.8, "npc"),just="left",Legend(col_fun = col_fun, title = "Correlation"))
draw(x = unit(0.75, "npc"), y = unit(0.8, "npc"), just="left",Legend(
  labels = unique(Br_color$Region[order(Br_color$subRG)]), title = "Brain Region", 
  legend_gp = gpar(fill = unique(Br_color$subColor[order(Br_color$subRG)]))))
draw(x = unit(0.75, "npc"), y = unit(0.6, "npc"), just="left",Legend(
  labels = unique(Br_color$Group[order(Br_color$subRG)]), title = "Broad Region", 
  legend_gp = gpar(fill = unique(Br_color$GpColor[order(Br_color$subRG)]))))

dev.off()



# 2 circle plot for pairwise intersection

intersect_matrix_circ <- intersect_matrix
for (i in 1:(nrow(intersect_matrix_circ)-1)){
  for (j in (i+1):nrow(intersect_matrix_circ)){
    print(c(i,j))
    intersect_matrix_circ[i,j] <-0
  }
}

ft_matrix_circ <- melt(intersect_matrix_circ)
ft_matrix_circ <- subset(ft_matrix_circ,value!=0)


col_fun = colorRamp2(seq(min(ft_matrix_circ$value), max(subset(ft_matrix_circ,Var1!=Var2)$value),6), rev(brewer.pal(6, "RdBu")))
col <- ifelse(ft_matrix_circ$Var1!=ft_matrix_circ$Var2,col_fun(ft_matrix_circ$value), "#FFFFFF00")

col_fun = colorRamp2(seq(min(ft_matrix_circ$value), max(subset(ft_matrix_circ,Var1!=Var2)$value),length=6), rev(brewer.pal(6, "RdBu")))
col <- ifelse(ft_matrix_circ$Var1!=ft_matrix_circ$Var2,col_fun(ft_matrix_circ$value), "#FFFFFF00")


# ft_matrix_circ$Var1 <- paste0(Br_color[ft_matrix_circ$Var1,"subRG"]," (",Br_color[ft_matrix_circ$Var1,"n_abc"],')')
# ft_matrix_circ$Var2 <- paste0(Br_color[ft_matrix_circ$Var2,"subRG"]," (",Br_color[ft_matrix_circ$Var2,"n_abc"],")")

Br_color$subRG <- factor(Br_color$subRG,levels = c("NEC","Limbic","ARC","HAB", "MDT","DRN","RMTG","VTA","BasGan"))

ft_matrix_circ$Var1 <- Br_color[as.character(ft_matrix_circ$Var1),"subRG"]
ft_matrix_circ$Var2 <- Br_color[as.character(ft_matrix_circ$Var2),"subRG"]




pdf(file = "../figures/abc_stat/f4_abc_circle_iso.pdf", width = 9, height = 6)

chordDiagram(ft_matrix_circ, self.link = 1, 
             order = sort(Br_color$subRG), 
             grid.col = Br_color$subColor[order(Br_color$subRG)],
             col=col, annotationTrack = c("grid"),
             annotationTrackHeight = mm_h(5) )

for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, 
              facing = "bending.inside", niceFacing = TRUE, col = "white")
}

draw(x = unit(0.82, "npc"), y = unit(0.4, "npc"), just="left",Legend(
  labels = unique(Br_color$Region[order(Br_color$subRG)]), title = "Brain Region", 
  legend_gp = gpar(fill = unique(Br_color$subColor[order(Br_color$subRG)]))))

draw(x = unit(0.82, "npc"), y = unit(0.6, "npc"),just="left",Legend(col_fun = col_fun, title = "Number of links"))

dev.off()




# proportion of unique links

total_link <- colSums(!is.na(X))

unique_n <- c()
for (i in 1:ncol(X)){
  
  unique_link <- length(which(!is.na(X[,i]) & rowSums(is.na(X[,-i]))==(ncol(X)-1)))
  unique_n <- c(unique_n,unique_link)
}

unique_prop <- round((unique_n/total_link)*100,2)
df_unique_p <- as.data.frame(cbind(colnames(X),unique_prop))
df_unique_p$Region <- as.character(Br_color[rownames(df_unique_p),'subRG' ] )
df_unique_p$unique_prop <- as.numeric(df_unique_p$unique_prop)

pdf(file = "../figures/abc_stat/f4_uniq_abc_proportion_iso.pdf", width = 5, height = 1.5)

ggplot(df_unique_p, aes(x=reorder(Region,-unique_prop),y=unique_prop))+
  geom_bar(stat = 'identity',fill="#1F78B4", width = 0.6)+theme_bw() +
  geom_text(aes(label = signif(unique_prop,3)), vjust = 1.5, colour = "white")+
  xlab('') +ylab('Unique links (%)')+theme_classic()+ 
  theme(legend.position = '',axis.line = element_line(colour = "black",linewidth = 0.2))
#axis.text.x = element_text(angle = 0, hjust = 1,vjust = 1))

dev.off()


# Number of links
N_EG$Region <- as.character(Br_color[N_EG$CellType,'subRG' ] )

sp1 <- ggplot(N_EG, aes(x=reorder(Region,-ID),y=ID))+
  geom_bar(stat = 'identity',fill="#1F78B4", width = 0.6)+theme_bw() +
  geom_text(aes(label = ID), vjust = 1.5, size=3)+
  xlab('') +ylab('Number of E-G links')+
  theme_bw()+ theme(legend.position = '',axis.line = element_line(colour = "black",linewidth = 0.2),
                                           axis.text.x =element_text(angle = 45, hjust = 1,vjust = 1)
  )

N_EPI$Region <- as.character(Br_color[N_EPI$CellType,'subRG' ] )

sp2 <- ggplot(N_EPI, aes(x=reorder(Region,-ID2),y=ID2))+
  geom_bar(stat = 'identity',fill="#1F78B4", width = 0.6)+theme_bw() +
  geom_text(aes(label = ID2), vjust = 1.5, size=3)+
  xlab('') +ylab('Number of E-P isoform links')+
  theme_bw()+ theme(legend.position = '',axis.line = element_line(colour = "black",linewidth = 0.2),
                    axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1)
  )



# the dirstribution of degree

## distance
ABC_distance <- c()
for( br in names(GRs)){
  distance <- as.data.frame(GRs[[br]][,c('ID2','ABC.Score','D2T','CellType')])[,-1:-5]
  distance <- distance[order(-distance$ABC.Score),]
  distance <- distance[!duplicated(distance$ID2),]
  ABC_distance <- rbind(ABC_distance,distance)
}
ABC_distance$Brain_region <- factor(Br_color[ABC_distance$CellType,'subRG'], 
                                    levels=Br_color$subRG)


p1 <- ggplot(ABC_distance, aes(D2T/1000, color=Brain_region))+
  geom_density( alpha=0.01)+xlim(c(-150,150)) + 
  xlab('Distance to promoter isoforms (Kb)')+theme_bw()+
  scale_color_manual('brain region',values = Br_color$color)+
  theme( legend.key.size = unit(0.2, 'cm'))
  #scale_color_brewer(brewer.pal(9, "Paired"))+
  #scale_fill_brewer(brewer.pal(9, "Paired"))
  

## 
ABC_link <- c()
for( br in names(GRs)){
  link <- as.data.frame(GRs[[br]][,c('ID2','eID','TargetGene','TargetTranscript','CellType','ABC.Score')])[,-1:-5]
  link <- link[order(-link$ABC.Score),]
  link <- link[!duplicated(link$ID2),]
  ABC_link <- rbind(ABC_link,link)
}
ABC_link$CellType <- factor(Br_color[ABC_link$CellType,'subRG'], levels=c("NEC","Limbic","ARC","HAB", "MDT","DRN","RMTG","VTA","BasGan"))

max(aggregate(data=ABC_link,ID2~CellType,length)$ID2)
min(aggregate(data=ABC_link,ID2~CellType,length)$ID2)

mean(aggregate(data=ABC_link,ID2~CellType,length)$ID2)
sd(aggregate(data=ABC_link,ID2~CellType,length)$ID2)/3

# how many promoter isoforms for each gene
gene_transcript <- aggregate(TargetTranscript ~ TargetGene + CellType, data = unique(ABC_link[,c('TargetTranscript','TargetGene','CellType')]), length) # Equivalent
gene_transcript <- within(gene_transcript, TargetTranscript[TargetTranscript>5] <- '5+')
gene_transcript$TargetTranscript <- as.character(gene_transcript$TargetTranscript)
p2 <-ggplot(gene_transcript, aes(  x=CellType, fill=TargetTranscript))+
  geom_bar(stat = 'count',position = "fill", width = 0.7)+
  scale_fill_manual(values = brewer.pal(8, "Paired")[c(1:4,7:8)])+theme_bw()+ xlab('') +ylab('Proportion')+theme(legend.position = '', axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))
  #theme(legend.title =  element_blank(), axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))

xx <- aggregate(.~TargetTranscript  + CellType , length ,data=gene_transcript)
xxx <- aggregate(TargetTranscript ~ CellType , length ,data=gene_transcript)
xxxx <- merge(xx,xxx, by="CellType") %>%mutate(p=TargetGene/TargetTranscript.y)
xxxxx <- subset(xxxx, TargetTranscript.x==1)
xxxxx$p <- 1-xxxxx$p

(mean(xxxxx$p))*100
(sd(xxxxx$p)/3)*100

# How many promoter isoforms does each enhancer regulate?
enh_transcript <- aggregate(TargetTranscript ~ eID + CellType, data = unique(ABC_link[,c('TargetTranscript','eID','CellType')]), length) # Equivalent
enh_transcript <- within(enh_transcript, TargetTranscript[TargetTranscript>5] <- '5+')
enh_transcript$TargetTranscript <- as.character(enh_transcript$TargetTranscript)
p3 <-ggplot(enh_transcript, aes(  x=CellType, fill=TargetTranscript))+geom_bar(stat = 'count',position = "fill", width = 0.7)+
  scale_fill_manual(values = brewer.pal(8, "Paired")[c(1:4,7:8)])+theme_bw()+ xlab('') +ylab('Proportion')+theme(legend.position = '', axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))

xx <- aggregate(.~TargetTranscript  + CellType , length ,data=enh_transcript)
xxx <- aggregate(TargetTranscript ~ CellType , length ,data=enh_transcript)
xxxx <- merge(xx,xxx, by="CellType") %>%mutate(p=eID/TargetTranscript.y)
xxxxx <- subset(xxxx, TargetTranscript.x==1)
xxxxx$p <- 1-xxxxx$p

(mean(xxxxx$p))*100
(sd(xxxxx$p)/3)*100


# How many enhancers for promoter isoforms ?
transcript_enh <- aggregate( eID ~  TargetTranscript + CellType, data = unique(ABC_link[,c('TargetTranscript','eID','CellType')]), length) # Equivalent
transcript_enh <- within(transcript_enh, eID[eID>5] <- '5+')
transcript_enh$eID <- as.character(transcript_enh$eID)
p4 <- ggplot(transcript_enh, aes(  x=CellType, fill=eID))+geom_bar(stat = 'count',position = "fill", width = 0.7)+
  scale_fill_manual(values = brewer.pal(8, "Paired")[c(1:4,7:8)])+theme_bw()+ xlab('') +ylab('Proportion')+
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1))


xx <- aggregate(.~eID  + CellType , length ,data=transcript_enh)
xxx <- aggregate(eID ~ CellType , length ,data=transcript_enh)
xxxx <- merge(xx,xxx, by="CellType") %>%mutate(p=TargetTranscript/eID.y)
xxxxx <- subset(xxxx,  eID.x==1)
xxxxx$p <- 1-xxxxx$p

(mean(xxxxx$p))*100
(sd(xxxxx$p)/3)*100

# the 

library(cowplot)
pdf('../figures/supp/sp1_abc_epi_stat.pdf', width = 12, height = 6)
cowplot::plot_grid(sp2,sp1,p1,
                   p2,p3,p4, align = 'h',axis='bt',
                   nrow=2,rel_widths = c(1,1, 1.25),labels = 'auto')
dev.off()

# pdf('./figures/abc_distance.pdf', width = 9, height = 4)
# cowplot::plot_grid(p1,p2,align = 'h',axis='bt')
# dev.off()
# 
# pdf('./figures/abc_degree.pdf', width = 9, height = 4)
# cowplot::plot_grid( p2,p3,p4, nrow=1,rel_widths = c(1,1, 1.25),labels = 'auto')
# dev.off()



