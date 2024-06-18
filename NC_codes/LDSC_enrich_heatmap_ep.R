library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# for enhancers regulating genes and isoforms, respectively
# for 5' and non 5' promoter isoforms
# LDSC enrichment analysis for enhancers with ABC activity
mytraits=list(`Neuropsychiatric and behavioral`=c("sz3",'asd', "bip2", #"bdAndSczShared", 
                                                  'mdd_without_23andMe',"adhd", "drinking", "smoking",
                                                  "insomn2", "neu2"),
              Neurological=c('alz2','als','ms','pd_without_23andMe','stroke'),
              `Non-brain related`=c('bmi','cad','cd','height','ibd','ra','dm2','uc'))

mytraits_df <- data.frame()
for (term in c('Neuropsychiatric and behavioral','Neurological','Non-brain related')){
  disease <- mytraits[[term]]
  df <- cbind(rep(term, length(disease)),disease)
  mytraits_df <- rbind( mytraits_df,df)
}

rownames(mytraits_df) <- mytraits_df$disease

#mytraits_df$disease_name <- c('SCZ','ASD','BIP','BD','MDD','ADHD','Drinking','Smoking','Insomnia','Neu','AD','ALS','MS','PD','Stroke',toupper(c('bmi','cad','cd')),'Height',toupper(c('ibd','ra')),'DM','UC')

get_qvalue_coef <- function(EP,trait,isoform){
  
  
  
  ldsc_output <- ifelse(EP=='Enhancer', paste0('/sc/arion/projects/roussp01a/liting/ldsc_1k_',isoform,'/meta-files/aggregatedPartInfo.tsv') ,
                        paste0('/sc/arion/projects/roussp01a/liting/ldsc_1k_Promoter_',isoform,'/meta-files/aggregatedPartInfo.tsv'))
  
  
  aggregatedPartInfo <- read.delim(ldsc_output) 
  
  
  aggregatedPartInfo <- subset(aggregatedPartInfo, sumstatID %in% mytraits_df$disease)
  aggregatedPartInfo$disease_type <- mytraits_df[aggregatedPartInfo$sumstatID,'V1']
  aggregatedPartInfo$sumstatName=gsub('[0-9]$','',gsub(' \\(without 23andMe\\)','',gsub(' \\(all\\)','',aggregatedPartInfo$sumstatName)))
  
  
  ## coef tau
  Coef_v <- dcast( aggregatedPartInfo[,c('sumstatName','annoID','normCoefficient')],
                   annoID~sumstatName,value.var='normCoefficient')
  rownames(Coef_v) <- Coef_v$annoID
  Coef_v <- Coef_v[,-1]
  Coef_v <- Coef_v[,unique(subset(aggregatedPartInfo,disease_type==trait)[,'sumstatName'])]
  
  
  ## q value
  #aggregatedPartInfo$q <- p.adjust(10^-aggregatedPartInfo$minus_log10_p_regression, method="BH")
  aggregatedPartInfo$q <- (10^-aggregatedPartInfo$minus_log10_p_regression)*9*23
  
  p_v <- dcast( aggregatedPartInfo[,c('sumstatName','annoID','q')],
                annoID~sumstatName,value.var='q')
  
  rownames(p_v) <- p_v$annoID
  p_v <- p_v[,-1]
  #p_v <- ifelse(p_v < 0.001,'***',ifelse(p_v < 0.01,'**', ifelse(p_v < 0.05,'*','')))
  p_v <- ifelse(p_v < 0.05,'+','')
  p_v <- p_v[,unique(subset(aggregatedPartInfo,disease_type==trait)[,'sumstatName'])]
  
  
  return(list(Coef_v=Coef_v,p_v=p_v))
}

Br_color <- read.table('/hpc/users/songl05/PF_25BR/data/color_br.txt',sep='\t',comment.char = '', header = T, row.names = 1)

isof <- c('5'="ABC 5'","n5"="ABC non-5'","nonreg"="Non ABC")
p_heatmap_byIsoform_row <- function(EP,trait,Isoforms ){
  
  p_v <- c()
  Coef_v <- c()
  
  for (isoform in Isoforms){
    
    p_v1 <- get_qvalue_coef(EP=EP,trait=trait, isoform=isoform)[['p_v']][rownames(Br_color),]
    Coef_v1 <- get_qvalue_coef(EP=EP,trait=trait, isoform=isoform)[['Coef_v']][rownames(Br_color),]

    p_v <- rbind(p_v, p_v1)
    Coef_v <- rbind(Coef_v,as.matrix(Coef_v1))
    
  }
  
  p_v <- as.data.frame(p_v)
  Coef_v <- as.data.frame(Coef_v)
  
  
  #if(min(Coef_v) < 0) { col_fun = colorRamp2(unique(c(seq(min(Coef_v),0,length=4),seq(0,max(Coef_v), length=4) )),brewer.pal(7, "RdBu")) }
  #if(min(Coef_v) > 0) { col_fun = colorRamp2(sort(c(seq(min(Coef_v),max(Coef_v),length=6))), (brewer.pal(6, "Blues"))) }
  { col_fun = colorRamp2(sort(c(seq(0,max(Coef_v),length=6))), (brewer.pal(6, "Blues"))) }
  #{ col_fun = colorRamp2(sort(c(seq(0,max(Coef_v),length=6))), (brewer.pal(6, "RdBu"))) }
  
  #ht_opt$TITLE_PADDING = unit(c(8, 8), "points")
  Coef_v[Coef_v < 0] <- NA
  print(rownames(p_v1))
  Heatmap(Coef_v,col=col_fun,
          
          cluster_rows = F,cluster_columns = F, 
          show_heatmap_legend=T,
          heatmap_legend_param = list(title = 'Coef tau'),
          row_split =rep(factor( Br_color[rownames(p_v1),'subRG'],levels = c(Br_color$subRG)), length(Isoforms)) ,
          cell_fun = function(j, i, x, y, width, height, fill) {
            #if(Coef_v[i, j] > -log10(0.05))
            grid.text( p_v[i, j], x, y, gp = gpar(fontsize = 7))},
          show_column_names = T, column_title = EP,
          show_row_names = F,row_title_rot = 0,
          column_title_gp = gpar(fill = "#E8F0F4", col = "black", border = "#E8F0F4"),
          left_annotation = rowAnnotation(show_legend=T,
                                          Region= rep(Br_color[rownames(p_v1),'Region'],length(Isoforms)),
                                          Group = factor(rep(isof[Isoforms],each=9),levels = c("ABC 5'","ABC non-5'","Non ABC")),
                                          col = list(Region = c("Dien" = "#FDBF6F", 
                                                                "BasGan" = "#33A02C", 
                                                                "MidBr" = "#FF7F00",
                                                                'Limbic'='#FB9A99',
                                                                "NEX"='#E31A1C'),
                                                     Group=c("ABC 5'"="#FEE090",
                                                             "ABC non-5'"="#4575B4",
                                                             "Non ABC"="#FB9A99")))
          
  )

  }



# Neuropsychiatric and behavioral
pdf(file = "/hpc/users/songl05/PF_25BR/figures/f4_ldsc_heatmap_NPsy_E.pdf", width = 6, height = 6)
p_heatmap_byIsoform_row(EP= 'Enhancer',trait= 'Neuropsychiatric and behavioral',Isoforms=c('5','n5','nonreg') )# 
dev.off()

pdf(file = "/hpc/users/songl05/PF_25BR/figures/f4_ldsc_heatmap_NPsy_P.pdf", width = 6, height = 6)
p_heatmap_byIsoform_row(EP= 'Promoter',trait= 'Neuropsychiatric and behavioral',Isoforms=c('5','n5','nonreg') )
dev.off()

# Neurological
pdf(file = "/hpc/users/songl05/PF_25BR/figures/ldsc_heatmap_Neu_E.pdf", width = 6, height = 6)
p_heatmap_byIsoform_row(EP= 'Enhancer',trait= 'Neurological',Isoforms=c('5','n5','nonreg') )# 
dev.off()

pdf(file = "/hpc/users/songl05/PF_25BR/figures/ldsc_heatmap_Neu_P.pdf", width = 6, height = 6)
p_heatmap_byIsoform_row(EP= 'Promoter',trait= 'Neurological',Isoforms=c('5','n5','nonreg') )
dev.off()

# Non-brain related
pdf(file = "/hpc/users/songl05/PF_25BR/figures/ldsc_heatmap_NonB_E.pdf", width = 6, height = 6)
p_heatmap_byIsoform_row(EP= 'Enhancer',trait= 'Non-brain related',Isoforms=c('5','n5','nonreg') )# 
dev.off()

pdf(file = "/hpc/users/songl05/PF_25BR/figures/ldsc_heatmap_NonB_P.pdf", width = 6, height = 6)
p_heatmap_byIsoform_row(EP= 'Promoter',trait= 'Non-brain related',Isoforms=c('5','n5','nonreg') )
dev.off()
