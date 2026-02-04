##show PAM50 expression difference
# BRCA_list$TCGA$exp
library(ggsci)
colpalettes<-c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
               pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
               pal_locuszoom("default")(7),pal_igv("default")(51),
               pal_uchicago("default")(9),pal_startrek("uniform")(7),
               pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
               pal_simpsons("springfield")(16),pal_gsea("default")(12))

colpalettes=colpalettes[!duplicated(colpalettes)]
library(circlize)
col_fun = colorRamp2(c(-3,0,3), c("blue", "white","red"))

dataname0="20260128"
for(database in c("TCGA","metabric","SCANB")){
  print(database)
  exp=BRCA_list[[database]]$exp
  clinical=BRCA_list[[database]]$clinical
  if(database=="TCGA"){clinical=clinical[!clinical$normalORnot%in%"normal",];exp=log(exp+1)}
  clinical=clinical[!(is.na(clinical$PAM50)|clinical$PAM50==""),]
  exp=exp[,clinical$ID]
  
  select_var <- order(apply(exp, 1, var), decreasing = TRUE)[1:1000]
  exp <- exp[select_var, ]
  
  
  ##clustering
  library(ComplexHeatmap)
  annot_df_col=clinical[,c("PAM50"),drop=F]
  column_color3=colpalettes[1:(length(as.character(unique(annot_df_col$PAM50))))]
  names(column_color3)=levels(annot_df_col$PAM50)
  
  anno_col=HeatmapAnnotation(df=annot_df_col,col=list(PAM50=column_color3))
  
  expMat=as.matrix(exp)
  expMat=t(scale(t(expMat)))
  heatmapresult=Heatmap(expMat, name = "scaledata", cluster_rows = T, cluster_columns = T, top_annotation = anno_col, use_raster=T,
                        column_split = annot_df_col$PAM50,
                        show_column_dend = F, show_row_dend = F, #row_split =geneloc$V2[match(rownames(expMat),geneloc$V1)],#row_gap = 0,
                        clustering_distance_rows='euclidean', show_row_names = F,show_column_names = F,
                        col =col_fun,
                        column_dend_height = unit(40, "mm"),#row_names_gp = gpar(fontsize=10),
                        row_km = 1)
  pdf(paste0(dataname0,database,".pdf"))
  h1.1=draw(heatmapresult)
  dev.off()
  
}


load("X:/analyze_data/FEdata/liuxiaoqin/20250815_BRCA_db/tabledataset/surv_df_list.RData")

genesurvdf=surv_df_list$gene
genesurvdf=lapply(genesurvdf,function(x){x[x$eventype%in%"OS"&x$variable%in%c("total","PAM50"),]})
# genesurvdf$Lancet2005=surv_df_list$gene$Lancet2005
genesurvdf=lapply(genesurvdf,function(x){x[,c("database","gene","variable","grouplevel","eventype","pval","mediantimehigh","mediantimelow","coef","coxP","highlow")]})
genesurvdf=genesurvdf[1:3]
writexl::write_xlsx(genesurvdf,path = "TCGA.metabric.SCANB.gene.logrank.xlsx")

genelist=list()

for(database in names(genesurvdf)){
  for(variable in c("total","PAM50","ER","PR","HER2","stage")){
    t1=genesurvdf$SCANB[genesurvdf$SCANB$variable%in%variable,]
    grouplevels=t1$grouplevel[!duplicated(t1$grouplevel)]
    tmplist1=lapply(genesurvdf,function(x){x[x$variable%in%variable&x$pval<0.05,]})
    for(grouplevel in grouplevels){
      tmplist2=lapply(tmplist1,function(x){x[x$grouplevel%in%grouplevel,]})
      tmplist3=lapply(tmplist2,function(x){list(good=x$gene[x$highlow==T],poor=x$gene[x$highlow==F])})
      genelist[[variable]][[grouplevel]]=tmplist3
    }    
  }
}


pdf("good.poor.intersect.pdf",width = 20)
b=fromList(t[grep("good",names(t))]);
upset(b,sets=colnames(b),set_size.show = T,keep.order = T,nintersects = NA)
b=fromList(t[grep("poor",names(t))]);
upset(b,sets=colnames(b),set_size.show = T,keep.order = T,nintersects = NA)
dev.off()


pdf("gene.intersect.TCGA.SCANB.pdf")
genelistintersect=list()
for(variable in names(genelist)){
  for(grouplevel in names(genelist[[variable]])){
    print(paste(variable,grouplevel))
    t=genelist[[variable]][[grouplevel]]
    genelistintersect[[variable]][[grouplevel]][["good"]]=intersect(intersect(t$TCGA$good,t$SCANB$good),t$metabric$good)
    genelistintersect[[variable]][[grouplevel]][["poor"]]=intersect(intersect(t$TCGA$poor,t$SCANB$poor),t$metabric$poor)
    # p1=ggvenn::ggvenn(lapply(t,function(x){x$good}),fill_color = c(brewer.pal(n = 8, name = "Set1")))+labs(title = paste(variable,grouplevel))
    # p2=ggvenn::ggvenn(lapply(t,function(x){x$poor}),fill_color = c(brewer.pal(n = 8, name = "Blues")))
    # print(p1+p2)
  }
}
dev.off()
genelistintersect[[variable]][[grouplevel]][["good"]]=intersect(t$TCGA$good,t$SCANB$good)
genelistintersect[[variable]][[grouplevel]][["poor"]]=intersect(t$TCGA$poor,t$SCANB$poor)
t=do.call("c",do.call("c",genelistintersect))
genelistintersect.list=t

t=lapply(t,function(x){data.frame(intersect=x)})
writexl::write_xlsx(t,path ="TCGA.metabric.SCANB.gene.intersect.xlsx")

##show some examples in three databases using showgenesurv.R




##function enrichment
t=do.call("c",do.call("c",genelistintersect))
markerlist=t

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

load("human.gene2id.RData")
enrichfunciton=list()
for(celltype in names(markerlist)){
  
  print(paste(celltype))
  genes=markerlist[[celltype]]
  
  geneID=gene2id$ENTREZID[match(genes,gene2id$SYMBOL,nomatch=0)]
  if(length(geneID)<5){next()}
  ego=enrichGO(gene = geneID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
  ekegg <-enrichKEGG(gene = geneID, organism = 'hsa', pAdjustMethod = "BH", qvalueCutoff = 0.05)
  
  if(!is.null(ego)){
    ego=ego@result[ego@result$pvalue<0.05,]
    enrichfunciton[[celltype]][["ego"]]=ego
  }
  if(!is.null(ekegg)){
    ekegg=ekegg@result[ekegg@result$pvalue<0.05,]
    ekegg$genesymbol=sapply(ekegg$geneID,function(x){a=unlist(strsplit(x,"/"));b=gene2id$SYMBOL[match(a,gene2id$ENTREZID,nomatch = 0)];return(paste0(b,collapse = "/"))})
    enrichfunciton[[celltype]][["ekegg"]]=ekegg
    
  }
}
function_BRCA_intersect3=list(fun=enrichfunciton,gene=markerlist)
save(function_BRCA_intersect3,file = "function_BRCA_intersect3.RData")
t=do.call("c",enrichfunciton)
writexl::write_xlsx(t,path = "function_BRCA_intersect3.xlsx")

load("function_BRCA_intersect3.RData")
showfun=list(
  total.total.good.ekegg=c("PPAR signaling pathway","p53 signaling pathway","Apoptosis","ECM-receptor interaction","Fatty acid degradation"),
  total.total.poor.ekegg=c("Cell cycle","HIF-1 signaling pathway","Endocrine resistance","Glycolysis / Gluconeogenesis","PI3K-Akt signaling pathway"),
  PAM50.LumA.good.ekegg=c("ECM-receptor interaction","Focal adhesion","Integrin signaling","Relaxin signaling pathway","Biosynthesis of unsaturated fatty acids"),
  PAM50.LumA.poor.ekegg=c("ABC transporters","Phosphatidylinositol signaling system","Protein processing in endoplasmic reticulum","Estrogen signaling pathway","Proteasome","cAMP signaling pathway","Central carbon metabolism in cancer","Pentose phosphate pathway"),
  PAM50.LumB.good.ekegg=c("PPAR signaling pathway","Biosynthesis of unsaturated fatty acids","Focal adhesion","One carbon pool by folate","Glycine, serine and threonine metabolism","Regulation of lipolysis in adipocytes","Th1 and Th2 cell differentiation","Fatty acid metabolism","Cellular senescence","Th17 cell differentiation","Adherens junction","Fatty acid elongation","AMPK signaling pathway","p53 signaling pathway"),
  PAM50.LumB.poor.ekegg=c("Antifolate resistance","Necroptosis","Autophagy - animal","Proteasome","Wnt signaling pathway","NOD-like receptor signaling pathway","AGE-RAGE signaling pathway in diabetic complications","Endocytosis","Retrograde endocannabinoid signaling","Mitophagy - animal","Nucleotide excision repair"),
  PAM50.Her2.good.ekegg=c("Cytokine-cytokine receptor interaction","Hematopoietic cell lineage","Th1 and Th2 cell differentiation","B cell receptor signaling pathway","Th17 cell differentiation","Antigen processing and presentation","Natural killer cell mediated cytotoxicity","T cell receptor signaling pathway","NF-kappa B signaling pathway","PD-L1 expression and PD-1 checkpoint pathway in cancer","Calcium signaling pathway","Peroxisome"),
  PAM50.Her2.poor.ekegg=c("Wnt signaling pathway","MAPK signaling pathway","Endocytosis","Autophagy - animal","ABC transporters","ErbB signaling pathway","Oxidative phosphorylation","Choline metabolism in cancer","AGE-RAGE signaling pathway in diabetic complications","Cellular senescence","PI3K-Akt signaling pathway","Autophagy - other","Platinum drug resistance"),
  PAM50.Basal.good.ekegg=c("Th1 and Th2 cell differentiation","Antigen processing and presentation","Th17 cell differentiation","Cell adhesion molecule (CAM) interaction","Cytokine-cytokine receptor interaction","Phagosome","TNF signaling pathway","Natural killer cell mediated cytotoxicity","T cell receptor signaling pathway","Efferocytosis","Chemokine signaling pathway"),
  PAM50.Basal.poor.ekegg=c("TGF-beta signaling pathway","PI3K-Akt signaling pathway","Focal adhesion","p53 signaling pathway","Rap1 signaling pathway","MicroRNAs in cancer","Signaling pathways regulating pluripotency of stem cells","Endocrine resistance","Notch signaling pathway","Integrin signaling","Cellular senescence"),
  PAM50.Normal.good.ekegg=c("FoxO signaling pathway","Butanoate metabolism","Apoptosis","Linoleic acid metabolism","Regulation of lipolysis in adipocytes","Longevity regulating pathway - multiple species","Drug metabolism - cytochrome P450"),
  PAM50.Normal.poor.ekegg=c("Cell cycle","Pentose phosphate pathway","p53 signaling pathway","Protein processing in endoplasmic reticulum","Glycolysis / Gluconeogenesis","DNA replication","Regulation of actin cytoskeleton","Proteoglycans in cancer","Central carbon metabolism in cancer")
)

pdf("function_BRCA_intersect3.show.pdf")
for(i in c("total.total","PAM50.LumA","PAM50.LumB","PAM50.Her2","PAM50.Basal","PAM50.Normal")){
  print(i)
  tlist=list()
  for(j in c("good","poor")){
    t=function_BRCA_intersect3$fun[[paste(i,j,sep = ".")]]$ekegg
    t=t[t$Description%in%showfun[[paste(i,j,"ekegg",sep = ".")]],]
    
    tlist[[j]]=t
    
  }
  df=do.call("rbind",tlist)
  df$type="good"
  df$type[grep("poor",rownames(df))]="poor"
  df$Count[df$type=="poor"]=-1*df$Count[df$type=="poor"]
  df=df[order(df$Count),]
  df$Description=factor(df$Description,levels = df$Description[!duplicated(df$Description)])
  p=ggplot(df,aes(Count,Description,fill=pvalue))+geom_bar(stat = "identity")+theme_classic()+
    scale_fill_gradient(low = "red",high = "blue")+
    labs(title = i)
  print(p)
}
dev.off()


###cibersort analysis
library(CIBERSORT)#读取包自带的LM22文件（免疫细胞特征基因文件）
sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
#样本表达矩阵文件
mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")

data(LM22)

data(mixed_expr)
results <- cibersort(sig_matrix = LM22, mixture_file = mixed_expr,perm = 1000,QN = F)

explist=lapply(BRCA_list[1:3],function(x){x$exp})
explist$metabric=exp(explist$metabric)
explist$SCANB=exp(explist$SCANB)
cibersort_res=list()
for(database in names(explist)){
  print(database)
  results=cibersort(sig_matrix = LM22, mixture_file = explist[[database]],perm = 1000,QN = F)
  cibersort_res[[database]]=results
}
save(cibersort_res,file = "cibersort_res.RData")

# genelistintersect.list score
genelistintersect.list.score=list()
library(GSVA)
for(database in c("TCGA","metabric","SCANB")){
  print(database)
  exp=BRCA_list[[database]]$exp
  meanscore=sapply(genelistintersect.list,function(x){apply(as.matrix(exp[x,]),2,mean)})
  ssGSEA_scores <- gsva(as.matrix(exp), genelistintersect.list, kcdf = "Gaussian",method = "ssgsea", verbose = TRUE)
  genelistintersect.list.score[[database]][["mean"]]=meanscore
  genelistintersect.list.score[[database]][["ssgsea"]]=ssGSEA_scores
  
}
save(genelistintersect.list.score,file = "genelistintersect.list.score.RData")


load("cibersort_res.RData")
cibersort_res=lapply(cibersort_res,as.data.frame)
cibersort_res_sig=lapply(cibersort_res,function(x){x[x$`P-value`<0.05&x$Correlation>0.25,]})
lapply(names(cibersort_res_sig),function(x){a=BRCA_list[[x]]$clinical;return(table(a$PAM50[match(rownames(cibersort_res_sig[[x]]),a$ID)]))})
# [[1]]
# LumA   LumB   Her2  Basal Normal 
# 425    138     63    142     33 
# [[2]]
# LumA        LumB        Her2       Basal claudin-low      Normal 
# 70          74          81          82         145          27 
# [[3]]
# LumA   LumB   Her2  Basal Normal 
# 2217   1417    908    691    982 
##celltype fraction comparison among subtypes in three dataset

siglist=cibersort_res_sig
siglist=lapply(names(siglist),function(x){a=BRCA_list[[x]]$clinical;cibersort_res_sig[[x]]$PAM50=a$PAM50[match(rownames(cibersort_res_sig[[x]]),a$ID)];return(cibersort_res_sig[[x]])})
names(siglist)=names(cibersort_res_sig)

pdf("cibersort_subtype.diff.pdf",width = 10,height = 10)
for(database in names(siglist)){
  print(database)
  t=siglist[[database]]

  
  df=reshape2::melt(t[,c(colnames(t)[1:22],"PAM50")])
  df$variable=factor(df$variable,levels = colnames(t)[1:22])
  df=df[df$PAM50%in%c(levels(siglist$TCGA$PAM50)),]
  df$PAM50=droplevels(df$PAM50)
  p=ggplot(df,aes(PAM50,value))+geom_violin(aes(fill=PAM50))+geom_boxplot(width=0.1,outlier.colour = NA)+
    facet_wrap(~variable,nrow = 5,scales = "free_y")+
    ggpubr::stat_compare_means(comparisons = as.list(data.frame(combn(levels(df$PAM50),2))),method = "t.test")+
    theme_classic()+scale_fill_manual(values = colpalettes)+labs(title = database)
  print(p)
}
dev.off()

##correlation between celltype fraction and geneset score in each subtype, show correlation heatmap

library(psych)
corlist=list()
for(database in names(siglist)){
  t=siglist[[database]]
  for(method in c("ssgsea","mean")){
    score=genelistintersect.list.score[[database]][[method]]
    if(method=="ssgsea"){score=t(score)}
    score=score[rownames(t),]
    
    cor.result <- corr.test(t[,1:22],score[,grep("total",colnames(score))],method = "pearson",adjust = "BH")
    corlist[[database]][[method]][['total']]=cor.result
    
    for(subtype in levels(siglist$TCGA$PAM50)){
      print(paste(database,method,subtype))
      t1=t[t$PAM50%in%subtype,]
      cor.result <- corr.test(t1[,1:22],score[rownames(t1),grep(subtype,colnames(score))],method = "pearson",adjust = "BH")
      corlist[[database]][[method]][[subtype]]=cor.result
    } 
  }  
}

toshow=list(total.macrophageM0=c("total","total.total.poor","Macrophages M0"),
            c("Her2","PAM50.Her2.good","T cells CD8"),
            c("Her2","PAM50.Her2.good","T cells CD4 memory activated"),
            c("Her2","PAM50.Her2.good","NK cells activated"),
            c("Her2","PAM50.Her2.good","Macrophages M2"),
            c("Her2","PAM50.Her2.good","Macrophages M0"),
            c("Her2","PAM50.Her2.poor","T cells CD4 memory activated"),
            c("Her2","PAM50.Her2.poor","Macrophages M2"),
            c("Her2","PAM50.Her2.poor","T cells CD8"),
            
            c("Basal","PAM50.Basal.good","Macrophages M1"),
            c("Basal","PAM50.Basal.good","T cells CD8"),
            c("Basal","PAM50.Basal.good","T cells CD4 memory activated"),
            c("Basal","PAM50.Basal.good","NK cells activated"),
            c("Basal","PAM50.Basal.good","Macrophages M0"),
            c("Basal","PAM50.Basal.good","Macrophages M2"),
            c("Basal","PAM50.Basal.poor","T cells CD4 memory activated"),
            
            c("Normal","PAM50.Normal.poor","T cells regulatory (Tregs)"),
            c("Normal","PAM50.Normal.poor","Monocytes")
            
            )
pdf("cibersort.prognosisgene.correlation.pdf",width = 15,height = 5)
for(i in 1:length(toshow)){
  plist=list()
  for(database in names(siglist)){
    t=siglist[[database]]
    score=t(genelistintersect.list.score[[database]]$ssgsea)
    score=score[rownames(t),]
    
    comb=cbind(t,score)
    if(toshow[[i]][[1]]=="total"){
      df=comb[,c(toshow[[i]][[2]],toshow[[i]][[3]])]
      colnames(df)=c("prognosisscore","cellpercentage")
      p=ggplot(df,aes(prognosisscore,cellpercentage))+geom_point()+theme_classic()+ggpubr::stat_cor()+
        stat_smooth(se=F,method="lm")+labs(title = paste(database,toshow[[i]][[2]],toshow[[i]][[3]]))
      # print(p)
      plist[[database]]=p
    }else{
      df=df=comb[comb$PAM50%in%toshow[[i]][[1]],c(toshow[[i]][[2]],toshow[[i]][[3]])]
      colnames(df)=c("prognosisscore","cellpercentage")
      p=ggplot(df,aes(prognosisscore,cellpercentage))+geom_point()+theme_classic()+ggpubr::stat_cor()+
        stat_smooth(se=F,method="lm")+labs(title = paste(database,toshow[[i]][[2]],toshow[[i]][[3]]))
      # print(p)
      plist[[database]]=p
    }
  }
  print(plist[[1]]+plist[[2]]+plist[[3]])
}
dev.off()



pdf("cibersort.prognosisgene.correlation.pdf")
library(ComplexHeatmap)
rlist=lapply(do.call("c",do.call("c",corlist)),function(x){x$r})
plist=lapply(do.call("c",do.call("c",corlist)),function(x){x$p})
for(i in names(rlist)){
  print(i)
  if(grepl("mean",i)){next()}
  expmat=rlist[[i]]
  expmat[plist[[i]]>0.05]=0
  # expmat[abs(expmat)<0.4]=0
  # expmat=expmat[apply(expmat,1,function(x){sum(is.na(x))})<2,]
  expmat=expmat[!is.na(expmat[,1]),]
  
  p=Heatmap(expmat,column_title  = i)
  print(p)
}
dev.off()


