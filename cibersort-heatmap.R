genelistintersect.list=list()

for(i in c("total.total.good","total.total.poor",paste0(rep(paste0("PAM50.",c("LumA","LumB","Her2","Basal","Normal")),each=2),c(".good",".poor")))){
  print(i)
  t=xlsx::read.xlsx("TCGA.metabric.SCANB.gene.intersect.xlsx",sheetName = i)
  genelistintersect.list[[i]]=t
}
genelistintersect.list=lapply(genelistintersect.list,function(x){x[,1]})

load("X:/analyze_data/FEdata/liuxiaoqin/20250815_BRCA_db/BRCA_list.RData")
explist=lapply(BRCA_list[1:3],function(x){x$exp})
explist$metabric=exp(explist$metabric)
explist$SCANB=exp(explist$SCANB)

genelistintersect.list.score=list()
library(GSVA)
for(database in c("TCGA","metabric","SCANB")){
  print(database)
  exp=BRCA_list[[database]]$exp
  meanscore=sapply(genelistintersect.list,function(x){apply(as.matrix(exp[x,]),2,mean)})
  ssGSEA_scores <- gsva(as.matrix(exp), genelistintersect.list, kcdf = "Gaussian",method = "ssgsea", verbose = TRUE)
  # genelistintersect.list.score[[database]][["mean"]]=meanscore
  genelistintersect.list.score[[database]][["ssgsea"]]=ssGSEA_scores
  
}
save(genelistintersect.list.score,file = "genelistintersect.list.score.RData")


load("X:/analyze_data/FEdata/liuxiaoqin/20260123_BRCA_metastasis/figures/01-pdfs/cibersort_res.RData")
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

##correlation between celltype fraction and geneset score in each subtype, show correlation heatmap

library(psych)
corlist=list()
for(database in names(siglist)){
  t=siglist[[database]]
  for(method in c("ssgsea")){#,"mean"
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



toshow=list(c("total","total.total.poor","Macrophages M0"),
            c("total","total.total.good","T cells CD4 memory resting"),
            c("total","total.total.good","T cells regulatory (Tregs)"),
            c("total","total.total.good","Macrophages M0"),
            
            c("LumA","PAM50.LumA.good","Macrophages M0"),
            c("LumA","PAM50.LumA.poor","Macrophages M0"),
            
            c("Her2","PAM50.Her2.poor","Macrophages M2"),
            c("Her2","PAM50.Her2.good","Macrophages M2"),
            c("Her2","PAM50.Her2.good","Macrophages M0"),
            c("Her2","PAM50.Her2.good","T cells CD8"),

            c("Basal","PAM50.Basal.good","Macrophages M0"),
            c("Basal","PAM50.Basal.poor","T cells CD4 memory activated"),
            c("Basal","PAM50.Basal.good","T cells CD8"),
            c("Basal","PAM50.Basal.good","Macrophages M1"),
            c("Basal","PAM50.Basal.good","Macrophages M2")
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
padjlist=lapply(do.call("c",do.call("c",corlist)),function(x){x$p.adj})
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


cordf=data.frame(matrix(0,ncol = 6))
colnames(cordf)=c("TCGA","metabric","SCANB","datarange","gsname","celltype")
# pvaluedf=data.frame(matrix(0,ncol = ))

for(i in 1:length(toshow)){
  type=toshow[[i]][1]
  corvalue=sapply(c("TCGA","metabric","SCANB"),function(x){rlist[[paste0(x,".ssgsea.",type)]][toshow[[i]][3],toshow[[i]][2]]})
  pvalue=sapply(c("TCGA","metabric","SCANB"),function(x){plist[[paste0(x,".ssgsea.",type)]][toshow[[i]][3],toshow[[i]][2]]})
  padjvalue=sapply(c("TCGA","metabric","SCANB"),function(x){padjlist[[paste0(x,".ssgsea.",type)]][toshow[[i]][3],toshow[[i]][2]]})
  df=data.frame(corvalue,pvalue,padjvalue)
  df=as.data.frame(t(df))
  df$datarange=toshow[[i]][1]
  df$gsname=toshow[[i]][2]
  df$celltype=toshow[[i]][3]
  cordf=rbind(cordf,df)
}
# rownames(cordf)=paste(cordf$gsname,cordf$celltype,sep = " vs ")
pvaluedf=as.matrix(cordf[grep("pvalue",rownames(cordf)),1:3])
padjvaluedf=as.matrix(cordf[grep("padjvalue",rownames(cordf)),1:3])
cordf=cordf[grep("corvalue",rownames(cordf)),]
rownames(cordf)=paste0(cordf$datarange,":",cordf$gsname, " vs ",cordf$celltype)

annotation_text <- matrix(
  ifelse(pvaluedf < 0.001, paste0(sprintf("%.3f", pvaluedf), "***"),
         ifelse(pvaluedf < 0.01, paste0(sprintf("%.3f", pvaluedf), "**"),
                ifelse(pvaluedf < 0.05, paste0(sprintf("%.3f", pvaluedf), "*"),
                       sprintf("%.3f", pvaluedf)))),
  nrow = nrow(pvaluedf),
  ncol = ncol(pvaluedf),
  dimnames = dimnames(pvaluedf[,1:3])
)
annotation_text <- matrix(
  ifelse(padjvaluedf < 0.001, paste0(sprintf("%.3f", padjvaluedf), "***"),
         ifelse(padjvaluedf < 0.01, paste0(sprintf("%.3f", padjvaluedf), "**"),
                ifelse(padjvaluedf < 0.05, paste0(sprintf("%.3f", padjvaluedf), "*"),
                       sprintf("%.3f", padjvaluedf)))),
  nrow = nrow(padjvaluedf),
  ncol = ncol(padjvaluedf),
  dimnames = dimnames(padjvaluedf[,1:3])
)
pdf("cor.heatmap.pdf",width = 4)
pheatmap::pheatmap(
  mat = cordf[,1:3],                # 相关性矩阵（热图颜色依据）
  display_numbers = annotation_text,  # 格子上标注p值
  number_color = "black",          # p值文字颜色（对比热图颜色更清晰）
  number_fontsize = 20,             # p值文字大小（根据矩阵大小调整）
  main = "Correlation Heatmap with P-adjust showing in cell",  # 标题
  color = colorRampPalette(c("#2166ac", "#b2182b"))(100),  # 红蓝配色（经典cor热图）
  border_color = "white",          # 格子边框（白色更整洁）
  cellwidth = 30,                  # 格子宽度（根据矩阵大小调整）
  # cellheight = 30,                 # 格子高度
  treeheight_row = 0,              # 关闭行聚类树（若需要聚类则设为默认值）
  treeheight_col = 0,              # 关闭列聚类树
  show_rownames = TRUE,            # 显示行名
  show_colnames = TRUE,            # 显示列名
  fontsize = 10,                   # 行列名字体大小
  legend_labels = c("-1", "0", "1"),  # 相关性图例标签
  legend_breaks = c(-1, 0, 1)      # 相关性图例断点
)
dev.off()
