load("X:/analyze_data/FEdata/liuxiaoqin/20250815_BRCA_db/BRCA_list.RData")

gene="TSPAN6"
database="TCGA"

showgenesurv(gene,database)

for(gene in c("CASP9","FN1")){
  for(database in c("TCGA","metabric","SCANB")){
    showgenesurv(gene,database)
  }
}

# gs=genelistintersect.list
gs=hsa_kegg_list_des[c("Cell cycle","Estrogen signaling pathway")]

for(database in c("TCGA","metabric","SCANB")){
  
  exp=BRCA_list[[database]]$exp
  scoreMat=sapply(gs,function(x){apply(as.matrix(exp[match(x,rownames(exp),nomatch = 0),]),2,median)})
  
  for(gsname in names(gs)){
    print(paste(database,gsname))
    showpathsurv(gsname=gsname,score=scoreMat[,gsname,drop=F],database=database)
    
  }
}

gs=genelistintersect.list[grep("Lum",names(genelistintersect.list))]
for(database in c("Lancet2005")){
  
  exp=BRCA_list[[database]]$exp
  scoreMat=sapply(gs,function(x){apply(as.matrix(exp[match(x,rownames(exp),nomatch = 0),]),2,mean)})
  
  for(gsname in names(gs)){
    print(paste(database,gsname))
    showpathsurv(gsname=gsname,score=scoreMat[,gsname,drop=F],database=database)
    
  }
}

for(database in c("TCGA","metabric","SCANB")){
  
  exp=BRCA_list[[database]]$exp
  scoreMat=sapply(gs,function(x){apply(as.matrix(exp[match(x,rownames(exp),nomatch = 0),]),2,mean)})
  clinical=BRCA_list[[database]]$clinical
  for(gsname in names(gs)){
    print(paste(database,gsname))
    df=data.frame(exp=scoreMat[,gsname])
    df=df[rownames(df)%in%clinical$ID,,drop=F]
    
    df=cbind(df,clinical[match(rownames(df),clinical$ID),])
    showexp(df,database=database,expvsvar,gsname)
    
  }
}

showexp=function(df,database,expvsvar,gsname){

  
  df0=df
  for(i in unique(do.call("c",expvsvar))){
    if(!i%in%colnames(df)){next()}
    print(i)
    df=df0
    if(database=="TCGA"&(i!="normalORnot")){
      df=df0[df0$normalORnot%in%c("primary","metastasis"),]
    }
    t=df[,c("exp",i)]
    colnames(t)[2]="type"
    t=t[!is.na(t$type),]
    t$type=droplevels(t$type)
    complist=list(c(levels(t$type)))
    if(length(levels(t$type))>2){
      complist=as.list(data.frame(combn(levels(t$type),2)))
    }
    
    p=ggplot(t,aes(type,exp,col=type))+geom_violin()+geom_boxplot(width=0.1,outlier.colour = NA)+theme_classic()+
      ggpubr::stat_compare_means(comparisons = complist)+
      labs(title=paste(i,gsname))+ylab(ifelse(length(gs)==1,"gene expression",paste0("geneset score by ",method)))+
      theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0),axis.title.x = element_blank(),
            # plot.background = element_rect(color = "black",fill = NA,linewidth = 1),
            legend.position = "none",plot.title = element_text(hjust = 0.5)
      )+scale_colour_brewer(palette = "Set1")
    # plist[[i]]=p
    # svg(paste0(output_path,"/",database,"-",i,"-",method,"-","geneset.score.svg"))#, width = 8, height = 6)  # width和height单位为英寸
    pdf(paste(database,gsname,i,"pdf",sep = "."))
    print(p)
    # 关闭SVG设备（必须执行，否则文件无法正常生成）
    dev.off()
    # svg_files=c(svg_files,paste0(timestamp,"/",database,"-",i,"-",method,"-","geneset.score.svg"))
  }
  

}


unicox=function(surv){
  a=survdiff(Surv(time = time,event = event)~score,data=surv)
  p.val=signif(1-pchisq(a$chisq,length(a$n)-1),3)
  
  fit=survfit(Surv(time = time,event = event)~score,data=surv)
  p=ggsurvplot(
    fit,                     # survfit object with calculated statistics.
    data = surv,             # data used to fit survival curves.
    risk.table = T,       # show risk table.
    pval = TRUE,             # show p-value of log-rank test.
    xlab = "Time in years",   # customize X axis label.
    risk.table.y.text = FALSE, # show bars instead of names in text annotations
    palette = "Set1",
    pval.method = T,surv.median.line = "hv"
  )
  # pvalue=as.numeric(gsub("p = ","",surv_pvalue(fit)$pval.txt))
  re=list(plot=p,pvalue=p.val)
  return(re)
}
expvsvar=list(TCGA=c("PAM50","ER","PR","HER2","stage","normalORnot","tissueOrigin"),
              metabric=c("PAM50","ER","PR","HER2","stage"),
              SCANB=c("PAM50","ER","PR","HER2","LN","stage","tissueOrigin"),
              Lancet2005=c("ER"),
              CCLE=c("metORnot","subtype")
)

expvsvar=list(TCGA=c("PAM50"),
              metabric=c("PAM50"),
              SCANB=c("PAM50"),
              Lancet2005=c("ER"),
              CCLE=c("metORnot","subtype")
)
showpathsurv(gsname=gsname,score=scoreMat[,gsname,drop=F],database=database)

showpathsurv=function(gsname,score,database){
  
  expdata=BRCA_list[[database]]$exp
  clinicaldata=BRCA_list[[database]]$clinical
  type="DRFS"
  # df=data.frame(exp=expdata[gsname,])
  df=data.frame(exp=score[,1])
  df=cbind(df,clinicaldata)
  
  
  library(survminer)
  library(survival)
  
  pdf(paste(database,gsname,"OS.pdf",sep = "."))
  clinical=df
  if(database=="TCGA"){
    clinical=df[!df$normalORnot%in%c("normal"),]
  }
  clinical$score=clinical$exp
  
  plist=list()
  # for(type in c(unique(do.call("c",survivaltype)))){}
  # if(!type%in%colnames(clinical)){next()}
  clinical$time=clinical[,paste0(type,".time")]
  clinical$event=clinical[,type]
  res.cut=surv_cutpoint(clinical,time = "time",event="event",variables = "score")
  res.cat=surv_categorize(res.cut)
  
  p=unicox(surv=res.cat)
  p1=p$plot+labs(title = paste(type,database,"total"))
  print(p1)
  
  tlist=list()
  for(i in unique(do.call("c",expvsvar))){
    if(!i%in%colnames(df)){next()}
    print(i)
    df=clinical
    if(database=="TCGA"&(i!="normalORnot")){
      df=df[df$normalORnot%in%c("primary","metastasis"),]
    }
    t0=df
    t0[,i]=droplevels(t0[,i])
    for(j in levels(t0[,i])){
      t=t0[t0[,i]%in%j,]
      
      t=t[!is.na(t$time),]
      if(nrow(t)<10){next()}
      res.cut=surv_cutpoint(t,time = "time",event="event",variables = "score")
      res.cat=surv_categorize(res.cut)
      p=unicox(surv=res.cat)
      p1=p$plot+labs(title = paste(type,database,i,j))
      print(p1)
    }
  }
  dev.off()
}


showgenesurv=function(gene,database){
  
  expdata=BRCA_list[[database]]$exp
  clinicaldata=BRCA_list[[database]]$clinical
  type="OS"
  df=data.frame(exp=expdata[gene,])
  df=cbind(df,clinicaldata)
  
  
  library(survminer)
  library(survival)
  
  pdf(paste(database,gene,"OS.pdf",sep = "."))
  clinical=df
  if(database=="TCGA"){
    clinical=df[!df$normalORnot%in%c("normal"),]
  }
  clinical$score=clinical$exp
  
  plist=list()
  # for(type in c(unique(do.call("c",survivaltype)))){}
  # if(!type%in%colnames(clinical)){next()}
  clinical$time=clinical[,paste0(type,".time")]
  clinical$event=clinical[,type]
  res.cut=surv_cutpoint(clinical,time = "time",event="event",variables = "score")
  res.cat=surv_categorize(res.cut)
  
  p=unicox(surv=res.cat)
  p1=p$plot+labs(title = paste(type,database,"total"))
  print(p1)
  
  tlist=list()
  for(i in unique(do.call("c",expvsvar))){
    if(!i%in%colnames(df)){next()}
    print(i)
    df=clinical
    if(database=="TCGA"&(i!="normalORnot")){
      df=df[df$normalORnot%in%c("primary","metastasis"),]
    }
    t0=df
    t0[,i]=droplevels(t0[,i])
    for(j in levels(t0[,i])){
      t=t0[t0[,i]%in%j,]
      
      t=t[!is.na(t$time),]
      if(nrow(t)<10){next()}
      res.cut=surv_cutpoint(t,time = "time",event="event",variables = "score")
      res.cat=surv_categorize(res.cut)
      p=unicox(surv=res.cat)
      p1=p$plot+labs(title = paste(type,database,i,j))
      print(p1)
    }
  }
  dev.off()
}
