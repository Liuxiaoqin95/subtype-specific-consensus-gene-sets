gs=genelistintersect.list
HRlist=list()
HRlist.continuous=list()
for(database in c("TCGA","metabric","SCANB")){
  
  exp=BRCA_list[[database]]$exp
  scoreMat=sapply(gs,function(x){apply(as.matrix(exp[match(x,rownames(exp),nomatch = 0),]),2,median)})
  
  for(gsname in names(gs)){
    print(paste(database,gsname))
    # t=showpathsurv(gsname=gsname,score=scoreMat[,gsname,drop=F],database=database)
    # t=showpathsurv.medianthreshold(gsname=gsname,score=scoreMat[,gsname,drop=F],database=database)
    t=showpathsurv.correctedforclinical(gsname=gsname,score=scoreMat[,gsname,drop=F],database=database)
    
    t$database=database
    t$gsname=gsname
    HRlist[[database]][[gsname]]=t
    
    # t=showpathsurv.continuous(gsname=gsname,score=scoreMat[,gsname,drop=F],database=database)
    # t$database=database
    # t$gsname=gsname
    # HRlist.continuous[[database]][[gsname]]=t
  }
}
gsname=gsname
score=scoreMat[,gsname,drop=F]
database=database

for(database in c("Lancet2005")){
  
  exp=BRCA_list[[database]]$exp
  scoreMat=sapply(gs,function(x){apply(as.matrix(exp[match(x,rownames(exp),nomatch = 0),]),2,median)})
  
  for(gsname in names(gs)[3:6]){
    print(paste(database,gsname))
    t=showpathsurv(gsname=gsname,score=scoreMat[,gsname,drop=F],database=database)
    # t=showpathsurv.medianthreshold(gsname=gsname,score=scoreMat[,gsname,drop=F],database=database)
    # t=showpathsurv.correctedforclinical(gsname=gsname,score=scoreMat[,gsname,drop=F],database=database)
    
    t$database=database
    t$gsname=gsname
    HRlist[[database]][[gsname]]=t
    
    # t=showpathsurv.continuous(gsname=gsname,score=scoreMat[,gsname,drop=F],database=database)
    # t$database=database
    # t$gsname=gsname
    # HRlist.continuous[[database]][[gsname]]=t
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
    palette = c("#377EB8", "#E41A1C"),#"Set1",
    pval.method = T,surv.median.line = "hv"
  )
  # pvalue=as.numeric(gsub("p = ","",surv_pvalue(fit)$pval.txt))
  re=list(plot=p,pvalue=p.val)
  return(re)
}

showpathsurv=function(gsname,score,database){
  
  expdata=BRCA_list[[database]]$exp
  clinicaldata=BRCA_list[[database]]$clinical
  # type="OS"
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
  
  # plist=list()
  # for(type in c(unique(do.call("c",survivaltype)))){}
  # if(!type%in%colnames(clinical)){next()}
  clinical$time=clinical[,paste0(type,".time")]
  clinical$event=clinical[,type]
  res.cut=surv_cutpoint(clinical,time = "time",event="event",variables = "score")
  res.cat=surv_categorize(res.cut)
  res.cat$score=factor(res.cat$score,levels = c("low","high"))
  p=unicox(surv=res.cat)
  p1=p$plot+labs(title = paste(type,database,"total"))
  
  
  model.1 <- coxph(Surv(time, event) ~ score, data = res.cat)
  cox_summary <- summary(model.1)
  HR <- round(exp(cox_summary$coefficients[, "coef"]), 3)
  CI_lower <- round(cox_summary$conf.int[, "lower .95"], 3)
  CI_upper <- round(cox_summary$conf.int[, "upper .95"], 3)
  p_value <- round(cox_summary$coefficients[, "Pr(>|z|)"], 3)
  
  # 构建标注文本（如：HR=1.892, 95%CI=1.234-2.890, p=0.0032）
  annotation_text <- paste0(
    "HR = ", HR, "\n",
    "95%CI = (", CI_lower, " - ", CI_upper, ")\n",
    "p = ", p_value
  )
  p1$plot <- p1$plot +
    annotate(
      "text",                  # 标注类型：文本
      x = 0, y = 0.3,        # 标注位置（根据你的图调整x/y）
      label = annotation_text, # 标注内容
      size = 3.5,              # 字体大小
      color = "black",         # 字体颜色
      hjust = 0                # 左对齐
    )
  print(p1)
  
  
  
  HRdf=data.frame(type="total",HR,CI_lower,CI_upper,p_value)
  
  
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
      res.cat$score=factor(res.cat$score,levels = c("low","high"))
      p=unicox(surv=res.cat)
      p1=p$plot+labs(title = paste(type,database,i,j))
      model.1 <- coxph(Surv(time, event) ~ score, data = res.cat)
      cox_summary <- summary(model.1)
      HR <- round(exp(cox_summary$coefficients[, "coef"]), 3)
      CI_lower <- round(cox_summary$conf.int[, "lower .95"], 3)
      CI_upper <- round(cox_summary$conf.int[, "upper .95"], 3)
      p_value <- round(cox_summary$coefficients[, "Pr(>|z|)"], 3)
      
      # 构建标注文本（如：HR=1.892, 95%CI=1.234-2.890, p=0.0032）
      annotation_text <- paste0(
        "HR = ", HR, "\n",
        "95%CI = (", CI_lower, " - ", CI_upper, ")\n",
        "p = ", p_value
      )
      p1$plot <- p1$plot +
        annotate(
          "text",                  # 标注类型：文本
          x = 0, y = 0.3,        # 标注位置（根据你的图调整x/y）
          label = annotation_text, # 标注内容
          size = 3.5,              # 字体大小
          color = "black",         # 字体颜色
          hjust = 0                # 左对齐
        )
      print(p1)
      
      HRdf=rbind(HRdf,data.frame(type=paste(i,j),HR,CI_lower,CI_upper,p_value))
    }
  }
  dev.off()
  return(HRdf)
}


showpathsurv.correctedforclinical=function(gsname,score,database){
  
  expdata=BRCA_list[[database]]$exp
  clinicaldata=BRCA_list[[database]]$clinical
  type="OS"
  # df=data.frame(exp=expdata[gsname,])
  df=data.frame(exp=score[,1])
  df=cbind(df,clinicaldata)
  
  
  library(survminer)
  library(survival)
  
  # pdf(paste(database,gsname,"OS.pdf",sep = "."))
  clinical=df
  if(database=="TCGA"){
    clinical=df[!df$normalORnot%in%c("normal"),]
  }
  clinical$score=clinical$exp
  
  # plist=list()
  # for(type in c(unique(do.call("c",survivaltype)))){}
  # if(!type%in%colnames(clinical)){next()}
  clinical$time=clinical[,paste0(type,".time")]
  clinical$event=clinical[,type]
  # res.cut=surv_cutpoint(clinical,time = "time",event="event",variables = "score")
  # res.cat=surv_categorize(res.cut)
  # res.cat$score=factor(res.cat$score,levels = c("low","high"))
  # p=unicox(surv=res.cat)
  # p1=p$plot+labs(title = paste(type,database,"total"))
  
  
  model.1 <- coxph(Surv(time, event) ~ score, data = clinical)
  if(database=="TCGA"){
    model.1 <- coxph(Surv(time, event) ~ score+age+stage+ER+PR+HER2, data = clinical)
  }else if(database=="metabric"){
    model.1 <- coxph(Surv(time, event) ~ score+age+stage+grade+LN+ER+PR+HER2+chemotherapy+hormonetherapy+radiotherapy, data = clinical)
  }else if(database=="SCANB"){
    model.1 <- coxph(Surv(time, event) ~ score+age+stage+grade+LN+ER+PR+HER2+TreatGroup, data = clinical)
  }else{
    print("ERROR, other database")
  }
  
  
  cox_summary <- summary(model.1)
  HR <- round(exp(cox_summary$coefficients[1, "coef"]), 3)
  CI_lower <- round(cox_summary$conf.int[1, "lower .95"], 3)
  CI_upper <- round(cox_summary$conf.int[1, "upper .95"], 3)
  p_value <- round(cox_summary$coefficients[1, "Pr(>|z|)"], 3)
  
  # # 构建标注文本（如：HR=1.892, 95%CI=1.234-2.890, p=0.0032）
  # annotation_text <- paste0(
  #   "HR = ", HR, "\n",
  #   "95%CI = (", CI_lower, " - ", CI_upper, ")\n",
  #   "p = ", p_value
  # )
  # p1$plot <- p1$plot +
  #   annotate(
  #     "text",                  # 标注类型：文本
  #     x = 0, y = 0.3,        # 标注位置（根据你的图调整x/y）
  #     label = annotation_text, # 标注内容
  #     size = 3.5,              # 字体大小
  #     color = "black",         # 字体颜色
  #     hjust = 0                # 左对齐
  #   )
  # print(p1)
  
  
  
  HRdf=data.frame(type="total",HR,CI_lower,CI_upper,p_value)
  
  
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
      # res.cut=surv_cutpoint(t,time = "time",event="event",variables = "score")
      # res.cat=surv_categorize(res.cut)
      # res.cat$score=factor(res.cat$score,levels = c("low","high"))
      # p=unicox(surv=res.cat)
      # p1=p$plot+labs(title = paste(type,database,i,j))
      
      model.1 <- coxph(Surv(time, event) ~ score, data = t)
      if(database=="TCGA"){
        model.1 <- coxph(Surv(time, event) ~ score+age+stage+ER+PR+HER2, data = t)
      }else if(database=="metabric"){
        model.1 <- coxph(Surv(time, event) ~ score+age+stage+grade+LN+ER+PR+HER2+chemotherapy+hormonetherapy+radiotherapy, data = t)
      }else if(database=="SCANB"){
        model.1 <- coxph(Surv(time, event) ~ score+age+stage+grade+LN+ER+PR+HER2+TreatGroup, data = t)
      }else{
        print("ERROR, other database")
      }
      
      cox_summary <- summary(model.1)
      HR <- round(exp(cox_summary$coefficients[1, "coef"]), 3)
      CI_lower <- round(cox_summary$conf.int[1, "lower .95"], 3)
      CI_upper <- round(cox_summary$conf.int[1, "upper .95"], 3)
      p_value <- round(cox_summary$coefficients[1, "Pr(>|z|)"], 3)
      
      # # 构建标注文本（如：HR=1.892, 95%CI=1.234-2.890, p=0.0032）
      # annotation_text <- paste0(
      #   "HR = ", HR, "\n",
      #   "95%CI = (", CI_lower, " - ", CI_upper, ")\n",
      #   "p = ", p_value
      # )
      # p1$plot <- p1$plot +
      #   annotate(
      #     "text",                  # 标注类型：文本
      #     x = 0, y = 0.3,        # 标注位置（根据你的图调整x/y）
      #     label = annotation_text, # 标注内容
      #     size = 3.5,              # 字体大小
      #     color = "black",         # 字体颜色
      #     hjust = 0                # 左对齐
      #   )
      # print(p1)
      
      HRdf=rbind(HRdf,data.frame(type=paste(i,j),HR,CI_lower,CI_upper,p_value))
    }
  }
  # dev.off()
  return(HRdf)
}


showpathsurv.medianthreshold=function(gsname,score,database){
  
  expdata=BRCA_list[[database]]$exp
  clinicaldata=BRCA_list[[database]]$clinical
  type="OS"
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
  
  # plist=list()
  # for(type in c(unique(do.call("c",survivaltype)))){}
  # if(!type%in%colnames(clinical)){next()}
  clinical$time=clinical[,paste0(type,".time")]
  clinical$event=clinical[,type]
  # res.cut=surv_cutpoint(clinical,time = "time",event="event",variables = "score")
  # res.cat=surv_categorize(res.cut)
  res.cat=clinical[,c("time","event","score")]
  medianvalue=median(res.cat$score)
  res.cat$score=plyr::mapvalues(res.cat$score>medianvalue,from = c(T,F),to=c("high","low"))
  res.cat$score=factor(res.cat$score,levels = c("low","high"))
  p=unicox(surv=res.cat)
  p1=p$plot+labs(title = paste(type,database,"total"))
  
  
  model.1 <- coxph(Surv(time, event) ~ score, data = res.cat)
  cox_summary <- summary(model.1)
  HR <- round(exp(cox_summary$coefficients[, "coef"]), 3)
  CI_lower <- round(cox_summary$conf.int[, "lower .95"], 3)
  CI_upper <- round(cox_summary$conf.int[, "upper .95"], 3)
  p_value <- round(cox_summary$coefficients[, "Pr(>|z|)"], 3)
  
  # 构建标注文本（如：HR=1.892, 95%CI=1.234-2.890, p=0.0032）
  annotation_text <- paste0(
    "HR = ", HR, "\n",
    "95%CI = (", CI_lower, " - ", CI_upper, ")\n",
    "p = ", p_value
  )
  p1$plot <- p1$plot +
    annotate(
      "text",                  # 标注类型：文本
      x = 0, y = 0.3,        # 标注位置（根据你的图调整x/y）
      label = annotation_text, # 标注内容
      size = 3.5,              # 字体大小
      color = "black",         # 字体颜色
      hjust = 0                # 左对齐
    )
  print(p1)
  
  
  
  HRdf=data.frame(type="total",HR,CI_lower,CI_upper,p_value)
  
  
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
      # res.cut=surv_cutpoint(t,time = "time",event="event",variables = "score")
      # res.cat=surv_categorize(res.cut)
      res.cat=t[,c("time","event","score")]
      medianvalue=median(res.cat$score)
      res.cat$score=plyr::mapvalues(res.cat$score>medianvalue,from = c(T,F),to=c("high","low"))
      res.cat$score=factor(res.cat$score,levels = c("low","high"))
      p=unicox(surv=res.cat)
      p1=p$plot+labs(title = paste(type,database,i,j))
      model.1 <- coxph(Surv(time, event) ~ score, data = res.cat)
      cox_summary <- summary(model.1)
      HR <- round(exp(cox_summary$coefficients[, "coef"]), 3)
      CI_lower <- round(cox_summary$conf.int[, "lower .95"], 3)
      CI_upper <- round(cox_summary$conf.int[, "upper .95"], 3)
      p_value <- round(cox_summary$coefficients[, "Pr(>|z|)"], 3)
      
      # 构建标注文本（如：HR=1.892, 95%CI=1.234-2.890, p=0.0032）
      annotation_text <- paste0(
        "HR = ", HR, "\n",
        "95%CI = (", CI_lower, " - ", CI_upper, ")\n",
        "p = ", p_value
      )
      p1$plot <- p1$plot +
        annotate(
          "text",                  # 标注类型：文本
          x = 0, y = 0.3,        # 标注位置（根据你的图调整x/y）
          label = annotation_text, # 标注内容
          size = 3.5,              # 字体大小
          color = "black",         # 字体颜色
          hjust = 0                # 左对齐
        )
      print(p1)
      
      HRdf=rbind(HRdf,data.frame(type=paste(i,j),HR,CI_lower,CI_upper,p_value))
    }
  }
  dev.off()
  return(HRdf)
}



showpathsurv.continuous=function(gsname,score,database){
  
  expdata=BRCA_list[[database]]$exp
  clinicaldata=BRCA_list[[database]]$clinical
  type="OS"
  # df=data.frame(exp=expdata[gsname,])
  df=data.frame(exp=score[,1])
  df=cbind(df,clinicaldata)
  
  
  library(survminer)
  library(survival)
  
  # pdf(paste(database,gsname,"OS.pdf",sep = "."))
  clinical=df
  if(database=="TCGA"){
    clinical=df[!df$normalORnot%in%c("normal"),]
  }
  clinical$score=clinical$exp
  
  # plist=list()
  # for(type in c(unique(do.call("c",survivaltype)))){}
  # if(!type%in%colnames(clinical)){next()}
  clinical$time=clinical[,paste0(type,".time")]
  clinical$event=clinical[,type]
  # res.cut=surv_cutpoint(clinical,time = "time",event="event",variables = "score")
  # res.cat=surv_categorize(res.cut)
  res.cat=clinical[,c("time","event","score")]
  # medianvalue=median(res.cat$score)
  # res.cat$score=plyr::mapvalues(res.cat$score>medianvalue,from = c(T,F),to=c("high","low"))
  # res.cat$score=factor(res.cat$score,levels = c("low","high"))
  # p=unicox(surv=res.cat)
  # p1=p$plot+labs(title = paste(type,database,"total"))
  
  
  model.1 <- coxph(Surv(time, event) ~ score, data = res.cat)
  cox_summary <- summary(model.1)
  HR <- round(exp(cox_summary$coefficients[, "coef"]), 3)
  CI_lower <- round(cox_summary$conf.int[, "lower .95"], 3)
  CI_upper <- round(cox_summary$conf.int[, "upper .95"], 3)
  p_value <- round(cox_summary$coefficients[, "Pr(>|z|)"], 3)
  
  # 构建标注文本（如：HR=1.892, 95%CI=1.234-2.890, p=0.0032）
  annotation_text <- paste0(
    "HR = ", HR, "\n",
    "95%CI = (", CI_lower, " - ", CI_upper, ")\n",
    "p = ", p_value
  )
  # p1$plot <- p1$plot +
  #   annotate(
  #     "text",                  # 标注类型：文本
  #     x = 0, y = 0.3,        # 标注位置（根据你的图调整x/y）
  #     label = annotation_text, # 标注内容
  #     size = 3.5,              # 字体大小
  #     color = "black",         # 字体颜色
  #     hjust = 0                # 左对齐
  #   )
  # print(p1)
  
  
  
  HRdf=data.frame(type="total",HR,CI_lower,CI_upper,p_value)
  
  
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
      # res.cut=surv_cutpoint(t,time = "time",event="event",variables = "score")
      # res.cat=surv_categorize(res.cut)
      res.cat=t[,c("time","event","score")]
      # medianvalue=median(res.cat$score)
      # res.cat$score=plyr::mapvalues(res.cat$score>medianvalue,from = c(T,F),to=c("high","low"))
      # res.cat$score=factor(res.cat$score,levels = c("low","high"))
      # p=unicox(surv=res.cat)
      # p1=p$plot+labs(title = paste(type,database,i,j))
      model.1 <- coxph(Surv(time, event) ~ score, data = res.cat)
      cox_summary <- summary(model.1)
      HR <- round(exp(cox_summary$coefficients[, "coef"]), 3)
      CI_lower <- round(cox_summary$conf.int[, "lower .95"], 3)
      CI_upper <- round(cox_summary$conf.int[, "upper .95"], 3)
      p_value <- round(cox_summary$coefficients[, "Pr(>|z|)"], 3)
      
      # 构建标注文本（如：HR=1.892, 95%CI=1.234-2.890, p=0.0032）
      annotation_text <- paste0(
        "HR = ", HR, "\n",
        "95%CI = (", CI_lower, " - ", CI_upper, ")\n",
        "p = ", p_value
      )
      # p1$plot <- p1$plot +
      #   annotate(
      #     "text",                  # 标注类型：文本
      #     x = 0, y = 0.3,        # 标注位置（根据你的图调整x/y）
      #     label = annotation_text, # 标注内容
      #     size = 3.5,              # 字体大小
      #     color = "black",         # 字体颜色
      #     hjust = 0                # 左对齐
      #   )
      # print(p1)
      
      HRdf=rbind(HRdf,data.frame(type=paste(i,j),HR,CI_lower,CI_upper,p_value))
    }
  }
  # dev.off()
  return(HRdf)
}

showpredictedRisksurv=function(gsname,score,database){
  
  expdata=BRCA_list[[database]]$exp
  clinicaldata=BRCA_list[[database]]$clinical
  type="OS"
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
  
  # plist=list()
  # for(type in c(unique(do.call("c",survivaltype)))){}
  # if(!type%in%colnames(clinical)){next()}
  clinical$time=clinical[,paste0(type,".time")]
  clinical$event=clinical[,type]
  # res.cut=surv_cutpoint(clinical,time = "time",event="event",variables = "score")
  # res.cat=surv_categorize(res.cut)
  res.cat=clinical[,c("time","event","score")]
  res.cat$score=plyr::mapvalues(res.cat$score,from=c(0,1),to=c("low","high"))
  res.cat$score=factor(res.cat$score,levels = c("low","high"))
  p=unicox(surv=res.cat)
  p1=p$plot+labs(title = paste(type,database,"total"))
  
  
  model.1 <- coxph(Surv(time, event) ~ score, data = res.cat)
  cox_summary <- summary(model.1)
  HR <- round(exp(cox_summary$coefficients[, "coef"]), 3)
  CI_lower <- round(cox_summary$conf.int[, "lower .95"], 3)
  CI_upper <- round(cox_summary$conf.int[, "upper .95"], 3)
  p_value <- round(cox_summary$coefficients[, "Pr(>|z|)"], 3)
  
  # 构建标注文本（如：HR=1.892, 95%CI=1.234-2.890, p=0.0032）
  annotation_text <- paste0(
    "HR = ", HR, "\n",
    "95%CI = (", CI_lower, " - ", CI_upper, ")\n",
    "p = ", p_value
  )
  p1$plot <- p1$plot +
    annotate(
      "text",                  # 标注类型：文本
      x = 0, y = 0.3,        # 标注位置（根据你的图调整x/y）
      label = annotation_text, # 标注内容
      size = 3.5,              # 字体大小
      color = "black",         # 字体颜色
      hjust = 0                # 左对齐
    )
  print(p1)
  
  
  
  HRdf=data.frame(type="total",HR,CI_lower,CI_upper,p_value)
  
  
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
      # res.cut=surv_cutpoint(t,time = "time",event="event",variables = "score")
      # res.cat=surv_categorize(res.cut)
      res.cat=t[,c("time","event","score")]
      res.cat$score=plyr::mapvalues(res.cat$score,from=c(0,1),to=c("low","high"))
      if(length(unique(res.cat$score))<2){next()}
      res.cat$score=factor(res.cat$score,levels = c("low","high"))
      
      p=unicox(surv=res.cat)
      p1=p$plot+labs(title = paste(type,database,i,j))
      model.1 <- coxph(Surv(time, event) ~ score, data = res.cat)
      cox_summary <- summary(model.1)
      HR <- round(exp(cox_summary$coefficients[, "coef"]), 3)
      CI_lower <- round(cox_summary$conf.int[, "lower .95"], 3)
      CI_upper <- round(cox_summary$conf.int[, "upper .95"], 3)
      p_value <- round(cox_summary$coefficients[, "Pr(>|z|)"], 3)
      
      # 构建标注文本（如：HR=1.892, 95%CI=1.234-2.890, p=0.0032）
      annotation_text <- paste0(
        "HR = ", HR, "\n",
        "95%CI = (", CI_lower, " - ", CI_upper, ")\n",
        "p = ", p_value
      )
      p1$plot <- p1$plot +
        annotate(
          "text",                  # 标注类型：文本
          x = 0, y = 0.3,        # 标注位置（根据你的图调整x/y）
          label = annotation_text, # 标注内容
          size = 3.5,              # 字体大小
          color = "black",         # 字体颜色
          hjust = 0                # 左对齐
        )
      print(p1)
      
      HRdf=rbind(HRdf,data.frame(type=paste(i,j),HR,CI_lower,CI_upper,p_value))
    }
  }
  dev.off()
  return(HRdf)
}

