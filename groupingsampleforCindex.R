## sample grouping for C-index, time-dependent AUC..
gs=genelistintersect.list

grouplist=list()

for(database in c("TCGA","metabric","SCANB")){
  
  exp=BRCA_list[[database]]$exp
  scoreMat=sapply(gs,function(x){apply(as.matrix(exp[match(x,rownames(exp),nomatch = 0),]),2,median)})
  
  for(gsname in names(gs)){
    print(paste(database,gsname))
    score=scoreMat[,gsname]
    
    clinicaldata=BRCA_list[[database]]$clinical
    type="OS"
    df=data.frame(exp=score)
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
    
    if(grepl("total",gsname)){
      res.cut=surv_cutpoint(clinical,time = "time",event="event",variables = "score")
      res.cat=surv_categorize(res.cut)
      res.cat$score=factor(res.cat$score,levels = c("low","high"))
      clinical$score=res.cat$score
      grouplist[[database]][[gsname]]=clinical
    }else{
      clinical=clinical[clinical$PAM50%in%(unlist(strsplit(gsname,"\\."))[2]),]
      res.cut=surv_cutpoint(clinical,time = "time",event="event",variables = "score")
      res.cat=surv_categorize(res.cut)
      res.cat$score=factor(res.cat$score,levels = c("low","high"))
      clinical$score=res.cat$score
      grouplist[[database]][[gsname]]=clinical
    }    
  }
}
save(grouplist,file = "grouplist.RData")

grouplist2=lapply(grouplist,function(x){x[grep("poor",names(x))]})

risklist2=lapply(names(risklist),function(x){a=BRCA_list[[x]]$clinical;a$sig.gene70=risklist[[x]]$sig.gene70;a$sig.endoPredict=risklist[[x]]$sig.endoPredict;return(a)})
names(risklist2)=names(risklist)

risklist2=lapply(risklist2,function(x){x$sig.gene70=plyr::mapvalues(x$sig.gene70,from = c(0,1),to=c("low","high"));x$sig.gene70=factor(x$sig.gene70,levels = c("low","high"));
x$sig.endoPredict=plyr::mapvalues(x$sig.endoPredict,from = c(0,1),to=c("low","high"));x$sig.endoPredict=factor(x$sig.endoPredict,levels = c("low","high"));
return(x)})



evaluationlist=list()
library(survcomp)
library(timeROC)
library(survIDINRI)
##C-index

Cindexlist=list()
for(database in names(grouplist2)){
  t1=grouplist2[[database]]
  
  for(gsname in names(t1)){
    print(paste(database,gsname))
    t2=t1[[gsname]]
    t2$gene70=risklist2[[database]]$sig.gene70[match(t2$ID,risklist2[[database]]$ID)]
    t2$endoPredict=risklist2[[database]]$sig.endoPredict[match(t2$ID,risklist2[[database]]$ID)]
    
    if(grepl("total",gsname)){
      print("step_step total")
    }else{
      t2=t2[t2$PAM50%in%(unlist(strsplit(gsname,"\\."))[2]),]
    }
    t2=t2[!is.na(t2$time),]
    t2$score=as.numeric(t2$score)
    t2$gene70=as.numeric(t2$gene70)
    t2$endoPredict=as.numeric(t2$endoPredict)
    
    surv_obj <- with(t2, Surv(time, event))
    
    
    data=t2
    
    # 计算每个模型的C指数
    cindex_my <- concordance.index(x = data$score, 
                                   surv.time = data$time, 
                                   surv.event = data$event,
                                   method = "noether") # Noether's SE提供置信区间
    
    cindex_mp <- concordance.index(x = data$gene70, 
                                   surv.time = data$time, 
                                   surv.event = data$event,
                                   method = "noether")
    
    cindex_ep <- concordance.index(x = data$endoPredict, 
                                   surv.time = data$time, 
                                   surv.event = data$event,
                                   method = "noether")
    Cindexlist[[database]][[gsname]]=c(cindex_my$c.index,cindex_mp$c.index,cindex_ep$c.index)
    # 查看结果
    # cindex_my$c.index  # C指数值
    # cindex_my$se        # 标准误
    # cindex_my$lower     # 置信区间下限
    # cindex_my$upper     # 置信区间上限
    
    # 两两比较：检验您的签名是否显著优于MammaPrint
    # comp_my_vs_mp <- cindex.comp(cindex_my, cindex_mp)
    # print(comp_my_vs_mp$p.value) # P值 < 0.05 表示您的模型显著更好
    
    # 比较您的签名 vs EndoPredict
    # comp_my_vs_ep <- cindex.comp(cindex_my, cindex_ep)
    # print(comp_my_vs_mp$p.value)
  }
}
evaluationlist[["Cindex"]]=Cindexlist
Cindex=do.call("rbind",do.call("c",Cindexlist))
colnames(Cindex)=c("thispaper","gene70","endoPredict")


##time-dependent AUC
AUClist=list()
for(database in names(grouplist2)){
  t1=grouplist2[[database]]
  
  for(gsname in names(t1)){
    print(paste(database,gsname))
    t2=t1[[gsname]]
    t2$gene70=risklist2[[database]]$sig.gene70[match(t2$ID,risklist2[[database]]$ID)]
    t2$endoPredict=risklist2[[database]]$sig.endoPredict[match(t2$ID,risklist2[[database]]$ID)]
    
    if(grepl("total",gsname)){
      print("step_step total")
    }else{
      t2=t2[t2$PAM50%in%(unlist(strsplit(gsname,"\\."))[2]),]
    }
    t2=t2[!is.na(t2$time),]
    t2$score=as.numeric(t2$score)
    t2$gene70=as.numeric(t2$gene70)
    t2$endoPredict=as.numeric(t2$endoPredict)
    
    
    data=t2
    
    # 定义感兴趣的时间点（例如：1年，3年，5年）
    times <- c(1,3,5) # 假设时间单位是天
    
    # 计算时间依赖AUC
    roc_results <- timeROC(T = data$time,
                           delta = data$event,
                           marker = data$score,
                           cause = 1, # 感兴趣的事件
                           times = times,
                           iid = TRUE) # 需要iid来计算标准误和置信区间
    
    # 查看特定时间点的AUC，例如5年
    # print(roc_results$AUC) 
    
    # 计算置信区间
    # confint(roc_results, level = 0.95)
    
    # 如果要比较两个模型的AUC（例如，您的模型 vs MammaPrint）
    # 需要为每个模型分别计算roc_results，然后使用survcomp包的函数
    roc_mp <- timeROC(T = data$time, delta = data$event, 
                      marker = data$gene70, cause = 1,
                      times = times, iid = TRUE)
    roc_ep <- timeROC(T = data$time, delta = data$event, 
                      marker = data$endoPredict, cause = 1,
                      times = times, iid = TRUE)
    # 使用 compareC 包或 survcomp 中的函数进行AUC比较
    # library(compareC)
    # compareC(data$time, data$event, data$score, data$gene70)
    
    a=roc_results$AUC;a=a[c("t=1","t=3","t=5")]
    b=gene70=roc_mp$AUC;b=b[c("t=1","t=3","t=5")]
    cc=endoPredict=roc_ep$AUC;cc=cc[c("t=1","t=3","t=5")]
    
    
    d=data.frame(thispaper=roc_results$AUC,gene70=roc_mp$AUC,endoPredict=roc_ep$AUC)
    AUClist[[database]][[gsname]]=d

  }
}
evaluationlist[["AUC"]]=AUClist
AUCs=do.call("rbind",do.call("c",AUClist))
AUCs$interesttime=sapply(rownames(AUCs),function(x){unlist(strsplit(x,".poor."))[2]})
AUCs$geneset=sapply(rownames(AUCs),function(x){paste(unlist(strsplit(x,"\\."))[2:3],collapse = ".")})
AUCs$database=sapply(rownames(AUCs),function(x){unlist(strsplit(x,"\\."))[1]})

# save()

# 拟合一个仅包含您的签名的Cox模型  ## risk已经是二分类变量，无法进行calibration分析
# cox_my <- coxph(Surv(time, event) ~ score, data = data)
# 
# # 获取基础风险函数，并计算每个患者在特定时间点（如5年）的预测生存概率
# # 使用 survfit 函数
# surv_prob_my <- survfit(cox_my, newdata = data)
# # 提取5年生存概率
# # 首先需要找到最接近5年的时间点
# time_point <- 5
# time_index <- findInterval(time_point, summary(surv_prob_my)$time)
# predicted_surv_5y <- summary(surv_prob_my, times = time_point)$surv
# 
# # 将患者按预测风险分组（例如，分成4-5组）
# data$pred_risk_group <- cut(predicted_surv_5y, 
#                             breaks = quantile(predicted_surv_5y, probs = c(0,median(unique(as.numeric(predicted_surv_5y))))), 
#                             include.lowest = TRUE)
# 
# # 计算每组观察到的生存概率（使用Kaplan-Meier法）
# observed_surv <- aggregate(data$time, by = list(data$pred_risk_group), 
#                            FUN = function(t) {
#                              km_fit <- survfit(Surv(t, data$event[data$pred_risk_group == group]) ~ 1)
#                              # 提取时间点最接近time_point的生存概率
#                              summary(km_fit, times = time_point)$surv
#                            })
# 
# # 绘制校准图
# plot(predicted_surv_5y, observed_surv, xlab = "Predicted 5-Year Survival", 
#      ylab = "Observed 5-Year Survival", main = "Calibration Plot")
# abline(0, 1, col = "red", lty = 2) # 添加对角线作为完美校准的参考线

##time-dependent NRI IDI
NRIlist=list()
IDIlist=list()
for(database in names(grouplist2)){
  t1=grouplist2[[database]]
  
  for(gsname in names(t1)){
    print(paste(database,gsname))
    t2=t1[[gsname]]
    t2$gene70=risklist2[[database]]$sig.gene70[match(t2$ID,risklist2[[database]]$ID)]
    t2$endoPredict=risklist2[[database]]$sig.endoPredict[match(t2$ID,risklist2[[database]]$ID)]
    
    if(grepl("total",gsname)){
      print("step_step total")
    }else{
      t2=t2[t2$PAM50%in%(unlist(strsplit(gsname,"\\."))[2]),]
    }
    t2=t2[!is.na(t2$time),]
    t2$score=as.numeric(t2$score)
    t2$gene70=as.numeric(t2$gene70)
    t2$endoPredict=as.numeric(t2$endoPredict)
    
    
    data=t2
    
    # survIDINRI 要求输入协变量矩阵
    # 将您的风险评分作为唯一的协变量（或者与其他临床变量组合）
    # 模型0 (旧模型): 例如，只包含临床变量
    # 模型1 (新模型): 在模型0的基础上 + 您的签名
    
    # 这里我们简化比较：假设旧模型是MammaPrint，新模型是您的签名
    # 您需要先构建一个包含两个模型风险评分的矩阵，或者拟合Cox模型得到线性预测值
    
    # 示例：比较两个单变量模型（您的签名 vs MammaPrint）
    # 首先，需要为每个模型拟合Cox并提取线性预测值（lp）
    cox_my <- coxph(Surv(time, event) ~ score, data = data)
    cox_mp <- coxph(Surv(time, event) ~ gene70, data = data)
    cox_ep <- coxph(Surv(time, event) ~ endoPredict, data = data)
    
    # 提取线性预测值作为风险评分
    data$lp_my <- predict(cox_my, type = "lp")
    data$lp_mp <- predict(cox_mp, type = "lp")
    data$lp_ep <- predict(cox_ep, type = "lp")
    
    # 准备 survIDINRI 的输入
    indata <- data.frame(time = data$time, event = data$event)
    covs0 <- as.matrix(data$lp_mp)   # 基础模型 (MammaPrint)
    covs00 <- as.matrix(data$lp_ep)   # 基础模型 (MammaPrint)
    covs1 <- as.matrix(data$lp_my)   # 新模型 (您的签名)
    
    # 定义评估的时间点，例如5年
    t0 <- 5
    
    # 运行IDI/NRI计算
    # npert 是扰动重抽样的次数，用于计算置信区间，可适当增加以提高精度
    mp=matrix(Inf,ncol = 4,nrow = 3)
    colnames(mp)=c("Est.","Lower","Upper","p-value")
    rownames(mp)=c(paste0("M",1:3))
    
    np=mp
    
    if(length(unique(data$gene70))==2){
      idi_result <- IDI.INF(indata = indata, 
                            covs0 = covs0, 
                            covs1 = covs1, 
                            t0 = t0, 
                            npert = 300, # 可以设置为1000
                            seed1 = 1234) # 设置种子以便结果可重复
      mp=IDI.INF.OUT(idi_result)
    }
    
    if(length(unique(data$endoPredict))==2){
      idi_result2 <- IDI.INF(indata = indata, 
                             covs0 = covs00, 
                             covs1 = covs1, 
                             t0 = t0, 
                             npert = 300, # 可以设置为1000
                             seed1 = 1234) # 设置种子以便结果可重复
      ep=IDI.INF.OUT(idi_result2)
    }
    
    # 输出结果
    
    
    # NRI：survIDINRI输出的m2是连续NRI，可以看作是事件组中风险评分提高的比例与非事件组中风险评分降低的比例之差。
    # 例如，NRI = 0.15 (95% CI: 0.05, 0.25) 表示，相对于MammaPrint，您的模型净正确地将15%的患者重新分类到了更准确的风险类别中。
    # 重点关注置信区间是否包含0。如果不包含0，说明改善是显著的。
    # IDI：m1是IDI值，衡量的是模型对事件组和非事件组平均预测风险差的改进。
    # IDI = 0.03 (95% CI: 0.01, 0.05) 表示，您的模型使事件组和非事件组的平均预测风险差异提高了3个百分点
    # 。即使这个数字看起来很小，在临床上也可能有重要意义，关键在于其置信区间是否包含0。
    # 可选：绘制图形展示结果
    # IDI.INF.GRAPH(idi_result)
    a=paste0(mp[,1],"(95%CI:",mp[,2],"-",mp[,3],")")
    b=paste0(ep[,1],"(95%CI:",ep[,2],"-",ep[,3],")")
    cc=data.frame(gene70=a[1:2],endoPredict=b[1:2])
    rownames(cc)=c("IDI","NRI")
    
    NRIlist[[database]][[gsname]]=cc["NRI",]
    IDIlist[[database]][[gsname]]=cc["IDI",]
    
  }
}

evaluationlist[["NRI"]]=NRIlist
evaluationlist[["IDI"]]=IDIlist
save(evaluationlist,file = "evaluationlist.RData")


evaluationlist2=lapply(evaluationlist,function(x){as.data.frame(do.call("rbind",do.call("c",x)))})
colnames(evaluationlist2$Cindex)=c("thispaper","gene70","endoPredict")

evaluationlist2=lapply(evaluationlist2,function(x){x$database=sapply(rownames(x),function(y){unlist(strsplit(y,"\\."))[1]});
x$geneset=sapply(rownames(x),function(y){paste(unlist(strsplit(y,"\\."))[2:3],collapse = ".")});
return(x)})

for(i in c("Cindex","AUC")){
  t=evaluationlist2[[i]]
  for(j in c("thispaper","gene70","endoPredict")){
    t[,j]=signif(t[,j],digits = 3)
  }
  evaluationlist2[[i]]=t
}

evaluationlist2$AUC$time=sapply(rownames(evaluationlist2$AUC),function(x){unlist(strsplit(x,"poor."))[2]})

writexl::write_xlsx(evaluationlist2,path = "evaluation.xlsx")



mat=evaluationlist2$Cindex[,1:3]
rownames(mat)=gsub(".poor","",rownames(mat))
max_mark <- matrix("", nrow = nrow(mat), ncol = 3)
for(i in 1:nrow(mat)) {
  row_data <- mat[i, ]
  if(all(is.na(row_data))) next  # 如果整行都是NA，跳过
  max_idx <- which.max(row_data)
  max_mark[i, max_idx] <- "*"  # 使用星星标记最大值
}

# 创建显示用的矩阵（将数值转换为字符，加上标记）
display_mat <- matrix("", nrow = nrow(mat), ncol = 3)
for(i in 1:nrow(mat)) {
  for(j in 1:ncol(mat)) {
    if(!is.na(mat[i, j])) {
      display_mat[i, j] <- paste0(round(mat[i, j], 2), max_mark[i, j])
    } else {
      display_mat[i, j] <- "NA"
    }
  }
}
colnames(display_mat) <- colnames(mat)
rownames(display_mat) <- rownames(mat)

# 使用pheatmap绘制
pdf("Cindex.heatmap.pdf",width = 4)
pheatmap::pheatmap(mat,cluster_rows = F,cluster_cols = F,
                   display_numbers = display_mat,  # 显示带有标记的数字
                   number_color = "black",
                   fontsize_number = 10,
                   main = "C-index Heatmap with Max Values Marked (*)",
                   # cellheight = 20,cellwidth = 30,
                   na_col = "lightgray")  # NA值显示为灰色
dev.off()