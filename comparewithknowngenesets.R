##compare with previous genesets
load("X:/analyze_data/FEdata/liuxiaoqin/20250815_BRCA_db/BRCA_list.RData")


library(genefu)
risklist=list()
data("sig.endoPredict")
data("sig.gene70")
data("sig.oncotypedx")
for(database in names(BRCA_list)[1:3]){
  print(database)
  exp=BRCA_list[[database]]$exp
  data0=t(exp)
  annot=data.frame(probe=colnames(data0))# EntrezGene.ID)
  annot$EntrezGene.ID=gene2id$ENTREZID[match(annot$probe,gene2id$SYMBOL)]
  annot=annot[!duplicated(annot$EntrezGene.ID),]
  annot=annot[!is.na(annot$EntrezGene.ID),]
  annot=annot[annot$EntrezGene.ID%in%sig.gene70$EntrezGene.ID,]
  
  data=data0[,annot$probe]
  colnames(data)=sig.gene70$probe[match(annot$EntrezGene.ID[match(colnames(data),annot$probe)],sig.gene70$EntrezGene.ID)]
  results <- gene70(
    data = data,          # 你的表达矩阵
    annot = annot,        # 你的探针注释，包含 EntrezGene.ID
    # do.mapping = TRUE,         # 启用跨平台基因映射
    std = "none",             # 标准化方法，通常选择 "scale"（均值0，标准差1）
    verbose = TRUE             # 打印处理信息，便于调试
  )
  risklist[[database]][["sig.gene70"]]=results$risk
  
  data=data0[,colnames(data0)%in%sig.endoPredict$symbol]
  
  colnames(data)=sig.endoPredict$probe.affy[match(colnames(data),sig.endoPredict$symbol)]
  annot=data.frame(probe=colnames(data))
  # annot=annot[annot$NCBI.gene.symbol%in%sig.endoPredict$symbol,]
  annot$EntrezGene.ID=sig.endoPredict$EntrezGene.ID[match(annot$probe,sig.endoPredict$probe.affy)]
  
  results <- endoPredict(
    data = data,          # 你的表达矩阵
    annot = annot,        # 你的探针注释，必须包含 EntrezGene.ID
    # do.mapping = TRUE,         # 强烈建议启用基因映射
    verbose = TRUE             # 打印处理信息，便于调试
  )
  risklist[[database]][["sig.endoPredict"]]=results$risk
  
  # sig.oncotypedx$probe.affy[(nrow(sig.oncotypedx)-2):nrow(sig.oncotypedx)]=paste0(1:3,"_at")
  # sig.oncotypedx$probe.affy=gsub("/",".",sig.oncotypedx$probe.affy)
  # sig.oncotypedx$probe.affy=gsub("-",".",sig.oncotypedx$probe.affy)
  # data=data0[,colnames(data0)%in%sig.oncotypedx$symbol]
  # # data=data[,!colnames(data)%in%c("RPLP0","GUSB","TFRC","GAPDH")]
  # colnames(data)=sig.oncotypedx$probe.affy[match(colnames(data),sig.oncotypedx$symbol)]
  # annot=data.frame(probe.affy=colnames(data))
  # annot$EntrezGene.ID=sig.oncotypedx$EntrezGene.ID[match(annot$probe.affy,sig.oncotypedx$probe.affy)]
  # results=oncotypedx(data=data, annot=annot,do.mapping = F,do.scaling = F, verbose = T)#do.mapping = T,
  
}
a=do.call("c",risklist)
save(risklist,file = "risklist.RData")


HRlist.pre=list()
for(database in names(risklist)){
  for(sigtype in c("sig.gene70","sig.endoPredict")){
    print(paste(database,sigtype))
    t=showpredictedRisksurv(gsname=sigtype,score=data.frame(score=risklist[[database]][[sigtype]]),database=database)
    t$database=database
    t$gsname=sigtype
    HRlist.pre[[database]][[sigtype]]=t
  }
}
# gsname=sigtype
# score=data.frame(score=risklist[[database]][[sigtype]])
# database=database

HRpred=do.call("rbind",do.call("c",HRlist.pre))
HRpred=HRpred[!HRpred$type%in%"PAM50 claudin-low",]
HRpred$gsname=paste0(HRpred$gsname,".poor")
# t1=grep("total",HRpred$type)
# newHRpred=HRpred[1,]
# for(i in t1){
#   thisgroup=HRpred[i:(i+5),]
#   m=matrix(NA,nrow = 1,ncol = ncol(HRpred))
#   m[1,1]=HRpred$gsname[i]
#   colnames(m)=colnames(HRpred)
#   thisgroup=rbind(m,thisgroup)
#   newHRpred=rbind(newHRpred,thisgroup)
# }
# newHRpred=newHRpred[-1,]
good=HRpred[grepl("good",HRpred$gsname),]
poor=HRpred[grepl("poor",HRpred$gsname),]
HRpredlist=lapply(c("good","poor"),function(x){HRpred[grepl(x,HRpred$gsname),]})
names(HRpredlist)=c("good","poor")
HRpredlist=HRpredlist[2]
HRpredlist2=lapply(HRpredlist,function(x){a=lapply(unique(x$database),function(y){x[x$database%in%y,]});names(a)=unique(x$database);return(a)})
# HRpredlist3=lapply(HRpredlist2,function(x){a=do.call("cbind",lapply(x,function(y){y[,c("HR","CI_lower","CI_upper","p_value","type","gsname")]}));return(a)
# })#a$type=x$TCGA$type;a$gsname=x$TCGA$gsname;return(a)
# HRpredlist3=lapply(HRpredlist2,function(x){do.call("cbind",x)})

tmp=HRpredlist2$poor
tmp=lapply(tmp,function(x){x$comb=paste(x$type,x$gsname);return(x)})
tmp2=merge(tmp$TCGA,tmp$metabric,by="comb",all = T)
tmp2=merge(tmp2,tmp$SCANB,by="comb",all=T)
a=apply(tmp2[,c("type.x","type.y","type")],1,function(x){unique(x[!is.na(x)])})
b=apply(tmp2[,c("gsname.x","gsname.y","gsname")],1,function(x){unique(x[!is.na(x)])})
tmp2$type=a
tmp2$gsname=b
tmp2=tmp2[,c("type","gsname","HR.x","CI_lower.x","CI_upper.x","p_value.x","HR.y","CI_lower.y","CI_upper.y","p_value.y","HR","CI_lower","CI_upper","p_value")]
colnames(tmp2)[3:14]=c(paste0(rep(c("TCGA.","metabric.","SCANB."),each=4),c("HR","CI_lower","CI_upper","p_value")))
tmp2$type=factor(tmp2$type,levels = c("total",paste0("PAM50 ",c("LumA","LumB","Her2","Basal","Normal"))))
tmp2$gsname=factor(tmp2$gsname,levels = c("sig.gene70.poor","sig.endoPredict.poor"))
tmp2=tmp2[order(tmp2$gsname,tmp2$type),]
HRpredlist3=list(poor=tmp2)

# HRpredlist3=lapply(HRpredlist3,function(x){x[,c(13:14,1:12)]})
# t=do.call("cbind",lapply(goodlist,function(x){x[,c("HRpred","CI_lower","CI_upper","p_value")]}))
# t$type=goodlist$TCGA$type
# t$gsname=goodlist$TCGA$gsname
HRpredlist4=lapply(HRpredlist3,function(x){t1=grep("total",x$type);

newHRpred=x[1,]
for(i in t1){
  thisgroup=x[i:(i+5),]
  m=matrix(NA,nrow = 1,ncol = ncol(x))
  m[1,1]=as.character(x$gsname[i])
  colnames(m)=colnames(x)
  thisgroup=rbind(m,thisgroup)
  newHRpred=rbind(newHRpred,thisgroup)
}
newHRpred=newHRpred[-1,]

})

writexl::write_xlsx(HRpred,path = "tmp.xlsx")


HRpredlist4=lapply(HRpredlist4,function(x){x$blank=paste(rep(" ", 20), collapse = " ");return(x)})
for(i in names(HRpredlist4)){
  t=HRpredlist4[[i]]
  for(j in c("TCGA.HR","TCGA.CI_lower","TCGA.CI_upper","TCGA.p_value","metabric.HR","metabric.CI_lower","metabric.CI_upper","metabric.p_value","SCANB.HR","SCANB.CI_lower","SCANB.CI_upper","SCANB.p_value")){
    t[,j]=as.numeric(t[,j])
  }
  HRpredlist4[[i]]=t
}

a=HRpredlist4$poor
for(i in c("TCGA","metabric","SCANB")){
  print(i)
  tmp=a[,paste0(i,".CI_upper")]
  tmp=tmp[!is.na(tmp)]
  tmp=tmp[!is.infinite(tmp)]
  a[is.infinite(a[,paste0(i,".CI_upper")]),paste0(i,".CI_upper")]=max(tmp)
}

# 绘制森林图
library(forestploter)
p <- forest(
  a[,c(1, 15)],  # 选择要在森林图中使用的数据列
  est = list(a$TCGA.HR,a$metabric.HR,a$SCANB.HR),
  lower = list(a$TCGA.CI_lower,a$metabric.CI_lower,a$SCANB.CI_lower),
  upper = list(a$TCGA.CI_upper,a$metabric.CI_upper,a$SCANB.CI_upper),
  ci_column = c(2),         # 指定CI列
  ref_line = 1,                # 添加参考线
  vert_line = c(0.5, 2),       # 添加垂直线
  nudge_y = 0.2,xlim = c(0,4)#,               # 垂直调整标签位置
  # theme = tm
)# 应用自定义主题
pdf("pred.forest.pdf")
p
dev.off()

HRpredlist5=list()
for(i in names(HRpredlist4)){
  t=HRpredlist4[[i]]
  t2=t[1,]
  for(j in names(gs)){
    tmp=t[t$gsname%in%j,]
    tmp=tmp[grep(unlist(strsplit(j,"\\."))[2],tmp$type),]
    t2=rbind(t2,tmp)
  }
  t2=t2[-1,]
  HRpredlist5[[i]]=t2
}
pdf("selected.forest.pdf")
a=HRpredlist5$good
forest(
  a[,c(1, 15)],  # 选择要在森林图中使用的数据列
  est = list(a$TCGA.HRpred,a$metabric.HRpred,a$SCANB.HRpred),
  lower = list(a$TCGA.CI_lower,a$metabric.CI_lower,a$SCANB.CI_lower),
  upper = list(a$TCGA.CI_upper,a$metabric.CI_upper,a$SCANB.CI_upper),
  ci_column = c(2),         # 指定CI列
  ref_line = 1,                # 添加参考线
  vert_line = c(0.5, 2),       # 添加垂直线
  nudge_y = 0.2,xlim = c(0,1)#,               # 垂直调整标签位置
  # theme = tm
)
a=HRpredlist5$poor
tmp=a$TCGA.CI_upper
tmp=tmp[!is.na(tmp)]
tmp=tmp[!is.infinite(tmp)]
a$TCGA.CI_upper[is.infinite(a$TCGA.CI_upper)]=max(tmp)
forest(
  a[,c(1, 15)],  # 选择要在森林图中使用的数据列
  est = list(a$TCGA.HRpred,a$metabric.HRpred,a$SCANB.HRpred),
  lower = list(a$TCGA.CI_lower,a$metabric.CI_lower,a$SCANB.CI_lower),
  upper = list(a$TCGA.CI_upper,a$metabric.CI_upper,a$SCANB.CI_upper),
  ci_column = c(2),         # 指定CI列
  ref_line = 1,                # 添加参考线
  vert_line = c(0.5, 2),       # 添加垂直线
  nudge_y = 0.2,xlim = c(0,10)#,               # 垂直调整标签位置
  # theme = tm
)
dev.off()






##compare with L1 Lasso
BRCA_list$TCGA$exp=log(BRCA_list$TCGA$exp+1)

# 加载必要的包
library(glmnet)      # Lasso-Cox 核心包
library(survival)    # 生存分析基础包
library(caret)       # 数据分割
library(ggplot2)     # 可视化

library(doParallel)

# 2. 设置并行计算核心数
# 在集群上，通常申请了多少个CPU核心就设置多少个
no_cores <- detectCores()  # 或手动指定，如 20
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# 3. 确认并行后端已注册
getDoParWorkers()  # 应显示核心数


cv_fitlist=list()
for(database in names(BRCA_list)[1:3]){
  print(database)
  expression_matrix=BRCA_list[[database]]$exp
  survival_data=BRCA_list[[database]]$clinical
  survival_data=survival_data[!is.na(survival_data$OS.time),]
  survival_data=survival_data[survival_data$OS.time>0,]
  expression_matrix=expression_matrix[,survival_data$ID]
  if(database=="TCGA"){
    survival_data=survival_data[survival_data$normalORnot%in%c("primary","metastasis"),]
    expression_matrix=expression_matrix[,survival_data$ID]
  }

  # 示例数据准备
  set.seed(123)
  
  # 1. 对齐样本顺序
  common_samples <- intersect(colnames(expression_matrix), 
                              survival_data$ID)
  expr_filtered <- expression_matrix[, common_samples]
  surv_filtered <- survival_data[match(common_samples, survival_data$ID), ]
  
  # 2. 转置表达矩阵（glmnet要求样本在行，基因在列）
  X <- t(expr_filtered)  # 现在维度: 样本 × 基因
  X=X[,apply(X,2,sum)>0]
  
  # 3. 创建生存对象
  y <- Surv(time = surv_filtered$OS.time, 
            event = surv_filtered$OS)
  
  # 4. 数据标准化（Lasso需要特征在同一尺度）
  X_scaled <- scale(X)  # 均值为0，标准差为1
  
  # 分割数据用于内部验证
  set.seed(123)
  train_index <- createDataPartition(surv_filtered$OS, 
                                     p = 0.7, 
                                     list = FALSE)
  
  X_train <- X_scaled[train_index, ]
  X_test <- X_scaled[-train_index, ]
  y_train <- y[train_index, ]
  y_test <- y[-train_index, ]
  
  
  # 执行带交叉验证的Lasso-Cox
  set.seed(123)
  cv_fit <- cv.glmnet(
    x = X_train,
    y = y_train,
    family = "cox",           # Cox比例风险模型
    alpha = 1,                 # alpha=1 表示Lasso L1惩罚
    nfolds = 10,               # 10折交叉验证
    type.measure = "deviance", # 评估指标：部分似然偏差
    parallel = TRUE,      # 启用并行
    standardize = FALSE        # 数据已经标准化过
  )
  cv_fitlist[[database]]=cv_fit
  # 查看交叉验证结果
  print(cv_fit)
  
  # 绘制交叉验证曲线
  pdf(paste(database,"cv.pdf"))
  plot(cv_fit)
  dev.off()
}
save(cv_fitlist,file = "cv_fitlist.RData")


for(database in names(cv_fitlist)){
  print(database)
  expression_matrix=BRCA_list[[database]]$exp
  survival_data=BRCA_list[[database]]$clinical
  survival_data=survival_data[!is.na(survival_data$OS.time),]
  survival_data=survival_data[survival_data$OS.time>0,]
  expression_matrix=expression_matrix[,survival_data$ID]
  if(database=="TCGA"){
    survival_data=survival_data[survival_data$normalORnot%in%c("primary","metastasis"),]
    expression_matrix=expression_matrix[,survival_data$ID]
  }
  
  # 示例数据准备
  set.seed(123)
  
  # 1. 对齐样本顺序
  common_samples <- intersect(colnames(expression_matrix), 
                              survival_data$ID)
  expr_filtered <- expression_matrix[, common_samples]
  surv_filtered <- survival_data[match(common_samples, survival_data$ID), ]
  
  # 2. 转置表达矩阵（glmnet要求样本在行，基因在列）
  X <- t(expr_filtered)  # 现在维度: 样本 × 基因
  X=X[,apply(X,2,sum)>0]
  
  # 3. 创建生存对象
  y <- Surv(time = surv_filtered$OS.time, 
            event = surv_filtered$OS)
  
  # 4. 数据标准化（Lasso需要特征在同一尺度）
  X_scaled <- scale(X)  # 均值为0，标准差为1
  
  # 分割数据用于内部验证
  set.seed(123)
  train_index <- createDataPartition(surv_filtered$OS, 
                                     p = 0.7, 
                                     list = FALSE)
  
  X_train <- X_scaled[train_index, ]
  X_test <- X_scaled[-train_index, ]
  y_train <- y[train_index, ]
  y_test <- y[-train_index, ]
  
  
  cv_fit=cv_fitlist[[database]]
  
  # 两种常用的λ选择策略
  # 1. lambda.min: 使偏差最小的λ（可能过拟合）
  lambda_min <- cv_fit$lambda.min
  cat("lambda.min:", lambda_min, "\n")
  
  # 2. lambda.1se: 最小偏差1个标准误内的最大λ（更稀疏，推荐）
  lambda_1se <- cv_fit$lambda.1se
  cat("lambda.1se:", lambda_1se, "\n")
  
  # 查看两种λ对应的非零系数数量
  coef_min <- coef(cv_fit, s = "lambda.min")
  coef_1se <- coef(cv_fit, s = "lambda.1se")
  
  n_genes_min <- sum(coef_min != 0)
  n_genes_1se <- sum(coef_1se != 0)
  
  cat("lambda.min 选择的基因数:", n_genes_min, "\n")
  cat("lambda.1se 选择的基因数:", n_genes_1se, "\n")
  
  # 使用lambda.1se提取非零系数
  selected_genes <- as.matrix(coef_1se)
  selected_genes <- selected_genes[selected_genes[, 1] != 0, , drop = FALSE]
  
  # 转换为数据框
  prognostic_genes <- data.frame(
    gene = rownames(selected_genes),
    coefficient = selected_genes[, 1],
    stringsAsFactors = FALSE
  )
  
  # 按系数绝对值排序
  prognostic_genes$abs_coef <- abs(prognostic_genes$coefficient)
  prognostic_genes <- prognostic_genes[order(prognostic_genes$abs_coef, 
                                             decreasing = TRUE), ]
  
  # 显示top基因
  print(head(prognostic_genes, 10))
  
  # 保存结果
  # write.csv(prognostic_genes, "lasso_selected_genes.csv", row.names = FALSE)
  
  
  # 计算每个样本的风险评分
  # 风险评分 = Σ(基因表达量 × Lasso系数)
  
  # 获取选中的基因列表
  genes_selected <- prognostic_genes$gene
  
  # 在原始数据中提取这些基因
  X_selected <- X_scaled[, genes_selected, drop = FALSE]
  
  # 计算风险评分
  risk_score <- X_selected %*% prognostic_genes$coefficient
  risk_score <- as.numeric(risk_score)
  names(risk_score) <- rownames(X_scaled)
  
  # 添加到生存数据
  surv_filtered$risk_score <- risk_score[surv_filtered$ID]
  
  # 根据中位数或最佳截点分组
  surv_filtered$risk_group <- ifelse(surv_filtered$risk_score > median(surv_filtered$risk_score, na.rm = TRUE),"High Risk", "Low Risk")
  
  # 1. Kaplan-Meier生存曲线
  library(survminer)
  
  fit_surv <- survfit(Surv(OS.time, OS) ~ risk_group, data = surv_filtered)
  
  ggsurvplot(
    fit_surv,
    data = surv_filtered,
    pval = TRUE,
    risk.table = TRUE,
    title = "Lasso Risk Score Stratification",
    xlab = "Time (months)",
    ylab = "Overall Survival Probability"
  )
  
  # 2. 时间依赖ROC曲线（评估预测准确性）
  library(timeROC)  ##R version 4.5
  
  roc_result <- timeROC(
    T = surv_filtered$OS.time,
    delta = surv_filtered$OS,
    marker = surv_filtered$risk_score,
    cause = 1,
    times = c(12, 36, 60),  # 1年，3年，5年
    iid = TRUE
  )
  
  # 查看AUC
  print(roc_result$AUC)
  
  # 3. 多因素Cox回归验证独立性
  cox_multi <- coxph(Surv(OS.time, OS) ~ risk_score + age + stage + grade, 
                     data = surv_filtered)
  summary(cox_multi)
}

# [1] "TCGA"
# lambda.min: 0.06017843 
# lambda.1se: 0.06017843 
# lambda.min 选择的基因数: 0 
# lambda.1se 选择的基因数: 0 
# [1] "metabric"
# lambda.min: 0.04232793 
# lambda.1se: 0.1418664 
# lambda.min 选择的基因数: 99 
# lambda.1se 选择的基因数: 0 
# [1] "SCANB"
# lambda.min: 0.01023283 
# lambda.1se: 0.0394325 
# lambda.min 选择的基因数: 376 
# lambda.1se 选择的基因数: 9