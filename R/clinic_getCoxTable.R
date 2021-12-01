
#' @title Get a model result of univariate and multivariate Cox regression Models
#' @description Get a model result of univariate and multivariate Cox regression Models
#' @param data a data frame containing cluster,time and status value
#' @param time.col colnames of time value
#' @param status.col  colnames of status value
#' @param is.multiple whether to do multiple cox regression
#' @param stepwise Whether use stepwise strategy to simple multiple cox model
#' @param direction method of stepwise.See \code{\link[MASS]{stepAIC}}.
#' @param cluster the cluster values like age and sex
#' @param control set the control value in a factor vector.For example, in T status,T1 should be a control in the contrast of T3 vs T1.
#' @param dig the decimal place of output numeric value.
#' @param names  part of file name.
#' @importFrom survival coxph cox.zph
#' @importFrom MASS stepAIC
#' @importFrom dplyr left_join
#' @seealso \code{\link[MASS]{stepAIC}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @export
getCoxTable <- function(data,
                        time.col,
                        status.col,
                        is.multiple = T,
                        stepwise = F,
                        direction="both",
                        cluster,
                        control,
                        dig=2,
                        names = "test1"){
  ## 包
  # nd <- c("survival","MASS","dplyr","readxl") # "broom",
  # Plus.library(nd)

  ## 提取数据
  data1 <- data[,c(time.col,status.col,cluster)]
  colnames(data1)[1:2] <- c("time","status")

  ## 重因子化
  for(i in 1:length(cluster)){
    control.i <- control[i]
    u.i <- unique(as.character(data1[,cluster[i]]))
    if(is.na(control.i)){
      #连续型变量，不需要因子化
      data1[,cluster[i]] <- data1[,cluster[i]]
    } else {
      #重因子化
      data1[,cluster[i]] <- factor(data1[,cluster[i]],levels = c(control.i,setdiff(u.i,control.i)))
    }

  }

  ### 单因素分析
  if(T){

    ## cox
    Univariate <- NULL;Univariate.ph <- NULL
    for(i in 1:length(cluster)){ # i=11
      c.i <- cluster[i]
      data1.i <- subset(data1,select = c("time","status",c.i))
      fit.i <- coxph(Surv(time,status)~., data = data1.i)
      result.i <- summary(fit.i)[["coefficients"]]
      Univariate <- rbind(Univariate,result.i)
      # ph假定
      if(!is.na(result.i[,'coef'])){
        ph.i <- cox.zph(fit.i)
        Univariate.ph <- rbind(Univariate.ph,ph.i$table)
      } else {
        # NA coef. Always due to lot of NA value in the data
        Univariate.ph <- rbind(Univariate.ph,
                               data.frame(row.names = rownames(result.i),
                                          chisq=NA,
                                          df=NA,
                                          p=NA))
      }
    }
    Univariate <- as.data.frame(Univariate)
    Univariate.ph <- as.data.frame(Univariate.ph)

    # result
    ci <- -round(qnorm((1-0.95)/2),2)
    Factor = rownames(Univariate)
    hr= Univariate$`exp(coef)`
    low.hr =  exp(Univariate$coef - ci*Univariate$`se(coef)`)
    up.hr =  exp(Univariate$coef + ci*Univariate$`se(coef)`)
    hr1 = paste(round(hr,dig),"(",round(low.hr,dig)," to ",round(up.hr,dig),")",sep = "")
    p.val = as.character(round2(Univariate$`Pr(>|z|)`,dig))
    df.u <- data.frame(
      Factor=Factor,
      hr=hr1,
      p=p.val
    )

  }

  ### 多因素分析
  if(is.multiple){
    ## 是否进行stepwise
    if(stepwise == T){
      ## 除去数据中的空值
      LuckyVerbose("运行MASS::stepAIC程序...")
      logic1 <- apply(data1,1,is.one.na)
      data2 <- data1[!logic1,];test1 <- table(logic1)[names(logic1) %in% T]
      if(T %in% logic1){LuckyVerbose(paste0("There are ",test1," rows containing at least 1 NA value."),levels = 2)}
      ## cox回归
      fit <- coxph(Surv(time,status) ~.,data = data2)
      ##stepwise
      step1 <- stepAIC(fit,direction=direction)
      fit2 <- coxph(step1[["formula"]],data = data2)
      LuckyVerbose("完成逐步AIC筛选法!")
    } else {
      fit2 <- coxph(Surv(time,status) ~.,data = data1)
    }

    ## multiple
    m <- summary(fit2)
    multivariate <- m[["coefficients"]]
    multivariate <- as.data.frame(multivariate)
    # multiple ph
    a <- base::tryCatch(multivariate.ph <- cox.zph(fit2), error = function(e)e, finally = NULL)
    if(is.null(a$message)){
      multivariate.ph <- as.data.frame(multivariate.ph$table)
    } else {
      multivariate.ph <- a$message
    }

    ## result
    Factor = rownames(multivariate)
    hr= multivariate$`exp(coef)`
    low.hr =  exp(multivariate$coef - ci*multivariate$`se(coef)`)
    up.hr =  exp(multivariate$coef + ci*multivariate$`se(coef)`)
    hr1 = paste(round(hr,dig),"(",round(low.hr,dig)," to ",round(up.hr,dig),")",sep = "")
    p.val = as.character(round2(multivariate$`Pr(>|z|)`,dig))
    df.m <- data.frame(
      Factor=Factor,
      hr=hr1,
      p=p.val
    )
  } else {
    df.m <- data.frame(Factor = NA)
    multivariate = NULL
    multivariate.ph = NULL
  }

  df1 <- left_join(df.u,df.m,by="Factor")
  df1 <- as.matrix(df1)

  ## 重命名行名
  rn <- df1[,"Factor"]
  newrownames <- function(rn,
                          cluster.i,
                          control.i){
    if(is.na(control.i)){
      #说明名字不需要更改
      name <-  cluster.i
    } else {
      #说明名字要更改
      p1 <- grep(cluster.i,rn)
      rn.i <- rn[p1]
      n2 <- Fastextra(rn.i,cluster.i,2)
      n1 <- rep(cluster.i,length(n2))
      name <- paste(n1,n2,sep = "_")
      name <- paste(name,control.i,sep = " vs. ")}
    return(as.character(name))
  }
  test1 <- data.frame(cluster = cluster,control = control)
  new.rn <- apply(test1,1,function(x)newrownames(rn,x[1],x[2]))
  df1[,"Factor"] <- unlist(new.rn)

  ## 进一步处理
  title <- c("Factor","HR(95%CI)","P","HR(95%CI)","P")[1:ncol(df1)]
  title1 <- c(" ","Univariate analysis"," ","Multivariable analysis"," ")[1:ncol(df1)]
  df2 <- rbind(title1,title,df1)
  colnames(df2) <- c("Factor","Univariate.HR(95%CI)","Univariate.P","Multivariable.HR(95%CI)","Multivariable.P")[1:ncol(df2)]

  ##输出结果
  write.csv(df2,paste0(names,"_Cox proportional hazards models.csv"),row.names = F)
  LuckyVerbose("Cox模型结果已保存在当前工作目录。")
  l <- list(
    Univariate=Univariate,
    Multivariate = multivariate,
    PHtest = list(Univariate = Univariate.ph,
                  Multivariate =  multivariate.ph)
  )
  return(l)
}





