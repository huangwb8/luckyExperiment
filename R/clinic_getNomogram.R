

#' @title Get nomogram,standard curves and return related data and formula
#' @description Get nomogram,standard curves and return related data and formula
#' @param data a data frame
#' @param time.col survival data colnames_time
#' @param status.col survival data colnames_status
#' @param cluster The vector of the marker values for which we want to compute the time-dependent ROC curves. Without loss of generality, the function assumes that larger values of the marker are associated with higher risks of events. If lower values of the marker are associated with higher risks of events, then reverse the association adding a minus to the marker values.
#' @param time.raw.type  one of "Day","Month","Year"
#' @param time.target.type  one of "Day","Month","Year"
#' @param time.knot time knot to create survival probability
#' @param dig Digits in nomogram plot.Default is 9.
#' @param nomo.width width of nomogram plot
#' @param nomo.height height of nomogram plot
#' @param nomo.xfrac saved size of nomogram plot
#' @param standard.cmethod cmethod of \code{\link[rms]{calibrate}}.One of 'hare' and 'KM'
#' @param standard.method method of \code{\link[rms]{calibrate}}.Default is 'boot'.
#' @param standard.nspot the number of predictive spots in standard curves.
#' @param standard.names type of prediction.Like "OS".
#' @param standard.ylim Default is c(0,1)
#' @param standard.width parameters of output plot
#' @param standard.height parameters of output plot
#' @param names # part of file names.
#' @seealso \code{\link[rms]{calibrate}};\code{\link[rms]{cph}};\code{\link[rms]{nomogram}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' library(lucky)
#' library(survival)
#' data(lung)
#' data=lung
#' time.col = "time"
#' status.col = "status"
#' colnames(data)
#' cluster = c("sex","age","ph.ecog","ph.karno") #
#' time.raw.type=c("Day","Month","Year")[1]
#' time.target.type=c("Day","Month","Year")[2]
#' time.knot = c(12,24)
#' dig=9
#' standard.nspot = c(4,3)
#' width = 15
#' height = 12
#' xfrac=0.45
#' names = "love1"
#' standard.names = "OS"
#' standard.ylim = c(0,1)
#' standard.width = 10
#' standard.height = 10
#'
#' ## draw a nomogram and output related data and formula
#' warnings(list1 <- getNomogram(data,
#'                               time.col,
#'                               status.col,
#'                               cluster,
#'                               time.raw.type=c("Day","Month","Year")[1],
#'                               time.target.type=c("Day","Month","Year")[2],
#'                               standard.cmethod = "KM",
#'                               time.knot = c(12,24),
#'                               standard.nspot = c(4,3),
#'                               standard.names = "OS",
#'                               names = "love1"))
#'
#' #more complex example
#' list1 <- getNomogram(data,
#'                      time.col,status.col,
#'                      cluster,
#'                      time.raw.type=c("Day","Month","Year")[1],
#'                      time.target.type=c("Day","Month","Year")[2],
#'                      time.knot = c(12,24),
#'                      dig=9,
#'                      nomo.width = 15,
#'                      nomo.height = 12,
#'                      nomo.xfrac=0.45,
#'                      standard.nspot = 3,
#'                      standard.names = "OS",
#'                      standard.ylim = c(0,1),
#'                      standard.width = 10,
#'                      standard.height = 10,
#'                      names = "love1")
#'
#' ## calculate a time-specified survival probability
#' S <- list1[["formula.nomo"]]
#' S(25) #25 months specified survival prediction.
#'
#' ## other information
#' View(list1[["data.nomo"]])
#' View(list1[["data.standard"]])
#' @export
getNomogram <- function(data,
                        time.col,status.col,
                        cluster,
                        time.raw.type=c("Day","Month","Year")[1],
                        time.target.type=c("Day","Month","Year")[2],
                        time.knot = c(12,24),
                        dig=9,
                        nomo.width = 15,
                        nomo.height = 12,
                        nomo.xfrac=0.45,
                        standard.cmethod = c('hare','KM')[2],
                        standard.method = "boot",
                        standard.nspot = 3,
                        standard.names = "OS",
                        standard.ylim = c(0,1),
                        standard.width = 10,
                        standard.height = 10,
                        names = "love1"){
  ## 包
  need <- c("rms","survival","nomogramEx","polspline")
  Plus.library(need)

  ## save path
  ### grobal option
  dir_output <- "./Nomogram"
  dir.create(dir_output,recursive = T,showWarnings = F)

  ## 选择data
  data1 <- as.data.frame(data[,c(time.col,status.col,cluster)],stringsAsFactors = F)
  colnames(data1)[1:2] <- c("time","status")

  ## 去除含有NA值的行
  test.NA <- apply(data1,1,is.one.na)
  data1 <- data1[!test.NA,]
  if(T %in% test.NA){
    x1 <- grep(T,test.NA)
    report.na <- paste0("There are some rows with NA value:",paste(x1,collapse = "_"),";they had been removed.")
    LuckyVerbose(report.na)
  }

  ## 统一time的单位
  # getTimeROC::convert.time
  data1$time <- convert.time(data1$time,
                             from.ts = time.raw.type,
                             to.ts = time.target.type)
  #report
  LuckyVerbose(paste0("The max of time is ",round(max(data1$time),2)," ",time.target.type,"s,and the min is ",round(min(data1$time),2)," ",time.target.type,"s.The time.knot parameter must be in this range."))

  ## 测试time.knot
  logic1 <- time.knot > max(data1$time)
  if(T %in% logic1){
    ## time.knot参数有误
    stop("Error!The time.knot is out of the time range!")
  } else {
    ## time.knot参数正确
    LuckyVerbose("Time.knot is in the right range.")
    #logic2 <- apply(data1,1,is.one.na);data1 <- data1[!logic2,]
    #units(data1$time) <- time.target.type #时间单位
    dd <<- datadist(data1);#dd对象一定要设置为全局变量才可以。
    options(datadist="dd")

    ## get cox fomula
    getCoxf <- function(data.i,
                        time = "time",
                        event = "status"){
      c.i <- setdiff(colnames(data.i),c(time,event))
      f1 <- paste0('`',c.i,"`") %>% paste(.,collapse = " + ") %>% paste0("Surv(",time,",",event,") ~ ",.)
      f2 <- as.formula(f1)
      return(f2)
    }
    f <- cph(getCoxf(data1),data = data1,x=TRUE,y=TRUE,surv=TRUE)
    surv <- Survival(f)

    ## 构建fun.list
    fun.list <- NULL
    for(i in time.knot){
      p <- match(i,time.knot)
      a <- alist(x = , i = 1, surv(i,x));a$i <- i
      fun.i <- as.function(a)
      fun.list <- c(fun.list,list(fun.i))
      names(fun.list)[p] <- time.knot[p]
    }

    ## 构建fun.label
    label1 <- paste0(time.target.type," Survival Probability")
    funlabel <- paste(time.knot,label1,sep = "-")
    names(fun.list) <- funlabel

    ## 绘制nomo图
    nomo <- nomogram(f,fun=fun.list)
    plot(nomo,xfrac=.45)
    pdf(paste0(dir_output,"/",names,"_nomogram.pdf"),width = nomo.width,height = nomo.height)
    plot(nomo,xfrac=nomo.xfrac)
    dev.off()
    LuckyVerbose("The nomograph had been saved in current workspace!")

    ## 绘制标准曲线
    if(length(standard.nspot)==1){
      standard.nspot <- rep(standard.nspot,length(time.knot))
    } else {
      standard.nspot <- standard.nspot
    }
    list.cal <- NULL # str(data1)
    for(i in 1:length(time.knot)){#i=1
      f1 <- cph(getCoxf(data1),data = data1,x=TRUE,y=TRUE,surv=TRUE,time.inc = time.knot[i])
      cal <- calibrate(f1,
                       cmethod=standard.cmethod,
                       method=standard.method,
                       u=time.knot[i],
                       m=floor(nrow(data1)/standard.nspot[i]),
                       B = 1000)
      list.cal <- c(list.cal,list(cal))
      n.i <- paste0(time.knot[i],"-",time.target.type)
      names(list.cal)[i] <- n.i

      #设置x轴的范围
      range <- cal[,"mean.predicted"]
      x.min <- ifelse(round(min(range),1) <= min(range),round(min(range),1),round(min(range),1)-0.1)
      x.max <- ifelse(round(max(range),1) <= max(range),round(max(range),1)+0.1,round(max(range),1))
      #设置x/y轴的标题
      xlab=paste0("Nomogram-Predicted Probability of ",time.knot[i],"-",time.target.type," ",standard.names)
      ylab=paste0("Actual ",time.knot[i],"-",time.target.type," ",standard.names,"(proportion)")

      #plot styple
      pdf(paste0(dir_output,"/",names,"_",n.i,"_predicted vs. observed standard plot.pdf"),width = standard.width,height = standard.height)
      par(cex=1.2) #全局设置
      plot(cal,lwd=2,lty=1,
           errbar.col=mycolor[4],
           xlim=c(x.min,x.max),
           ylim=standard.ylim,
           xlab=xlab,ylab=ylab,
           col=mycolor[1],
           subtitles = F)
      lines(cal[,c("mean.predicted",standard.cmethod)],type="b",lwd=2,col=mycolor[1], pch=16)
      abline(0,1,lty=2,lwd=2,col="black")
      par(las = 0) #默认设置
      dev.off()

    }
    LuckyVerbose("Complete the standard validation plots!")

    ## 公式
    f.list <- nomogramEx(nomo=nomo,np=length(time.knot),digit=dig)
    names(f.list)[2:length(f.list)] <- c(cluster,funlabel)

    ## 预测值
    data2 <- data1
    for(i in 1:length(funlabel)){
      data2$x <- surv(time.knot[i],f$linear.predictors)
      colnames(data2)[ncol(data2)] <- funlabel[i]
    }

    ## 输出结果
    result <- list(
      nomo = nomo,
      data.nomo = data2,
      data.standard = list.cal,
      formula.nomo = function(x)surv(x,f$linear.predictors),
      index.nomo = f.list)
    saveRDS(result,paste0(dir_output,"/",names,"_resultNomogram.rds"))
    LuckyVerbose("All done!")
    return(result)
  }
}

####====================Assistant Function===================####
# getCoxfomula











