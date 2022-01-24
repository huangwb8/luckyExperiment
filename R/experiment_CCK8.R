


# Class "LuckyCCK8"
#' @name luckyExperiment-classes
#' @title LuckyCCK8-class
#' @docType class
#' @slot Repeat repeat parameters
#' @slot Data datasets
#' @slot Plot ggplot
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @keywords classes
#' @exportClass LuckyCCK8
setClass(
  Class = "LuckyCCK8",
  slots = c(
    Repeat = "list",
    Data = "list",
    Plot = "list"
  )
)


#' @title Fast methods for CCK8 data analysis
#' @description Fast methods for CCK8 data analysis
#' @param data a data frame
#' @param cell.type the colname of cell type
#' @param group the colname of group
#' @param empty.group the colname of the reference group without cells in the 96-well board.
#' @param stimulate.time the colname of a stimulating factor.For example,different dose of drugs;different treatments.
#' @param observe.time colname.observation time after CCK8 treatment.Time nodes like 1h,2h,4h are recommanded.
#' @param value colname of OD value
#' @param x.lab the title of x axis
#' @param y.lab the title of y axis
#' @param title the title of plot
#' @param size the size of plot
#' @param verbose whether do report
#' @importFrom Rmisc summarySE
#' @importFrom plyr ddply summarise
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot facet_wrap geom_line geom_point geom_errorbar labs theme_bw theme
#' @importFrom ggpubr compare_means
#' @return a \code{LuckyCCK8} object
#' @seealso \code{\link{plot}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' data <- readxl::read_xlsx("E:/iProjects/exosome/experiment/exosome_1A/CCK8/CCK8_HGC27的合适种植数量及观察时间_20190424.xlsx", na = "NA")
#'
#' object <- lucky::FastCCK8(data)
#' View(object$Data$statistics$two)
#' View(object$Data$statistics$all)
#'
#' select <- summary(object,observe_time = "2h")
#' @export
FastCCK8 <- function(data,
                     cell.type = "cell",
                     group = "group",
                     empty.group = "0",
                     stimulate.time = "stimulate.time",
                     observe.time = "observe.time",
                     value = "od",
                     x.lab = "time",
                     y.lab = "OD value(450nm)",
                     title = "",
                     size = 20,
                     method.two = "t.test",
                     method.all = "anova",
                     verbose = T){
  ## 加载包
  # nd <- c("ggplot2","plyr","ggpubr");Plus.library(nd)

  ## 数据预处理
  data$group <- factor(data$group,levels = as.character(unique(data$group)))
  select <- c(cell.type,group,stimulate.time,observe.time,value)
  data2 <- data[select]
  colnames(data2) <- c("cell","group","stimulate.time","observe.time","value")
  data3 <- Rmisc::summarySE(data2, measurevar="value", groupvars=c("cell","group","stimulate.time","observe.time"))

  ## ggplot for quality control
  p_qc <- ggplot(data3,aes(x = stimulate.time,y = value,color = group,group=group)) +
    facet_wrap(.~ observe.time) +
    geom_line(size = 1.5) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin=value-se, ymax=value+se),width=0.1,size=1.5) +
    labs(title = "Quality Control") +
    theme_bw() +
    theme(
      plot.title = element_text(size = (size/20)*18,face = "bold",hjust = 0.5),
      axis.title = element_text(size = (size/20)*15,face = "bold"),
      axis.text = element_text(size = (size/20)*12,face = "bold"),
      legend.title=element_blank(),
      legend.position = "bottom",
      legend.text =element_text(size = (size/20)*12,face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = (size/20)*12,face = "bold"))
  if(verbose) win.graph(10,10);print(p_qc)

  ## standard data via 0 group
  data_real <- data2[!data2$group %in% empty.group,]
  data_reference <- data2[data2$group %in% empty.group,]
  data_adj <- exp_autoref(data_real,data_reference)

  ## standard data plot
  data_adj_2 <- Rmisc::summarySE(data_adj, measurevar="adj.value", groupvars=c("cell","group","stimulate.time","observe.time")) #   str(data_adj_2) # ?Rmisc::summarySE
  data_adj_2$group <- factor(data_adj_2$group,levels = as.character(unique(data_adj_2$group)))
  p_sd <- ggplot(data_adj_2,aes(x = stimulate.time,y = adj.value,color = group,group=group)) +
    facet_wrap(.~ observe.time) +
    geom_line(size = 1.5) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin=adj.value-se, ymax=adj.value+se),width=0.1,size=1.5) +
    labs(x = x.lab,y = y.lab,title = title) +
    theme_bw() +
    theme(
      plot.title = element_text(size = (size/20)*18,face = "bold",hjust = 0.5),
      axis.title = element_text(size = (size/20)*15,face = "bold"),
      axis.text = element_text(size = (size/20)*12,face = "bold"),
      legend.title=element_blank(),
      legend.position = "bottom",
      legend.text =element_text(size = (size/20)*12,face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = (size/20)*12,face = "bold"));
  if(verbose) win.graph(10,10);print(p_sd)

  ## statistics
  res_two_test <- ddply(.data = data_adj, .variables = c("observe.time","stimulate.time"), .fun = function(x)two_test(x,method.two))
  res_two_test <- res_two_test[-match(".y.",colnames(res_two_test))]
  res_all_test <- ddply(.data = data_adj, .variables = c("observe.time","stimulate.time"), .fun = function(x)all_test(x,method.all))
  res_all_test <- res_all_test[-match(".y.",colnames(res_all_test))]

  ## Output Data
  l <- list(
    Repeat = list(
      cell.type = cell.type,
      group = group,
      empty.group = empty.group,
      stimulate.time = stimulate.time,
      observe.time = observe.time,
      value = value
    ),
    Data = list(
      metaData = data,
      adjData = data_adj,
      plotData = data_adj_2,
      statistics = list(two = res_two_test,
                        all = res_all_test)
    ),
    Plot = list(
      qcPlot = p_qc,
      realPlot = p_sd
    )
  )
  class(l) <- "LuckyCCK8"
  return(l)
}


#### S3
#' @name summary
#' @aliases summary
#' @docType methods
#' @rdname summary-methods
#' @title summary method
#' @description plot method for LuckyCCK8 Object
#' @return New plot of CCK8 result
#' @inheritParams FastCCK8
#' @importFrom stats4 summary
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot facet_wrap geom_line geom_point geom_errorbar labs theme_bw theme
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @exportMethod summary
setMethod(
  "summary",
  signature(object="LuckyCCK8"),
  function (object,
            observe_time = "2h",
            x.lab = "time",
            y.lab = "OD value(450nm)",
            title = "",
            size = 20,
            legend.position = "bottom",
            verbose = T){
    ## select a special data for plotting
    ds <- object$Data$plotData
    ds_2 <- dplyr::filter(ds,observe.time %in% observe_time)
    ## new plot
    p <- ggplot(ds_2,aes(x = stimulate.time,y = adj.value,color = group,group=group)) +
      geom_line(size = 1.5) +
      geom_point(size = 2) +
      geom_errorbar(aes(ymin=adj.value-se, ymax=adj.value+se),width=0.1,size=1.5) +
      labs(x = x.lab,y = y.lab,title = title) +
      theme_bw() +
      theme(
        plot.title = element_text(size = (size/20)*18,face = "bold",hjust = 0.5),
        axis.title = element_text(size = (size/20)*15,face = "bold"),
        axis.text = element_text(size = (size/20)*12,face = "bold"),
        legend.title=element_blank(),
        legend.position = legend.position,
        legend.text =element_text(size = (size/20)*12,face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = (size/20)*12,face = "bold"));
    if(verbose) print(p)
    return(p)
  }
)



####=============Assistant Functions==============####
## 根据给定的参考数据集进行标准化
#' @importFrom plyr ddply summarise
#' @importFrom dplyr filter
exp_autoref <- function(data_real,
                        data_reference){

  ## mean of reference data
  g <- setdiff(colnames(data_reference),"value")
  df.r <- ddply(.data = data_reference,
                .variables = g,
                plyr::summarise,
                mean = mean(value))

  ## relative value of real data
  get1 <- function(vt,df.r){
    # vt = data_real[1,]
    ns <- names(vt)
    vt2 <- as.character(as.matrix(vt))
    ref.i <- dplyr::filter(df.r,cell %in% vt2[1],stimulate.time %in% vt2[3],observe.time %in% vt2[4])
    ref.i <- as.numeric(ref.i[,"mean"])
    adjust.value <- as.numeric(vt[ns %in% "value"]) - ref.i
    return(adjust.value)
  }
  data_real$adj.value <- apply(data_real,1,
                               function(x)get1(x,df.r))
  return(data_real)
}

## 计算不同组的两两比较（t检验）
two_test <- function(data_adj_i,method.two){
  # data_adj_i = data_adj # head(data_adj_i)
  stat_res <- compare_means(adj.value~group,data=data_adj_i,method = method.two,paired = F)
  stat_res <- as.data.frame(stat_res,stringsAsFactors = F)
  return(stat_res)
}

all_test <- function(data_adj_i,method.all){
  # data_adj_i = data_adj_2 # head(data_adj_i)
  stat_res <- compare_means(adj.value~group,data=data_adj_i,method = method.all)
  stat_res <- as.data.frame(stat_res,stringsAsFactors = F)
  return(stat_res)
}









