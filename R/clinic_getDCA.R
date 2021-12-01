

#' @title Fast way to draw net benefit curves
#' @description Fast way to draw net benefit curves
#' @param data a data frame with survival data(time & event)
#' @param group.by group colname
#' @param time.col colname of time value
#' @param status.col Numeric.\code{status.col} must be a binary parameter with 0(event not happened yet) or 1(event happened) value.
#' @param clusters list of characters.Some of the colnames of the data representing the value involving the glm model.
#' @param plot.confidence.intervals whether plot confidence intervals
#' @param plot.type Character. One/some of c("time","cluster","disperse"). if \code{plot.type="time"},then multiple clusters would be ploted for every time; if \code{plot.type="cluster"}, then multiple times would be ploted for every cluster; if \code{plot.type="disperse"}, then every single plot would be ploted indivitually.
#' @param plot.title.position List. The title position of the plot.Only available when \code{plot.type="time"} or \code{plot.type="cluster"}. Its length must be the same of \code{time.knot}
#' @param plot.label List. The name of list are the paste0 of elements in one of \code{clusters} with "-".The content of the list are the real names that would be printed in the plot legend/title.
#' @param names part of save files or plots.
#' @inheritParams getNomogram
#' @inheritParams rmda::decision_curve
#' @inheritParams rmda::plot_decision_curve
#' @importFrom plyr adply
#' @importFrom rmda decision_curve plot_decision_curve
#' @return DCA plots and related model list
#' @details The length of \code{population.prevalence} should be equal to \code{time.knot}. Default of \code{population.prevalence} is \code{c(NA,NA)}, \cr which means we use an estimated population prevalence based on raw dataset.
#' @seealso \code{\link[rmda]{plot_decision_curve}}.
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## data preparation
#' data("myeloma",package = "survminer");
#' data <- myeloma;rm(myeloma)
#' data$group <- ifelse(data$molecular_group %in% c( "Cyclin D-1", "Cyclin D-2", "Proliferation"),"proliferation","others")
#'
#' ## Here is a useful example:
#'
#' # group by
#' a <- getDCA(data,
#'             group.by = "group",
#'             time.col = "time",
#'             status.col = "event",
#'             clusters = list(c("TP53","WHSC1"),
#'                             c("TP53","WHSC1","DEPDC1")),
#'             time.raw.type=c("Day","Month","Year")[2],
#'             time.target.type=c("Day","Month","Year")[3],
#'             time.knot = c(3,5),
#'             plot.type = c("time","cluster","disperse"),
#'             plot.title.position = list("3" = c(0.85,0.24),
#'                                        "5" = c(0.85,0.24)),
#'             plot.label = list("TP53-WHSC1" = "simpleModel",
#'                   "TP53-WHSC1-DEPDC1" = "complexModel"),
#'             names="project")
#'
#' # none group
#' a <- getDCA(data,
#'             group.by = NULL,
#'             time.col = "time",
#'             status.col = "event",
#'             clusters = list(c("TP53","WHSC1"),
#'                             c("TP53","WHSC1","DEPDC1")),
#'             time.raw.type=c("Day","Month","Year")[2],
#'             time.target.type=c("Day","Month","Year")[3],
#'             time.knot = c(3,5),
#'             plot.type = c("time","cluster","disperse"),
#'             plot.title.position = list("3" = c(0.85,0.24),
#'                                        "5" = c(0.85,0.24)),
#'             plot.label = list("TP53-WHSC1" = "simpleModel",
#'                   "TP53-WHSC1-DEPDC1" = "complexModel"),
#'             names="project")
#'
#' ## Of cause, you can select a sigle plot type like "time", which
#' ## is the most comman calling:
#' a <- getDCA(data,
#'             group.by = "group",
#'             time.col = "time",
#'             status.col = "event",
#'             clusters = list(c("TP53","WHSC1"),
#'                             c("TP53","WHSC1","DEPDC1")),
#'             time.raw.type=c("Day","Month","Year")[2],
#'             time.target.type=c("Day","Month","Year")[3],
#'             time.knot = c(3,5),
#'             plot.type = c("time","cluster","disperse")[1],
#'             plot.title.position = list("3" = c(0.85,0.24),
#'                                        "5" = c(0.85,0.24)),
#'             plot.label = list("TP53-WHSC1" = "simpleModel",
#'                   "TP53-WHSC1-DEPDC1" = "complexModel"),
#'             names="project")
#' @export
getDCA <- function(data,
                   group.by = NULL,
                   time.col = "time",
                   status.col = "event",
                   clusters,
                   time.raw.type=c("Day","Month","Year")[1],
                   time.target.type=c("Day","Month","Year")[3],
                   time.knot = c(3,5),
                   family = binomial(link ='logit'),
                   thresholds= seq(0,1, by = 0.01),
                   confidence.intervals = 0.95,
                   study.design = 'cohort',
                   population.prevalence = c(NA,NA),
                   cost.benefit.axis =FALSE,
                   col= grDevices::rainbow(n=8,v=0.8),
                   plot.confidence.intervals=FALSE,
                   standardize = FALSE,
                   plot.type = c("time","cluster","disperse")[1],
                   plot.title.position = list("3" = c(0.85,0.24),
                                              "5" = c(0.85,0.24)),
                   plot.label,
                   names="project"){

  ### grobal option
  old_par <- par()
  dir_output <- "./decisionCurveAnalysis"
  dir.create(dir_output,recursive = T,showWarnings = F)

  ### time convertion
  data.x <- as.data.frame(data,stringsAsFactors = F)
  data.x[,time.col] <- convert.time(data.x[,time.col],
                             from.ts = time.raw.type,
                             to.ts = time.target.type)
  # report
  LuckyVerbose(paste0("The max of time is ",round(max(data.x[,time.col]),2)," ",time.target.type,"s,and the min is ",round(min(data.x[,time.col]),2)," ",time.target.type,"s. Note: time.knot must be in the interval!"))

  ### dca data accession
  L <- list()
  for(n in 1:length(clusters)){ # n=1

    ## select one cluster
    cluster <- clusters[[n]]
    cluster.name <- paste0(cluster,collapse = "_")

    ## data filtering
    data1 <- data.x[,c(time.col,status.col,cluster,group.by)]
    colnames(data1)[1:2] <- c("time","status")

    ## delete na value
    test.NA <- apply(data1,1,is.one.na)
    data1 <- data1[!test.NA,]
    if(T %in% test.NA){
      x1 <- grep(T,test.NA)
      LuckyVerbose("There are some rows with NA value:",paste(x1,collapse = "_"),";they had been removed.")
    }

    ### get group name
    if(is.null(group.by)){
      group.name="NonGroup"
      group_vector=rep("NonGroup",nrow(data1))
    } else {
      group.name=unique(as.character(data1[,group.by]))
      group_vector=as.character(data1[,group.by])
      data1 <- data1[-match(group.by,colnames(data1))]
    }

    ## get outcome data
    for(gn in group.name){

      # data for every group
      data.gn <- data1[group_vector %in% gn,]

      # get dca data
      for(i in 1:length(time.knot)){ # i=1
        knot.i <- time.knot[i]
        outcome.i_name <- paste0("outcome_",knot.i,"_",time.target.type)

        ## rawData
        data.i <- adply(data.gn,1,
                        function(x)getOutcome(x,knot.i = knot.i),
                        .id =outcome.i_name )
        colnames(data.i)[ncol(data.i)] <- outcome.i_name
        L[[outcome.i_name]][[gn]][[cluster.name]][["Data"]] <- data.i


        ## DCA data
        LuckyVerbose("=======")
        LuckyVerbose(cluster.name,": DCA for ",outcome.i_name,"...")
        data.i_2 <- data.i[,-grep("time|status",colnames(data.i))]
        population.prevalence.i <- getPopPre(data.i_2,population.prevalence[i])
        if(study.design == "cohort"){
          LuckyVerbose("study.design = cohort")
          dca.i <- decision_curve(
            getGlmFormula(data.i_2,outcome.i_name),
            data = data.i_2,
            family = family,
            thresholds= thresholds,
            confidence.intervals = confidence.intervals,
            study.design = study.design
          )
        } else {
          LuckyVerbose("study.design = case-control")
          dca.i <- decision_curve(
            getGlmFormula(data.i_2,outcome.i_name),
            data = data.i_2,
            family = family,
            thresholds= thresholds,
            confidence.intervals = confidence.intervals,
            study.design = study.design,
            population.prevalence = population.prevalence.i
          )
        }
        L[[outcome.i_name]][[gn]][[cluster.name]][["DCA"]] <- dca.i
      }
    }

  }

  ### get DCA plot
  nameTime <- names(L)
  nameGroup <- names(L[[1]])
  nameCluster <- names(L[[1]][[1]])

  ## time
  if("time" %in% plot.type){
    LuckyVerbose("=======")
    LuckyVerbose("o(*￣▽￣*)o: multiple clusters would be ploted for every time...")

    ## get data
    plotData <- list()
    for(i in 1:length(nameTime)){ # i=1
      for(j in 1:length(nameGroup)){
        for(k in 1:length(nameCluster)){
          nGC <- paste0(nameGroup[j],'--',nameCluster[k])
          plotData[[nameTime[i]]][[nGC]]<- L[[nameTime[i]]][[nameGroup[j]]][[nameCluster[k]]][["DCA"]]
        }
      }
    }

    ## DCA plot
    par(cex=2,lwd=3)
    pdf(paste0(dir_output,"/DCAplot_bytime_",names,".pdf"),7,7)
    for(i in 1:length(nameTime)){ # i=1
      plotData_i <- plotData[[i]]

      ## curve name
      t.i<- Fastextra(names(plotData)[i],"_",2)
      title <- paste0(t.i,"-",time.target.type)
      curve.names <- getDCALabel2(names(plotData_i),plot.label)
      curve.names2 <- sort(curve.names)
      plotData_i2 <- plotData_i[Fastmatch(curve.names2,curve.names)]
      curve.names3 <- rmNonGroupLabel(curve.names2)

      ## label position
      tp <- plot.title.position[[t.i]]

      plot_decision_curve(plotData_i2,
                          curve.names=curve.names3, # 改动
                          cost.benefit.axis =cost.benefit.axis,
                          col= col,
                          confidence.intervals=plot.confidence.intervals,
                          standardize = standardize)
      text(tp[1],tp[2],title)
    }
    dev.off()
    par(cex=old_par$cex,lwd=old_par$lwd)
  }

  ## cluster
  if("cluster" %in% plot.type){
    LuckyVerbose("=======")
    LuckyVerbose("(*▔＾▔*) : multiple time would be ploted for every cluster...")

    ## get data
    plotData <- list()
    for(i in 1:length(nameTime)){ # i=1
      for(j in 1:length(nameGroup)){
        for(k in 1:length(nameCluster)){
          nTG <- paste0(nameGroup[j],"--",nameTime[i])
          plotData[[nameCluster[k]]][[nTG]]<- L[[nameTime[i]]][[nameGroup[j]]][[nameCluster[k]]][["DCA"]]
        }
      }
    }

    ## DCA plot
    par(cex=2,lwd=3)
    pdf(paste0(dir_output,"/DCAplot_bycluster_",names,".pdf"),7,7)
    for(i in 1:length(nameCluster)){ # i=1
      plotData_i <- plotData[[i]]
      nTime.i <- paste0(Fastextra(Fastextra(names(plotData_i),"--",2),"_",2)," ",time.target.type)
      nGroup.i <- Fastextra(names(plotData_i),"--",1)
      curve.names <- paste(nGroup.i,nTime.i,sep = ": ")
      curve.names2 <- sort(curve.names)
      plotData_i <- plotData_i[Fastmatch(curve.names2,curve.names)]
      curve.names3 <- rmNonGroupLabel(curve.names2)
      plot_decision_curve(plotData_i,
                          curve.names=curve.names3,
                          cost.benefit.axis =cost.benefit.axis,
                          col= col,
                          confidence.intervals=plot.confidence.intervals,
                          standardize = standardize)

      text(plot.title.position[[1]][1],
           plot.title.position[[1]][2],
           getDCALabel(names(plotData)[i],plot.label))
    }
    dev.off()
    par(cex=old_par$cex,lwd=old_par$lwd)
  }

  ## disperse
  if("disperse" %in% plot.type){
    LuckyVerbose("=======")
    LuckyVerbose("ㄟ(▔,▔)ㄏ : indivitual ploting...")

    ## get data
    plotData <- list() # i <- j <- k <- 1
    for(i in 1:length(nameTime)){ # i=1
      for(j in 1:length(nameGroup)){
        for(k in 1:length(nameCluster)){
          n <- paste(nameTime[i],nameGroup[j],nameCluster[k],sep  = "--")
          plotData[[n]]<- L[[nameTime[i]]][[nameGroup[j]]][[nameCluster[k]]][["DCA"]]
        }
      }
    }

    ## DCA plot
    par(cex=2,lwd=3)
    pdf(paste0(dir_output,"/DCAplot_bydisperse_",names,".pdf"),7,7)
    for(i in 1:length(plotData)){ # i=1
      plotData_i <- plotData[[i]]
      nCluster.i <- Fastextra(names(plotData)[i],"--",3)
      nC <- getDCALabel(nCluster.i,plot.label)

      nTime.i <- Fastextra(names(plotData)[i],"--",1)
      nT <- paste0(Fastextra(nTime.i,"_",2),"-",time.target.type)

      nG <- Fastextra(names(plotData)[i],"--",2)

      curve.names <- paste0(nG,": ",nT," ",nC)
      curve.names2 <- rmNonGroupLabel(curve.names)

      # DCA plot
      plot_decision_curve(plotData_i,
                          curve.names = curve.names2,
                          cost.benefit.axis =cost.benefit.axis,
                          col= "red",
                          confidence.intervals=plot.confidence.intervals,
                          standardize = standardize)
    }
    dev.off()
    par(cex=old_par$cex,lwd=old_par$lwd)
  }

  if(!all(plot.type %in% c("time","cluster","disperse"))){
    LuckyVerbose('O__O"...please select a right plot.type: "time","cluster","disperse" ')
  }

  ### save data
  saveRDS(L,file = paste0(dir_output,"/resultDCA.rds"))
  LuckyVerbose("(～￣▽￣)～ All done!")
  return(L)
}

