

#' @title ggplot style for decision estimate analysis
#' @description ggplot style for decision estimate analysis
#' @param load.existed.res whether to use existed rds result in the work space
#' @param plot.ntbft.ymin the ymin of plot.ntbft. Default is \code{ymin=-0.05}
#' @param plot.ntbft.ymax the ymax of plot.ntbft. Default is \code{ymax=max(nb[,c(nolines,nobands)],na.rm=T)}
#' @param plot.legend.position the position of plot legend
#' @param parameter.label the label of parameter label.
#' @inheritParams getDCA
#' @importFrom boot boot boot.ci
#' @importFrom reshape2 melt
#' @importFrom dplyr arrange
#' @details beta version
#' @return ggplot style plot and result
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' ## This is a simulative process and available only with CORRECT VARIABLES
#' @export
getDCA2 <- function(data,
                    group.by = NULL,
                    time.col = "time",
                    status.col = "event",
                    clusters,
                    time.raw.type=c("Day","Month","Year")[1],
                    time.target.type=c("Day","Month","Year")[3],
                    time.knot = c(3,5),
                    study.design = 'cohort',
                    population.prevalence = c(NA,NA),
                    load.existed.res = F,
                    plot.ntbft.by = c("group","model")[1],
                    plot.ntbft.ymin=NULL,
                    plot.ntbft.ymax=NULL,
                    width = 9,
                    height = 7,
                    add.bootc = T,
                    add.cv = F,
                    add.ci = T,
                    plot.legend.position="right",
                    plot.type = c("time","cluster","disperse")[1],
                    plot.label,
                    parameter.label = c("all"="All positive","none" = "All negative"),
                    plot.family = "Arial",
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
  path_resRDS <- paste0(dir_output,"/result_getDCA2_",names,".rds")
  logic <- length(list.files(path = dir_output,pattern = basename(path_resRDS))) == 0
  if(logic | !load.existed.res ){
    LuckyVerbose("RDS object of model isn't found. Create new one, please wait a minute... ")
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
      for(gn in group.name){ # gn = "proliferation" # gn = "NonGroup"

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
          L[[outcome.i_name]][[gn]][[cluster.name]][["data"]] <- data.i


          ## DCA data
          LuckyVerbose("=======")
          LuckyVerbose(cluster.name,": DCA for ",outcome.i_name,"...")
          data.i_2 <- data.i[,-grep("time|status",colnames(data.i))]
          ppi <- getPopPre(data.i_2,population.prevalence[i])

          ## Calculate Net benefit
          fit <- getGlmFormula(data.i_2,outcome.i_name)
          LuckyVerbose("Calculate Net benefit...")
          nb.i <- ntbft(data = data.i_2,
                        outcome = outcome.i_name,
                        frm=fit,
                        exterdt=NULL,
                        pred=NULL,
                        xstart=0.01,xstop=0.99,step=0.01,
                        type="treated",
                        study.design = study.design,
                        population.prevalence = ppi)
          L[[outcome.i_name]][[gn]][[cluster.name]][["ntbft"]] <- nb.i
          L[[outcome.i_name]][[gn]][[cluster.name]][["fit"]] <- fit

          ## Bootstrap method to correct overfitting
          #library(boot)
          LuckyVerbose("Bootstrap method to correct overfitting...")
          rstls <- boot(data=data.i_2,
                        statistic=diffnet,
                        R=500,
                        outcome = outcome.i_name,
                        frm=getGlmFormula(data.i_2,outcome.i_name),
                        xstart=0.01,xstop=0.99,step=0.01,
                        type="treated",
                        study.design = study.design,
                        population.prevalence = ppi)
          L[[outcome.i_name]][[gn]][[cluster.name]][["bootstrap"]] <- rstls

          ## Cross validation to correct overfitting
          LuckyVerbose("Cross validation to correct overfitting...")
          set.seed(2019)
          cv <- ntbft.cv(data=data.i_2,
                         outcome = outcome.i_name,
                         frm=getGlmFormula(data.i_2,outcome.i_name),
                         n_folds=10,R=500,
                         xstart=0.01,xstop=0.99,step=0.01,
                         type="treated",
                         study.design = study.design,
                         population.prevalence = ppi)
          L[[outcome.i_name]][[gn]][[cluster.name]][["cv"]] <- cv


          ## Confidence interval for the decision curve
          LuckyVerbose("Confidence interval for the decision curve...")
          set.seed(126)
          boot.cfint<-boot(data=data.i_2,
                           statistic=boot.confint,
                           R=500,
                           outcome = outcome.i_name,
                           frm=getGlmFormula(data.i_2,outcome.i_name),
                           xstart=0.01,xstop=0.99,step=0.01,
                           type="treated",
                           study.design = study.design,
                           population.prevalence = ppi)
          nbci<-NULL
          for(i in 1:length(boot.cfint$t0)){
            nbci<-rbind(nbci,boot.ci(boot.cfint,
                                     type="perc",index=i)$percent)
          }
          L[[outcome.i_name]][[gn]][[cluster.name]][["nbci"]] <- nbci

        }
      }


    }
    saveRDS(L,path_resRDS)
  } else {
    LuckyVerbose("RDS object of model found！Load it.")
    L <- readRDS(path_resRDS)
  }

  ###==========================get DCA plot======================###
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
          plotData[[nameTime[i]]][[nGC]]<- L[[nameTime[i]]][[nameGroup[j]]][[nameCluster[k]]]
        }
      }
    }

    ## DCA plot
    cairo_pdf(paste0(dir_output,"/ggDCA_bytime_plot.ntbft.by_",plot.ntbft.by,"_",names,".pdf"),width = width,height = height,family = plot.family,onefile = T)
    for(i in 1:length(nameTime)){ # i=1
      plotData_i <- plotData[[i]]

      ## merge data
      dca.i <- extraNtbif(plotData_i)

      ## ggplot for DCA plot
      p.i <- plot.ntbft(nb = dca.i,
                        group="group",
                        model="model",
                        plot.type = "time",
                        by = plot.ntbft.by,
                        plot.label = plot.label,
                        add.bootc = add.bootc,
                        add.cv = add.cv,
                        add.ci = add.ci,
                        ymin = plot.ntbft.ymin,
                        ymax = plot.ntbft.ymax,
                        nolines=2:4,
                        nobands=c(7,8),
                        legpos=plot.legend.position,
                        parameter.label = parameter.label)
      p.i <- p.i + ggtitle(paste0(Fastextra(nameTime[i],"_",2),"-",time.target.type)) + theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 18)); print(p.i)
    }
    dev.off()
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
          nGT <- paste0(nameGroup[j],'--',nameTime[i])
          plotData[[nameCluster[k]]][[nGT]]<- L[[nameTime[i]]][[nameGroup[j]]][[nameCluster[k]]]
        }
      }
    }
    names(plotData) <- getDCALabel(names(plotData),plot.label)


    ## DCA plot
    cairo_pdf(paste0(dir_output,"/ggDCA_bycluter_plot.ntbft.by_",plot.ntbft.by,"_",names,".pdf"),width = width,height = height,family = plot.family,onefile = T)
    for(i in 1:length(nameCluster)){ # i=1
      plotData_i <- plotData[[i]]

      ## merge data
      dca.i <- extraNtbif(plotData_i)

      ## ggplot for DCA plot
      p.i <- plot.ntbft(nb = dca.i,
                        group="group",
                        model="model",
                        plot.type = "cluster",
                        by = plot.ntbft.by,
                        #by = "model",
                        plot.label = plot.label,
                        add.bootc = add.bootc,
                        add.cv = add.cv,
                        add.ci = add.ci,
                        ymin = plot.ntbft.ymin,
                        ymax = plot.ntbft.ymax,
                        nolines=2:4,
                        nobands=c(7,8),
                        legpos=plot.legend.position,
                        parameter.label = parameter.label)

      p.i <- p.i + ggtitle(names(plotData)[i]) + theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 18)) ; print(p.i)
    }
    dev.off()
  }

  ## others
  if("disperse" %in% plot.type){
    LuckyVerbose("=======")
    LuckyVerbose("ㄟ(▔,▔)ㄏ : indivitual ploting...")

    ## get data
    plotData <- list()
    for(i in 1:length(nameTime)){ # i=1
      for(j in 1:length(nameGroup)){
        for(k in 1:length(nameCluster)){
          nGCT <- paste0(nameGroup[j],'--',nameCluster[k],'--',nameTime[i])
          plotData[[nGCT]]<- L[[nameTime[i]]][[nameGroup[j]]][[nameCluster[k]]]
        }
      }
    }

    ## DCA plot
    cairo_pdf(paste0(dir_output,"/ggDCA_bydisperse_",names,".pdf"),width = width,height = height,family = plot.family,onefile = T)
    for(i in 1:length(plotData)){ # i=1
      plotData_i <- plotData[i]

      ## merge data
      dca.i <- extraNtbif(plotData_i)
      dca.i$model <- getDCALabel(as.character(dca.i$model),plot.label)

      ## ggplot for DCA plot
      p.i <- plot.ntbft(nb = dca.i,
                        group="group",
                        model="model",
                        plot.type = "cluster",
                        by = "model",
                        plot.label = plot.label,
                        add.bootc = add.bootc,
                        add.cv = add.cv,
                        add.ci = add.ci,
                        ymin = plot.ntbft.ymin,
                        ymax = plot.ntbft.ymax,
                        nolines=2:4,
                        nobands=c(7,8),
                        legpos=plot.legend.position,
                        parameter.label = parameter.label)
      n.i <- names(plotData)[i]
      nG <- getGroupName(Fastextra(n.i,"--",1))
      nT <- gsub("outcome_","",Fastextra(n.i,"--",3)) %>% gsub("_","-",.)
      nC <- getDCALabel(Fastextra(n.i,"--",2),plot.label)
      title.i <- ifelse(is.null(nG),paste0(nT," ",nC),paste0(nG,nT," ",nC))
      p.i <- p.i + ggtitle(title.i) + theme(plot.title = element_text(hjust = 0.5,face = "bold",size = 18)); print(p.i)
    }
    dev.off()
   }

  if(!all(plot.type %in% c("time","cluster","disperse"))){
    LuckyVerbose('O__O"...please select a right plot.type: "time","cluster","disperse" ')
  }

  ###=========================Output======================###
  return(L)

  }



