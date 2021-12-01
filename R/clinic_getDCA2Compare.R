
#' @title Comparision of getDCA2 result
#' @description Comparision of getDCA2 result
#' @param object the result of \code{\link{getDCA2}}
#' @param glop_start the start of threshold probability
#' @param glop_end the end of threshold probability
#' @param glop_step the step of threshold probability
#' @param seed Numeric with 2 values. The first seed is for the boot of \code{nbdiff}, and the second one is for the boot of \code{areadiff}.
#' @inheritParams getDCA2
#' @seealso \code{\link{getDCA2}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @export
getDCA2Compare <- function(object,
                           glop_start = 0.01,
                           glop_end = 0.1,
                           glop_step = 0.01,
                           seed = c(127,128),
                           names="project"){

  ## output dir
  dir_output <- "./decisionCurveAnalysis"
  dir.create(dir_output,recursive = T,showWarnings = F)

  ## get names of time,group and cluster
  nameTime <- names(object)
  nameGroup <- names(object[[1]])
  nameCluster <- names(object[[1]][[1]])

  ## get compare data for every time-group combination
  compareData <- list()
  for(i in 1:length(nameTime)){
    for(j in 1:length(nameGroup)){
      nGT <- paste0(nameTime[i],"--",nameGroup[j])
      compareData[[nGT]] <- object[[nameTime[i]]][[nameGroup[j]]]
    }
  }

  ## calculate for every model
  x <- list()
  for(i in 1:length(compareData)){ # i=1
    t.i <- names(compareData)[i]
    n.i2 <- getGroupName(Fastextra(t.i,"--",2))
    n.i1 <- gsub("outcome_","",Fastextra(t.i,"--",1))
    n.i <- ifelse(is.null(n.i2),n.i1,paste0(n.i1,"--",n.i2))
    LuckyVerbose("==========",t.i,"==========")
    x[[t.i]] <- getDCACompare_one(compareData[[i]],
                                  glop_start = glop_start,
                                  glop_end = glop_end,
                                  glop_step = glop_step,
                                  seed = seed,
                                  name = n.i)
  }

  ## output data
  l <- list(
    Parameter = list(
      glop_start = glop_start,
      glop_end = glop_end,
      glop_step = glop_step,
      seed = seed,
      names=names
    ),
    Data = x
  )
  saveRDS(l,paste0(dir_output,"/result_getDCA2Compare_",names,"_glop_start-",glop_start,"_glop_end-",glop_end,".rds"))
  return(l)
}

