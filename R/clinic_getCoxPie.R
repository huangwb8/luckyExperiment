


#' @title Get pie plot of Cox chisq proportion for variants in Cox regression model
#' @description Get pie plot of Cox chisq proportion for variants in Cox regression model
#' @param pie.color a character vector with color value and variable names. Default is missing.
#' @param Others.color the color of others. It works when \code{global=T}
#' @param legend.name the name of legened
#' @param legend.position the position of legend
#' @param is.title whether plot title
#' @param plot.size the size of saved plots
#' @param dc the digit counts of proportion
#' @inheritParams getDCA
#' @inheritParams survival::cox.zph
#' @importFrom ggplot2 ggplot aes geom_bar coord_polar scale_fill_manual labs theme element_text element_blank
#' @importFrom dplyr arrange
#' @importFrom survival coxph cox.zph
#' @return pie data and plots
#' @seealso \code{\link[survival]{cox.zph}}; \code{\link[survival]{cox.ph}}
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @examples
#' data("myeloma",package = "survminer");
#' data <- myeloma;rm(myeloma)
#' show.color(mycolor)
#' pie.color <- c(TP53="#FB8072", WHSC1="#BEBADA", DEPDC1="#80B1D3", Others="#A65628")
#' pie.label = c("TP53" = "tp53","WHSC1" = "whsc1","DEPDC1"="depdc1")
#' a <- getCoxPie(data,
#'                time.col = "time",
#'                status.col = "event",
#'                clusters = list(c("TP53","WHSC1"),
#'                                c("TP53","WHSC1","DEPDC1")),
#'                time.raw.type = c("Day", "Month","Year")[2],
#'                time.target.type = c("Day", "Month", "Year")[3],
#'                global=T,
#'                pie.color = pie.color,
#'                pie.label = pie.label,
#'                plot.size = 7,
#'                names = "project")
#' @export
getCoxPie <- function(data,
                      time.col = "time",
                      status.col = "event",
                      clusters,
                      time.raw.type = c("Day", "Month","Year")[1],
                      time.target.type = c("Day", "Month", "Year")[3],
                      global=T,
                      pie.color,
                      Others.color = mycolor[27],
                      pie.label,
                      dc = 1,
                      legend.name = "Variable",
                      legend.position = "right",
                      is.title = T,
                      plot.size = 7,
                      names = "project"){

  ### grobal option
  old_par <- par()
  dir_output <- "./CoxPiePlot"
  dir.create(dir_output,recursive = T,showWarnings = F)

  ### time convertion
  data.x <- as.data.frame(data,stringsAsFactors = F)
  # getDCA::convert.time
  data.x[,time.col] <- convert.time(data.x[,time.col],
                                    from.ts = time.raw.type,
                                    to.ts = time.target.type)
  # report
  #LuckyVerbose(paste0("The max of time is ",round(max(data.x[,time.col]),2)," ",time.target.type,"s,and the min is ",round(min(data.x[,time.col]),2)," ",time.target.type,"s. Note: time.knot must be in the interval!"))

  ### build cox regression model

  # pdf name
  if(global){
    nPDF <- paste0(dir_output,"/CoxPiePlot_global_",names,".pdf")
  } else {
    nPDF <- paste0(dir_output,"/CoxPiePlot_",names,".pdf")
  }

  # pie color
  if(missing(pie.color)){
    r <- c(4,10,3,5);pie.color <- mycolor[c(r,setdiff(1:length(mycolor),r))]
    # fill colours
    var <- unique(unlist(clusters))
    fill_color <- pie.color[1:length(var)];
    names(fill_color) <- var
    fill_color <- c(fill_color,Others.color);names(fill_color)[length(fill_color)] <- "Others"
  } else {
    fill_color <- pie.color
  }

  cairo_pdf(nPDF,plot.size,plot.size,family = "Arial",onefile = T)
  l <- list()
  for(i in 1:length(clusters)){ # i=1
    c.i <- clusters[[i]]
    nC <- paste0(c.i,collapse = "_")
    data.i <- data.x[,c(time.col,status.col,c.i)]
    colnames(data.i)[1:2] <- c("time","event")

    ## test the proportional hazards assumption for a Cox regression model fit
    fit <- coxph(getCoxfomula(data.i,c.i),data = data.i)
    res.fit <- cox.zph(fit, transform="km", global=T);
    df <- as.data.frame(res.fit$table,stringsAsFactors = F)
    l[[nC]][["fit"]] <- res.fit

    ## pie data
    if(global){
      LuckyVerbose(nC,": Use global chisq as whole...")
      Others <- df$chisq[rownames(df) %in% "GLOBAL"] - sum(df$chisq[!rownames(df) %in% "GLOBAL"])
      df.pie <- data.frame(Var = "Others", Value = Others)
      for(i in 1:(nrow(df)-1)){
        df.pie.i <- data.frame(Var = rownames(df)[i],
                               Value = df$chisq[i])
        df.pie <- rbind(df.pie,df.pie.i)
      }
    } else {
      LuckyVerbose(nC,": Use sum of variant chisq as whole...")
      df.pie <- NULL
      for(i in 1:(nrow(df)-1)){
        df.pie.i <- data.frame(Var = rownames(df)[i],
                               Value = df$chisq[i])
        df.pie <- rbind(df.pie,df.pie.i)
      }
    }
    df.pie$Por <- df.pie$Value/sum(df.pie$Value)
    df.pie <- arrange(df.pie,desc(Value))
    l[[nC]][["data"]] <- df.pie

    ## ggplot2 pie

    # title
    if(is.title){
      title <- paste0(c.i,collapse = "+")
    } else {
      title <- NULL
    }

    # fill scale
    rVar <- factorPieData(df.pie,c.i) # print(rVar)
    fc <- fill_color[Fastmatch(rVar,names(fill_color))]
    # pie label
    if(missing(pie.label)){
      myLabel <- paste(rVar,
                       paste0(round(df.pie$Por[Fastmatch(rVar,df.pie$Var)]*100,dc),"%"),
                       sep = ": ")
    } else {
      myLabel <- paste(pie.label[Fastmatch(rVar,names(pie.label))],
                       paste0(round(df.pie$Por[Fastmatch(rVar,df.pie$Var)]*100,dc),"%"),
                       sep = ": ")
    }

    # pie plot
    p <- ggplot(df.pie,aes(x="",y=Value,fill=Var)) +
      geom_bar(stat = "identity") +
      coord_polar(theta = "y") +
      scale_fill_manual(breaks = rVar,
                        values = fc,
                        labels = myLabel,
                        name = legend.name) +
      labs(x = "", y = "", title = title) +
      theme(plot.title = element_text(hjust = 0.5,face = "bold",size=15),
            legend.text = element_text(face = "bold",size=12),
            legend.title = element_text(face = "bold",size=12),
            legend.position = legend.position,
            axis.text.x = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank())
    print(p)
  }
  dev.off()

  ## output data
  saveRDS(l,paste0(dir_output,"/CoxPiePlot_",names,".rds"))
  LuckyVerbose("All done!")
  return(l)
}

###====================Assistant function========================###

## get cox fomula for every cluster
#' @export
getCoxfomula <- function(data.i,
                         c.i){

  f1 <- paste0("Surv(time,event) ~ ",paste(c.i,collapse = " + "))
  f2 <- as.formula(f1)
  return(f2)
}

## specify level in fill scale
factorPieData <- function(df.pie,c.i){
  var <- as.character(df.pie[,"Var"])
  if("Others" %in% var){
    l <- c(c.i,"Others")
  } else {
    l <- c.i
  }
  return(l)
}























