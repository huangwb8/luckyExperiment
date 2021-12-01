
## Assistant sub-functions for getDCA series functions

####==========================getDCA=========================####


## get outcome data
getOutcome <- function(vt,
                       time="time",
                       status="status",
                       knot.i = 3){

  ## test
  if(F){
    vt <- data1[1,]
  }

  ## get target data
  time.i <- as.numeric(vt[names(vt) %in% time])
  status.i <- as.numeric(vt[names(vt) %in% status])

  if(status.i == 0){
    ## 无论时间是否长于time knot，由于event未发生，因此计为0
    return(0)
  } else {
    ## event已经发生。此时根据时间来判断
    if(time.i < knot.i){
      # event在time knot内发生，计为1
      return(1)
    } else {
      # event在time knot之后发生，计为0
      return(0)
    }
  }
  ## End
}

## time convertion
#1 year = (12+1/6) month = 365 day;1 month = 30 day
convert.time <- function(time,from.ts,to.ts){
  if(from.ts == "Day"){
    day <- time
    month <- time/30
    year <- time/365
  } else {
    if(from.ts == "Month"){
      day <- time*30
      month <- time
      year <- time/(12+1/6)
    } else {
      if(from.ts == "Year"){
        day <- time*365
        month <- time*(12+1/6)
        year <- year
      } else {
        stop("Error!not a right type of time!")
      }
    }
  }
  #输出结果
  time1 <- data.frame(
    Day=day,
    Month=month,
    Year=year
  )
  time2 <- time1[,to.ts]
  return(time2)
}

## get glm formula
getGlmFormula <- function(data.i_2,
                          outcome.i_name){
  others <- setdiff(colnames(data.i_2),outcome.i_name)

  f1 <- paste0(outcome.i_name," ~ ",paste0(others,collapse = " + "))
  f2 <- as.formula(f1)
  return(f2)
}


## get label for DCA plot
getDCALabel <- function(label,plot.label){
  # label=c("TP53_WHSC1","TP53_WHSC1_DEPDC1")
  n <- NULL
  for(i in 1:length(label)){ # i=1
    label.i <- sort(Fastextra(label[i],"_"))
    for(j in 1:length(plot.label)){ # j=1
      label.j <- sort(Fastextra(names(plot.label)[j],"-"))
      if(identical(label.i,label.j)){
        # completely matched
        n <- c(n,plot.label[[j]])
      }
    }
  }
  return(n)

}

getDCALabel2 <- function(label,plot.label){
  # label = names(plotData_i)
  label1 <- Fastextra(label,"--",1)
  label2 <- Fastextra(label,"--",2)
  label2.2 <- getDCALabel(label2,plot.label)
  return(paste(label1,label2.2,sep=": "))
}

## get group name
getGroupName <- function(gn){

  if(gn %in% c("NoneGroup","NonGroup")){
    return(NULL)
  } else {
    return(paste0(gn,": "))
  }
}

## estimate population.prevaluence from cohort
getPopPre <- function(data.i,population.prevalence.i){
  if(is.na(population.prevalence.i)){
    vt <- data.i[,grep("outcome",colnames(data.i))]
    r <- 1.1 * (sum(vt == 1)/length(vt))
  } else {
    r <- population.prevalence.i
  }
  return(r)
}

## romove NonGroup label from cureve.names
rmNonGroupLabel <- function(curve.names2){
  if(is.one.true(grepl("NonGroup",curve.names2))){
    x <- gsub("NonGroup: ","",curve.names2)
    return(x)
  } else {
    return(curve.names2)
  }
}


####=======================getDCA2=========================####

## Net benefit function
ntbft<-function(data,
                outcome,
                frm=NULL,exterdt=NULL,
                pred=NULL,xstart=0.01,
                xstop=0.99,step=0.01,
                type="treated",
                study.design = 'cohort',
                population.prevalence = 0.3){
  pt<-seq(from=xstart,to=xstop,by=step)
  lpt<-length(pt)
  if(type=="treated") coef<-cbind(rep(1,lpt),rep(0,lpt))
  if(type=="untreated") coef<-cbind(rep(0,lpt),rep(1,lpt))
  if(type=="overall") coef<-cbind(rep(1,lpt),rep(1,lpt))
  if(type=="adapt") coef<-cbind(1-pt,pt)
  response<-as.vector(t(data[outcome]))
  if(is.data.frame(exterdt)) response<-as.vector(t(exterdt[outcome]))

  ## event.rate
  if(study.design == 'cohort'){
    event.rate <- mean(response)
  } else {
    event.rate <- population.prevalence
  }

  ## Pred value
  nball<-event.rate-(1-event.rate)*pt/(1-pt)
  nbnone<-1-event.rate-event.rate*(1-pt)/pt
  if(is.null(pred)){
    model<-glm(frm,data=data,family=binomial("logit"))
    pred<-model$fitted.values
    if(is.data.frame(exterdt))
      pred<-predict(model,newdata=exterdt,type="response")
  }

  # pred and response should be of the same length
  N<-length(pred)
  nbt<-rep(NA,lpt)
  nbu<-rep(NA,lpt)
  for(t in 1:lpt){
    tp<-sum(pred>=pt[t] & response==1)
    fp<-sum(pred>=pt[t] & response==0)
    fn<-sum(pred<pt[t] & response==1)
    tn<-sum(pred<pt[t] & response==0)
    nbt[t]<-tp/N-fp/N*(pt[t]/(1-pt[t]))
    nbu[t]<-tn/N-fn/N*((1-pt[t])/pt[t])
  }
  nb<-data.frame(pt)
  names(nb)<-"threshold"
  nb["all"]<-coef[,1]*nball
  nb["none"]<-coef[,2]*nbnone
  nb["pred"]<-coef[,1]*nbt+coef[,2]*nbu
  return(nb)
}


## Bootstrap method to correct overfitting
diffnet<-function(data,ii,
                  outcome,frm,
                  xstart=0.01,xstop=0.99,
                  step=0.01,type="treated",
                  study.design = 'cohort',
                  population.prevalence = 0.3){
  dd<-data[ii,]
  nb<-ntbft(data=dd,outcome=outcome,frm=frm,
            xstart=xstart,xstop=xstop,
            step=step,type=type,
            study.design = study.design,
            population.prevalence = population.prevalence)
  nb0<-ntbft(data=dd,outcome=outcome,frm=frm,
             exterdt=data,xstart=xstart,
             xstop=xstop,step=step,type=type,
             study.design = study.design,
             population.prevalence = population.prevalence)
  diff<-nb$pred-nb0$pred
  return(diff)
}


## Cross validation to correct overfitting
ntbft.cv <- function(data,outcome,frm,
                     n_folds=10,R=200,
                     xstart=0.01,xstop=0.99,
                     step=0.01,type="treated",
                     study.design = 'cohort',
                     population.prevalence = 0.3){
  cv<-NULL
  for(i in 1:R){
    n_train<-nrow(data)
    folds_i<-sample(rep(1:n_folds,length.out=n_train))
    pred<-rep(NA,n_train)
    for(k in 1:n_folds){
      test_i<-which(folds_i==k)
      train_dt<-data[-test_i,]
      test_dt<-data[test_i,]
      model<-glm(frm,data=train_dt,family=binomial("logit"))
      pred[test_i]<-predict(model,newdata=test_dt,type="response")
    }
    cv<-cbind(cv,ntbft(data=data,outcome=outcome,
                       pred=pred,xstart=xstart,xstop=xstop,
                       step=step,type=type,
                       study.design = study.design,
                       population.prevalence = population.prevalence)$pred)
  }
  return(rowMeans(cv))
}


## Confidence interval for the decision curve
boot.confint <- function(data,ii,
                         outcome,frm,
                         xstart=0.01,xstop=0.99,
                         step=0.01,type="treated",
                         study.design = 'cohort',
                         population.prevalence = 0.3){
  dd<-data[ii,]
  nb<-ntbft(data=dd,outcome=outcome,frm=frm,
            xstart=xstart,xstop=xstop,
            step=step,type=type,
            study.design = study.design,
            population.prevalence = population.prevalence)
  return(nb$pred)
}


## extra multiple ntbif data
#' @export
extraNtbif <- function(plotData_i){

  df <- NULL
  for(i in 1:length(plotData_i)){ # i=1


    p.i <- plotData_i[[i]]

    ## Basedata
    df.i <- p.i$ntbft
    group.i <- Fastextra(names(plotData_i)[i],"--",1)
    model.i <- Fastextra(names(plotData_i)[i],"--",2)

    ## Bootstrap prediction
    pred.bootc <- df.i$pred - rowMeans(t(p.i$bootstrap$t))

    ## Cross validation prediction
    pred.cv <- p.i$cv

    ## Confidence interval
    # View(p.i$nbci)
    low <- p.i$nbci[,4]
    up <- p.i$nbci[,5]

    ## Merge all data to a data frame
    df.i2 <- cbind(df.i,
                   pred.bootc = pred.bootc,
                   pred.cv = pred.cv,
                   low = low,
                   up = up,
                   group=group.i,
                   model=model.i)
    df <- rbind(df,df.i2)
  }
  return(df)

}


## get DCA plot
plot.ntbft <- function(nb,
                       group="group",
                       model="model",
                       plot.type = c("time","cluster")[1],
                       by = c("group","model")[1],
                       plot.label,
                       add.bootc = T,
                       add.cv = F,
                       add.ci = T,
                       nolines=2:4,
                       nobands=c(7,8),
                       ymin=NULL,
                       ymax=NULL,
                       legpos="right",
                       parameter.label){

  ###====================change model names======================###
  # need: getDCA_assistant function: getDCA

  if(is.null(ymin)) ymin = -0.05
  if(is.null(ymax)) ymax = max(nb[,c(nolines,nobands)],na.rm=T)

  colnames(nb)[Fastmatch(c("all","none"),colnames(nb))] <- as.character(parameter.label[c("all","none")])

  # cluster name have to be changed based on plot.label
  if(plot.type == "time"){
    nb$model <- as.character(nb$model)
    for(i in unique(nb$model)){
      nb$model[nb$model %in% i] <- getDCALabel(i,plot.label)
    }
    cf.name = "Model"
  }

  if(plot.type == "cluster"){
    nb$model <- as.character(nb$model)
    for(i in unique(nb$model)){
      nb$model[nb$model %in% i] <- gsub("outcome_","",i)
    }
    cf.name = "Time"

  }

  ## select data
  # colnames(nb) |pred.bootc|pred.cv
  id.vars=c("threshold","group","model")
  expr <- "^low$|^up$"
  if(!add.bootc) expr <- paste0(expr,"|^pred.bootc$")
  if(!add.cv) expr <- paste0(expr,"|^pred.cv$")
  nb.melt <- reshape2::melt(nb[,-grep(expr,colnames(nb))],
                            id.vars= id.vars,
                            value.name="Netbenefit",
                            variable.name="Parameters")
  nb.melt <- dplyr::arrange(nb.melt,group,model)

  #View(filter(nb.melt,Parameters %in% "all",model %in% "complexModel"))
  #View(filter(nb.melt,Parameters %in% "all",model %in% "simpleModel"))

  ## ggplot

  # line
  lty_limits = c(as.character(parameter.label["all"]),as.character(parameter.label["none"]),"pred","pred.bootc","pred.cv");
  lty_values = c(2,2,1,3,4)
  lty_size=c(1,1.2,1.2,1,1)
  s <- lty_limits %in% unique(as.character(nb.melt$Parameters))
  lty_limits <- lty_limits[s]
  lty_values <- lty_values[s]
  lty_size <- lty_size[s]

  # group by by="group" by = "model"
  if(by=="group"){
    # confidence interval
    if(!add.ci){
      ribbon <- NULL
      lab <- labs(x = "Threshold probability",y = "Net benefit",colour=cf.name)
    } else {
      ribbon <- geom_ribbon(data=nb,aes(x = threshold,ymin = low,ymax = up,fill = model),linetype=2,alpha=0.1)
      lab <- labs(x = "Threshold probability",y = "Net benefit",colour=cf.name ,fill = cf.name )
    }

    # expression of geom_line
    expr_geom_line <- geom_line(aes(x=threshold,y=Netbenefit,
                                    colour = model,
                                    linetype = Parameters,
                                    size = Parameters))

    # expression of facet
    if(length(unique(nb.melt$group))>1){
      expr_facet <- facet_grid(. ~ group)
    } else {
      expr_facet <- NULL
    }

    # all & none
    if(plot.type == "time"){
      # 由于按group分组，所以不同的model必然共享一个none & all线
      x.an <- filter(nb.melt,Parameters %in% as.character(parameter.label))
      # expr_an <- ggplot(x.an,aes(x=threshold,y=Netbenefit,linetype = Parameters,size = Parameters)) +  geom_line(colour = "black")
      expr_an <- geom_line(data = x.an,aes(x=threshold,y=Netbenefit,linetype = Parameters,size = Parameters),colour = "black")
    }
    if(plot.type == "cluster"){
      # 由于按group分组，所以不同的model必然共享一个none线
      x.an <- filter(nb.melt,Parameters %in% as.character(parameter.label["none"]))
      # expr_an <- ggplot(x.an,aes(x=threshold,y=Netbenefit,linetype = Parameters,size = Parameters)) +  geom_line(colour = "black")
      expr_an <- geom_line(data = x.an,aes(x=threshold,y=Netbenefit,linetype = Parameters,size = Parameters),colour = "black")
    }

  }
  if(by=="model"){
    if(length(unique(nb.melt$group)) >1){

      ## multiple groups

      # confidence interval
      if(!add.ci){
        ribbon <- NULL
        lab <- labs(x = "Threshold probability",y = "Net benefit",colour="Group")
      } else {
        ribbon <- geom_ribbon(data=nb,aes(x = threshold,ymin = low,ymax = up,fill = group),linetype=2,alpha=0.1)
        lab <- labs(x = "Threshold probability",y = "Net benefit",colour="Group",fill = "Group")
      }

      # expression of geom_line
      expr_geom_line <- geom_line(aes(x=threshold,y=Netbenefit,colour = group,linetype = Parameters,size = Parameters))

    } else {

      ## Single groups

      # confidence interval
      if(!add.ci){
        ribbon <- NULL
        lab <- labs(x = "Threshold probability",y = "Net benefit")
      } else {
        ribbon <- geom_ribbon(data=nb,aes(x = threshold,ymin = low,ymax = up),fill = mycolor[4],linetype=2,alpha=0.1)
        lab <- labs(x = "Threshold probability",y = "Net benefit")
      }

      # expression of geom_line
      expr_geom_line <- geom_line(aes(x=threshold,y=Netbenefit,linetype = Parameters,size = Parameters),colour = mycolor[4])
    }
    expr_facet <- facet_grid(. ~ model)

    # expression of facet
    if(length(unique(nb.melt$model))>1){
      expr_facet <- facet_grid(. ~ model)
    } else {
      expr_facet <- NULL
    }

    # all & none
    # 由于按model分组，所以不同的group通常不共享一个all线，但共享none线
    x.an <- filter(nb.melt,Parameters %in% as.character(parameter.label["none"]))
    # expr_an <- ggplot(x.an,aes(x=threshold,y=Netbenefit,linetype = Parameters,size = Parameters)) +  geom_line(colour = "black")
    expr_an <- geom_line(data = x.an,aes(x=threshold,y=Netbenefit,linetype = Parameters,size = Parameters),colour = "black")


  }

  # plot
  p <- ggplot(nb.melt) +
    expr_geom_line +
    scale_linetype_manual(limits = lty_limits,values = lty_values) +
    scale_size_manual(values=lty_size) +
    ribbon +  expr_an + lab +
    scale_y_continuous(limits=c(ymin,ymax))+
    expr_facet +
    theme_bw() +
    theme(axis.title = element_text(face = "bold",size = 15),
          axis.text = element_text(face = "bold",size = 12),
          legend.title = element_text(face = "bold",size = 15),
          legend.text =element_text(face = "bold",size = 12),
          legend.position = legpos,
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(face = "bold",size = 15))
  #panel.grid =element_blank(),

  return(p)

}


####==================getDCA2Compare====================####

#' @importFrom boot boot boot.ci
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme  element_blank element_text
getDCACompare_one <- function(l,
                              glop_start = 0.01,
                              glop_end = 0.1,
                              glop_step = 0.01,
                              seed = c(127,128),
                              name = "3-Year"){

  ## value names
  val.name <- unique(unlist(Fastextra(names(l),"_")))

  ## data containing all value names
  df <- NULL
  for(i in 1:length(l)){ # i=2
    data.i <- l[[i]]$data
    s.i <- grepl("outcome_",colnames(data.i)) | colnames(data.i) %in% val.name
    data.i2 <- data.i[s.i]
    if(is.null(df)){
      df <- data.i2
    } else {
      add.id <- !colnames(data.i2) %in% colnames(df)
      if(is.one.true(add.id)) df <- cbind(df,data.i2[add.id])
    }
  }
  outcome <- colnames(df)[grep("outcome_",colnames(df))]

  ## get two combination
  cb <- combn(names(l),2)
  compare <- list()
  for(i in 1:ncol(cb)){ # i=1
    cb.i <- cb[,i]
    model1 <- l[[cb.i[1]]]$fit
    model2 <- l[[cb.i[2]]]$fit
    n.cbi <- paste0(cb.i,collapse = " vs. ")

    LuckyVerbose("======")
    LuckyVerbose(n.cbi," :")

    ## boot diff
    set.seed(seed[1])
    boot.diff <- boot(data=df,
                      statistic=nbdiff,
                      R=500,
                      outcome=outcome,
                      frm1=model1,
                      frm2=model2,
                      xstart=0.01,
                      xstop=0.99,
                      step=0.01,
                      type="treated")

    ## get pvalue for every sub boot
    LuckyVerbose("Boot separate P value...")
    pvalue <- NULL
    for(i in 1:length(boot.diff$t0)){
      pvalue <- c(pvalue,mean(abs(boot.diff$t[,i]-boot.diff$t0[i])>abs(boot.diff$t0[i])))
    }
    portion <- paste0(round(100*mean(pvalue<=0.05),2),"%")

    ## plot p value
    LuckyVerbose("Get separate P value plot...")
    if(T){
      df.pvalue <- data.frame(x=0.01*(1:99),
                              y=pvalue,
                              cohort = ifelse(pvalue <=0.05,"P.val<=0.05","P.val>0.05"))
      p.i <- ggplot(df.pvalue,aes(x=x,y=y,colour = cohort)) +
        geom_point(size = 3) +
        labs(title = paste0(name," : ",n.cbi),
             x = "Boot Seeds",
             y = "P value",
             subtitle = paste0("p.value<=0.05"," ",portion)) +
        theme_bw() +
        theme(panel.grid =element_blank(),
              plot.title = element_text(face = "bold",size = 15,hjust = 0.5),
              plot.subtitle = element_text(face = "bold",size = 12,hjust = 0.5),
              axis.title = element_text(face = "bold",size = 15),
              axis.text = element_text(face = "bold",size = 12),
              legend.title = element_text(face = "bold",size = 12),
              legend.text =element_text(face = "bold",size = 12),
              legend.position = "bottom")
      win.graph(8,8);print(p.i)
    }

    ## get single p value
    LuckyVerbose("Boot global P value...")
    set.seed(seed[2])
    boot.area <- boot(data=df,
                      statistic=areadiff,
                      R=1000,
                      outcome=outcome,
                      frm1=model1,
                      frm2=model2,
                      xstart=glop_start,
                      xstop=glop_end,
                      step=glop_step,
                      type="treated")
    glopvalue <- mean(abs(boot.area$t-boot.area$t0)>abs(boot.area$t0))
    LuckyVerbose("The result shows that the P value for the difference in A-NBC between in the threshold range from",glop_start,"to",glop_end,"is",glopvalue,"(we have thus a significant difference between the two models on that threshold range if glopvalue < 0.05).",levels = 2)

    ## output data
    compare[[n.cbi]][["boot.diff"]] <- boot.diff
    compare[[n.cbi]][["pvalue.plot"]] <- p.i
    compare[[n.cbi]][["boot.area"]] <- boot.area
    compare[[n.cbi]][["glopvalue"]] <- glopvalue
  }

  ## output
  return(compare)
}

nbdiff<-function(data,ii,outcome,
                 frm1,frm2,
                 xstart=0.01,
                 xstop=0.99,step=0.01,
                 type="treated"){
  dd<-data[ii,]
  nb1<-ntbft(data=dd,outcome=outcome,frm=frm1,
             xstart=xstart,xstop=xstop,
             step=step,type=type)
  nb2<-ntbft(data=dd,outcome=outcome,frm=frm2,
             xstart=xstart,xstop=xstop,
             step=step,type=type)
  nb.diff<-nb2$pred-nb1$pred
  return(nb.diff)
}


areadiff<-function(data,ii,outcome,
                   frm1,frm2,
                   xstart=0.01,xstop=0.99,
                   step=0.01,type="treated"){
  dd<-data[ii,]
  nb1<-ntbft(dd,outcome=outcome,frm=frm1,
             xstart=xstart,xstop=xstop,
             step=step,type=type)
  nb2<-ntbft(dd,outcome=outcome,frm=frm2,
             xstart=xstart,xstop=xstop,
             step=step,type=type)
  area1<-0; area2<-0
  for(i in 1:(nrow(nb1)-1)){
    area1<-area1+(nb1$pred[i]+nb1$pred[i+1])*step/2
    area2<-area2+(nb2$pred[i]+nb2$pred[i+1])*step/2
  }
  return(area2-area1)
}



