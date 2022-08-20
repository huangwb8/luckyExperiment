




# 2022-01-14 Test S4 in R

# Reference: http://adv-r.had.co.nz/S4.html

# Reference of DOSE

# https://github.com/YuLab-SMU/DOSE/blob/master/R/00-AllClasses.R
# https://github.com/YuLab-SMU/DOSE/blob/master/R/AllGenerics.R

#' @name Animal-class
#' @title Animal-class
#' @docType class
#' @description Animal class
#' @slot Time Time
#' @slot Group Group
#' @slot ID ID
#' @slot weight weight
#' @slot tumorVolume tumorVolume
#' @author Weibin Huang<\email{654751191@@qq.com}>
#' @seealso \code{\link{animal}}
#' @exportClass Animal
#' @keywords classes
setClass(Class = "Animal",
         slot = (
           list(
             Time = "character",
             Group = "character",
             ID = "character",
             weight = "numeric",
             tumorVolume = "numeric"
           )
         ))

#' @title Create an Animal object from a data frame
#' @description Create an Animal object from a data frame
#' @param x a data frame
#' @param tumorSize a specified tumor size vector for animals. Always looks like: '7.15,8.56;4.85,9.65;'
#' @return Tumor volume for \code{\link{Animal-class}}
#' @author Weibin Huang<\email{654751191@@qq.com}
tumor_volume <- function(x,
                         tumorSize = "tumor_size_mm") {
  v <- x[, tumorSize]

  tv <- NULL

  for (i in 1:length(v)) {
    # i=86

    v.i <- gsub(' ', '', v[i]) # (8.58*7.61^2 + 4.85*4.55^2)*0.5
    v.i2 <- Fastextra(v.i, split = ';')
    volume_i2 <- 0
    for (para in v.i2) {
      # para = "7.61,8.58"
      p.j <- as.numeric(Fastextra(para, split = ','))

      # volume = [1/2 x (the major axis) x (the minor axis)^2]
      volume_i2 <- volume_i2 + 0.5 * max(p.j) * (min(p.j)) ^ 2
    }

    tv <- c(tv, volume_i2)

  }

  return(tv)

}


#' @title Create an Animal object from a data frame
#' @description Create an Animal object from a data frame
#' @inheritParams tumor_volume
#' @param Time Time.Character.
#' @param Group The group of animals
#' @param ID The ID of animals
#' @param weight The weight of animals
#' @return \code{\link{Animal-class}}
#' @author Weibin Huang<\email{654751191@@qq.com}
#' @export
animal <- function(x,
                   Time = "Time",
                   Group = "Group",
                   ID = "ID",
                   tumorSize = "tumor_size_mm",
                   weight = "weigth_g") {
  targetCol <- c(Time, Group, ID, weight)

  x2 <- x[targetCol]

  colnames(x2) <- c('Time', 'Group', 'ID', 'weight')

  x2$ID <- as.character(x2$ID)

  x2$tumorVolume <- tumor_volume(x, tumorSize = tumorSize)

  # class(x2) <- 'Animal'

  # Make a new Animal object
  x3 <- new(
    'Animal',
    Time = as.character(x2$Time),
    Group = as.character(x2$Group),
    ID = as.character(x2$ID),
    weight = x2$weight,
    tumorVolume = x2$tumorVolume
  )
  return(x3)

}


#' @rdname Animal-method
#' @name Animal-method
#' @title Animal-method
#' @description \code{as.data.frame} method for \code{Animal} class
#' @param object \code{\link{Animal-class}}
#' @return a data frame
#' @seealso \code{\link{Animal}}
#' @author Weibin Huang<\email{654751191@@qq.com}
#' @export
as.data.frame.Animal <- function(object) {
  d <- data.frame(
    Time = object@Time,
    Group = object@Group,
    ID = object@ID,
    weight = object@weight,
    tumorVolume = object@tumorVolume,
    stringsAsFactors = F
  )
  return(d)
}

#' @rdname Animal-method
#' @name Animal-method
#' @description \code{plot} method for \code{Animal} class
#' @param type One of \code{tumorVolume} and \code{weight}
#' @param list.color A color vector with names of \code{Group} slot.
#' @param plot.size The size of the \code{ggplot} plot.
#' @param verbose Whether report process and results
#' @inheritParams ggplot2::theme
#' @inheritParams ggrepel::geom_label_repel
#' @import ggplot2
#' @importFrom ggrepel geom_label_repel
#' @importFrom graphics plot
#' @importFrom cowplot plot_grid
#' @examples
#' list.color <- mycolor[1:7]
#' names(list.color) <- paste('G',1:7,sep = '')
#' @export
plot.Animal <- function(object,
                        type = c('tumorVolume', 'weight')[1],
                        list.color = {
                          list.color <- mycolor[1:7]
                          names(list.color) <- paste('G', 1:7, sep = '')
                          list.color
                        },
                        legend.position = 'right',
                        label.size = 0.15,
                        plot.size = 25,
                        verbose = T) {
  d <- as.data.frame(object)

  # Data preparation
  all_type <- c('weight', 'tumorVolume')
  if (type %in% all_type) {
    if (verbose) message('Plot ', type, '...')
    d2 <- na.omit(d[c("Time", "Group", "ID", type)])
    colnames(d2)[ncol(d2)] <- 'type'
    targetTime <-
      ifelse(is.null(levels(d2$Time)), rev(unique(d2$Time))[1], rev(levels(d2$Time))[1])
    d2$Label <- ifelse(d2$Time %in% targetTime, d2$ID, NA)
    d2$Group <-
      factor(d2$Group, levels = names(list.color))

    # line charts
    p1 <- ggplot(d2, aes(
      x = Time,
      y = type,
      color = Group,
      group = ID
    )) +
      geom_line(size = 1.5) +
      geom_point(size = 2) +
      scale_color_manual(breaks = names(list.color),
                         values = as.character(list.color)) +
      geom_label_repel(
        aes(label = Label),
        label.size = label.size,
        nudge_x = 1,
        show.legend = F,
        na.rm = T
      ) +
      labs(y = type) +
      theme_bw() +
      theme(
        plot.title = element_text(
          size = (plot.size / 20) * 18,
          face = "bold",
          hjust = 0.5
        ),
        axis.title = element_text(size = (plot.size / 20) * 15, face = "bold"),
        axis.text = element_text(size = (plot.size / 20) * 12, face = "bold"),
        legend.title = element_blank(),
        legend.position = legend.position,
        legend.text = element_text(size = (plot.size / 20) * 12, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = (plot.size / 20) * 12, face = "bold")
      )

    # Boxplot
    dodge <- position_dodge(width = 0.8)
    p2 <- ggplot(d2,
                 aes(
                   x = Group,
                   y = type,
                   fill = Group,
                   group = Group
                 )) +
      stat_boxplot(
        geom = 'errorbar',
        width = 0.3,
        size = 1.5,
        position = dodge
      ) +
      geom_boxplot(
        outlier.shape = NA,
        width = 0.6,
        size = 1.5,
        color = 'black',
        position = dodge
      ) +
      geom_point(
        position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
        aes(group = Group),
        alpha = 0.5,
        size = 4
      ) +
      scale_fill_manual(breaks = names(list.color),
                        values = as.character(list.color)) +
      labs(x = NULL, y = type) +
      facet_wrap(. ~ Time, scales = 'fixed') +
      theme_bw() +
      theme(
        axis.text = element_text(
          size = plot.size / 15 * 12,
          colour = "black",
          face = "bold"
        ),
        axis.title.x = element_text(
          size = plot.size,
          colour = "black",
          face = "bold"
        ),
        axis.title.y = element_text(
          size = plot.size,
          colour = "black",
          face = "bold"
        ),
        legend.text = element_text(
          size = plot.size / 15 * 12,
          colour = "black",
          face = "bold"
        ),
        legend.title = element_text(
          size = plot.size / 15 * 12,
          colour = "black",
          face = "bold"
        ),
        legend.position = legend.position,
        strip.text = element_text(
          size = plot.size / 15 * 12,
          colour = "black",
          face = "bold"
        )
      )

  } else {
    if (verbose)
      message('Please enter right type:  ',
              paste0(all_type, collapse = ', '),
              '.')
  }


  # Output results
  l <- list(lineGraph = p1,
            boxPlot = p2)
  if (verbose){
    print(plot_grid(
      l$lineGraph,
      l$boxPlot,
      ncol = 1,
      labels = LETTERS[1:2],
      label_size = plot.size / 25 * 20,
      rel_heights = c(2,1)
    ))
    message('Done!')
  }
  return(l)
}

setGeneric("statistics", function(object, ...) {
  standardGeneric("statistics")
})

#' @rdname Animal-method
#' @name Animal-method
#' @description \code{statistics} method for \code{Animal} class
#' @param targetTime Strings like '20220120'. Default is \code{NULL}, which means that every time would be considerated.
#' @inheritParams ggpubr::compare_means
#' @importFrom dplyr filter
#' @importFrom ggpubr compare_means
#' @importFrom DT datatable
#' @importFrom rstatix t_test
#' @importFrom rstatix wilcox_test
#' @exportMethod statistics
setMethod(
  "statistics",
  signature(object='Animal'),
  function(object,
           type = c('tumorVolume', 'weight')[1],
           targetTime = NULL,
           method = "wilcox.test",
           paired = FALSE,
           p.adjust.method = 'BH',
           verbose = T){

    # Data
    d <- as.data.frame(object)
    if(!is.null(targetTime)){
      d <- filter(d, Time %in% targetTime)
    }
    d2 <- d[c("Time", "Group", "ID", type)]; colnames(d2)[ncol(d2)] <- 'type'

    # Comparison of Means
    # res <- compare_means(type ~ Group,
    #                      data = d2,
    #                      method = "wilcox.test",
    #                      group.by = 'Time',
    #                      paired = FALSE,
    #                      p.adjust.method = 'BH')


    # P value function
    fun_comparision <- function(method){
      if(method == "wilcox.test"){
        if(verbose) message('Use Wilcox test.')
        rstatix::wilcox_test
      } else if(method == "t.test"){
        if(verbose) message('Use t test.')
        rstatix::t_test
      } else {
        if(verbose) message('Wrong p value method. Use t test for default.')
        rstatix::t_test
      }
    }
    fun_comparision_target <- fun_comparision(method)


    # P value
    time <- unique(d2$Time)
    df.p <- data.frame()
    for(i in 1:length(time)){ # i=1
      d2.i <- d2[d2$Time %in% time[i],]
      df.p_i <- fun_comparision_target(
        data = d2.i,
        formula = type ~ Group,
        p.adjust.method = p.adjust.method,
        paired = paired,
        alternative = "two.sided",
        mu = 0,
        conf.level = 0.95,
        detailed = FALSE
      )
      df.p_i2 <- cbind(Time = time[i], df.p_i)
      df.p <- rbind(df.p, df.p_i2)
    }

    if(verbose) print(DT::datatable(df.p))

    # Output
    return(df.p)
  }
)



