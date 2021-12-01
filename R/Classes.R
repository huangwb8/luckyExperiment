

##' Class "LuckyCCK8"
##' @name LuckyCCK8-class
##' @docType class
##' @slot Repeat repeat parameters
##' @slot Data datasets
##' @slot Plot ggplot
##' @author Weibin Huang<\email{654751191@@qq.com}>
##' @keywords classes
##' @exportClass LuckyCCK8
setClass(
  Class = "LuckyCCK8",
  slots = c(
    Repeat = "list",
    Data = "list",
    Plot = "list"
  )
)

if(F){
  prototype = list(
    Repeat = list(
      cell.type = "character",
      group = "character",
      empty.group = "character",
      stimulate.time = "character",
      observe.time = "character",
      value = "character"
    ),
    Data = list(
      metaData = "data.frame",
      adjData = "data.frame",
      plotData = "data.frame"
    )
  )
}










