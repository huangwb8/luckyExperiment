% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clinic_getDCA.R
\name{getDCA}
\alias{getDCA}
\title{Fast way to draw net benefit curves}
\usage{
getDCA(
  data,
  group.by = NULL,
  time.col = "time",
  status.col = "event",
  clusters,
  time.raw.type = c("Day", "Month", "Year")[1],
  time.target.type = c("Day", "Month", "Year")[3],
  time.knot = c(3, 5),
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = 0.95,
  study.design = "cohort",
  population.prevalence = c(NA, NA),
  cost.benefit.axis = FALSE,
  col = grDevices::rainbow(n = 8, v = 0.8),
  plot.confidence.intervals = FALSE,
  standardize = FALSE,
  plot.type = c("time", "cluster", "disperse")[1],
  plot.title.position = list(`3` = c(0.85, 0.24), `5` = c(0.85, 0.24)),
  plot.label,
  names = "project"
)
}
\arguments{
\item{data}{a data frame with survival data(time & event)}

\item{group.by}{group colname}

\item{time.col}{colname of time value}

\item{status.col}{Numeric.\code{status.col} must be a binary parameter with 0(event not happened yet) or 1(event happened) value.}

\item{clusters}{list of characters.Some of the colnames of the data representing the value involving the glm model.}

\item{time.raw.type}{one of "Day","Month","Year"}

\item{time.target.type}{one of "Day","Month","Year"}

\item{time.knot}{time knot to create survival probability}

\item{family}{a description of the error distribution and link function to pass to 'glm' used for model fitting. Defaults to binomial(link = 'logit') for logistic regression.}

\item{thresholds}{Numeric vector of high risk thresholds to use when plotting and calculating net benefit values.}

\item{confidence.intervals}{Numeric (default 0.95 for 95\% confidence bands) level of bootstrap confidence intervals to plot. Set as NA or 'none' to remove confidence intervals. See details for more information.}

\item{study.design}{Either 'cohort' (default) or 'case-control' describing the study design used to obtain data. See details for more information.}

\item{population.prevalence}{Outcome prevalence rate in the population used to calculate decision curves when study.design = 'case-control'.}

\item{cost.benefit.axis}{logical (default TRUE) indicating whether to print an additional x-axis showing relative cost:benefit ratios in addition to risk thresholds.}

\item{col}{vector of color names to be used in plotting corresponding to the 'predictors' given. Default colors will be chosen from rainbow(..., v = .8). See details for more information on plot parameters.}

\item{plot.confidence.intervals}{whether plot confidence intervals}

\item{standardize}{logical (default TRUE) indicating whether to use the standardized net benefit (NB/disease prevalence) or not.}

\item{plot.type}{Character. One/some of c("time","cluster","disperse"). if \code{plot.type="time"},then multiple clusters would be ploted for every time; if \code{plot.type="cluster"}, then multiple times would be ploted for every cluster; if \code{plot.type="disperse"}, then every single plot would be ploted indivitually.}

\item{plot.title.position}{List. The title position of the plot.Only available when \code{plot.type="time"} or \code{plot.type="cluster"}. Its length must be the same of \code{time.knot}}

\item{plot.label}{List. The name of list are the paste0 of elements in one of \code{clusters} with "-".The content of the list are the real names that would be printed in the plot legend/title.}

\item{names}{part of save files or plots.}
}
\value{
DCA plots and related model list
}
\description{
Fast way to draw net benefit curves
}
\details{
The length of \code{population.prevalence} should be equal to \code{time.knot}. Default of \code{population.prevalence} is \code{c(NA,NA)}, \cr which means we use an estimated population prevalence based on raw dataset.
}
\examples{
## data preparation
data("myeloma",package = "survminer");
data <- myeloma;rm(myeloma)
data$group <- ifelse(data$molecular_group \%in\% c( "Cyclin D-1", "Cyclin D-2", "Proliferation"),"proliferation","others")

## Here is a useful example:

# group by
a <- getDCA(data,
            group.by = "group",
            time.col = "time",
            status.col = "event",
            clusters = list(c("TP53","WHSC1"),
                            c("TP53","WHSC1","DEPDC1")),
            time.raw.type=c("Day","Month","Year")[2],
            time.target.type=c("Day","Month","Year")[3],
            time.knot = c(3,5),
            plot.type = c("time","cluster","disperse"),
            plot.title.position = list("3" = c(0.85,0.24),
                                       "5" = c(0.85,0.24)),
            plot.label = list("TP53-WHSC1" = "simpleModel",
                  "TP53-WHSC1-DEPDC1" = "complexModel"),
            names="project")

# none group
a <- getDCA(data,
            group.by = NULL,
            time.col = "time",
            status.col = "event",
            clusters = list(c("TP53","WHSC1"),
                            c("TP53","WHSC1","DEPDC1")),
            time.raw.type=c("Day","Month","Year")[2],
            time.target.type=c("Day","Month","Year")[3],
            time.knot = c(3,5),
            plot.type = c("time","cluster","disperse"),
            plot.title.position = list("3" = c(0.85,0.24),
                                       "5" = c(0.85,0.24)),
            plot.label = list("TP53-WHSC1" = "simpleModel",
                  "TP53-WHSC1-DEPDC1" = "complexModel"),
            names="project")

## Of cause, you can select a sigle plot type like "time", which
## is the most comman calling:
a <- getDCA(data,
            group.by = "group",
            time.col = "time",
            status.col = "event",
            clusters = list(c("TP53","WHSC1"),
                            c("TP53","WHSC1","DEPDC1")),
            time.raw.type=c("Day","Month","Year")[2],
            time.target.type=c("Day","Month","Year")[3],
            time.knot = c(3,5),
            plot.type = c("time","cluster","disperse")[1],
            plot.title.position = list("3" = c(0.85,0.24),
                                       "5" = c(0.85,0.24)),
            plot.label = list("TP53-WHSC1" = "simpleModel",
                  "TP53-WHSC1-DEPDC1" = "complexModel"),
            names="project")
}
\seealso{
\code{\link[rmda]{plot_decision_curve}}.
}
\author{
Weibin Huang<\email{654751191@qq.com}>
}
