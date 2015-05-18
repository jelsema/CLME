

#' @title 
#' Fibroid Growth Study
#' 
#' @description
#' This data set contains a subset of the data from the Fibroid Growth Study.
#'
#' \tabular{rll}{
#' [,1] \tab  ID  \tab ID for subject. \cr
#' [,2] \tab fid  \tab ID for fibroid (each women could have multiple fibroids). \cr
#' [,3] \tab lfgr \tab log fibroid growth rate. See details. \cr
#' [,4] \tab age  \tab age category Younger, Middle, Older. \cr
#' [,5] \tab loc  \tab location of fibroid, corpus, fundus, or lower segment. \cr
#' [,6] \tab bmi  \tab body mass index of subject. \cr
#' [,7] \tab preg \tab parity, whether the subject had delivered a child. \cr
#' [,8] \tab race \tab race of subject (Black or White only). \cr
#' [,9] \tab vol  \tab initial volume of fibroid. \cr
#' }
#' 
#' @details
#' The response variable \code{lfgr} was calculated as the change in log fibroid volume, 
#' divided by the length of time between measurements. The growth rates were averaged to produce
#'  a single value for each fibroid, which was scaled to represent a 6-month percent change in volume.
#' 
#' @references
#' Peddada, Laughlin, Miner, Guyon, Haneke, Vahdat, Semelka, Kowalik, Armao, Davis, and Baird(2008).
#'  Growth of Uterine Leiomyomata Among Premenopausal Black and White Women.
#'   Proceedings of the National Academy of Sciences of the United States of America, 105(50),
#'    19887-19892. URL \url{http://www.pnas.org/content/105/50/19887.full.pdf}.
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name fibroid
#' @usage data(fibroid)
#' @format A frame containing 240 observations on 9 variables.
NULL



#' @title 
#' Experiment on mice
#' 
#' @description
#' This data set contains the data from an experiment on 24 Sprague-Dawley rats from Cora et al (2012).
#'
#' \tabular{rll}{
#' [,1]  \tab id   \tab ID for rat (factor). \cr
#' [,2]  \tab time \tab time period (in order, 0 , 6, 24, 48, 72, 96 hours). \cr
#' [,3]  \tab temp \tab storage temperature reference (\code{''Ref''}) vs. room temperature (\code{''RT''}). \cr
#' [,4]  \tab sex  \tab sex, male (\code{''Male''}) vs. female (\code{''Female''}). Coded as \code{''Female''=1}. \cr
#' [,5]  \tab wbc  \tab white blood cell count (\eqn{10^3 / \mu L}{10^3 / mu L}). \cr
#' [,6]  \tab rbc  \tab red blood cell count )\eqn{10^6 / \mu L}{10^6 / mu L}). \cr
#' [,7]  \tab hgb  \tab hemoglobin concentration (g/dl). \cr
#' [,8]  \tab hct  \tab hematocrit (\%). \cr
#' [,9]  \tab spun \tab (HCT \%). \cr
#' [,10] \tab mcv  \tab MCV, a measurement of erythrocyte volume (fl). \cr
#' [,11] \tab mch  \tab mean corpuscular hemoglobin (pg). \cr     % ????
#' [,12] \tab mchc \tab mean corpuscular hemoglobin concentration (g/dl). \cr
#' [,13] \tab plts \tab platelet count (\eqn{10^3 / \mu L}{10^3 / mu L}). \cr
#' }
#' 
#' @details
#' The response variable \code{lfgr} was calculated as the change in log fibroid volume, 
#' divided by the length of time between measurements. The growth rates were averaged to produce
#'  a single value for each fibroid, which was scaled to represent a 6-month percent change in volume.
#' 
#' @references
#' Cora M, King D, Betz L, Wilson R, and Travlos G (2012). 
#' Artifactual changes in Sprauge-Dawley rat hematologic parameters after storage of samples at 3 C and 21 C.
#'  Journal of the American Association for Laboratory Animal Science, 51(5), 616-621. 
#'  URL \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3447451/}.
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name rat.blood
#' @usage data(rat.blood)
#' @format A frame containing 241 observations on 13 variables.
NULL














