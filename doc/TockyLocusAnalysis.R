## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----tockyprep, include=TRUE--------------------------------------------------
library(TockyPrep)
library(TockyLocus)

## ----files, include=TRUE------------------------------------------------------
# Example data load
# Define the base path
file_path <- system.file("extdata", package = "TockyLocus")

# Define files
negfile <- "Timer_negative.csv"
samplefiles <- list.files(file_path, pattern = "sample_", full.names = FALSE)
samplefiles <- setdiff(samplefiles, file.path(file_path, negfile))

## ----prep, include=TRUE-------------------------------------------------------
# Preprocessing data
prep <- prep_tocky(path = file_path, samplefile = samplefiles, negfile = negfile, interactive = FALSE)

## ----timer_transform, include=TRUE--------------------------------------------
# Normalizing and transforming data
x <- timer_transform(prep, blue_channel = 'Timer.Blue', red_channel = 'Timer.Red', select = FALSE, verbose = FALSE)

## ----class, include=TRUE------------------------------------------------------
class(x)

## ----sample_definition1, include = TRUE---------------------------------------
sample_definition <- read.csv(file.path(file_path, 'sampledef.csv'))
sample_definition <- as.data.frame(sample_definition)
head(sample_definition)
# Normalizing and transforming data

x <- sample_definition(x, sample_definition = sample_definition, interactive = FALSE)


## ----sample_definition_alternative, include=FALSE-----------------------------
#x <- sample_definition(x, output_dir = 'outuput', interactive = TRUE)


## ----plotAngleDensity, fig.width=4, fig.height=4------------------------------
# Visualizing the results
plotAngleDensity(x)

## ----TockyLocus, include=TRUE-------------------------------------------------
x <- TockyLocus(x)


## ----PlotTockyLocusLegend, fig.width=5, fig.height=2, include=TRUE------------
plotTockyLocusLegend(mar_par = c(2, 2, 4, 2))


## ----plot_tocky_locus, fig.width=8, fig.height=8, include=TRUE----------------
plot_tocky_locus(x, n = 4)


## ----PlotTockyLocus, fig.width=8, fig.height=4, include=TRUE------------------
plotTockyLocus(x, verbose = FALSE)


## ----PlotTockyLocus_Locus, fig.width=6, fig.height=4, include=TRUE------------
plotTockyLocus(x, group_by = FALSE,  verbose = FALSE)


