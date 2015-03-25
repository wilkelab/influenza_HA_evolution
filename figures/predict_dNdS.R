rm(list = ls())
library(ggplot2)
library(grid)
library(cowplot)
library(igraph)

setwd('~/Google Drive/Data/influenza_HA_evolution/')
df <- read.table('manuscript/numbering_table.csv', head=T, sep=',', stringsAsFactors = F)