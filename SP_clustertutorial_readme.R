library(reshape2); 
library(vegan); 
library(dplyr)
library(ade4); 
library(plotly)
library(compositions); 
library(pracma); 
library(DESeq2); 
library(fpc); 
library(tidyverse)
library(purrr)
library(cluster)
library(RColorBrewer)
library(ape)
setwd("/Users/taylorsehein/github/WHOI_repo")
raw <- read.delim("18SiTag_SP_forCluster.txt")
head(raw)
colnames(raw)
str(raw)
summary(raw)
count.filtered <- raw
names(count.filtered)
seq_counts<-count.filtered[1:341]
tax_key<-count.filtered[c(1,341)]; head(tax_key[1:2,])