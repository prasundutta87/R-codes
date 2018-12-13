##This script takes in SNP IDs from heatmap clusters, takes their corresponding clumped SNPs from "clump_output.clumped"
##and converts the resulting intermediate file to 0-based coordinate bed format producing bed file for each cluster

library(stringr) ##for str_sort
library(tidyr) #for separate_rows
library(data.table) ##for rbindlist,fread,fwrite
library(dplyr) ##for left_join
options(stringsAsFactors = F)

setwd("~/Thesis_PhD/Graphs_tables_for_thesis/xpehh/")
filenames=str_sort(list.files(pattern = "cluster_\\d_SNPs.txt"),numeric = T)
clumped=fread("clump_output.clumped")

for(i in 1:7){
  output_filename=strsplit(filenames[i],".txt")[[1]][1]
  cluster=fread(filenames[i])
  dt=left_join(cluster,select(clumped,SNP,SP2),by="SNP")
  dt$SP2[dt$SP2=="NONE"]=" " ##make second column empty wherever first column has the word 'NONE'
  dt$SP2=paste0(dt$SP2,",",dt$SNP) ##add data from first column to the second column in every row
  dt=separate_rows(dt,SP2,sep = ",") 
  dt = dt[-which(dt$SP2 == " "), ] ##remove rows where second column is empty as they got generated because of an additional comma that got added in the empty cells from where 'NONE' was removed
  dt$SNP=as.numeric(sapply(dt$SNP,function(SNP) strsplit(SNP,"_")[[1]][2]))
  dt$SP2_new=as.numeric(sapply(dt$SP2,function(SP2) strsplit(SP2,"_|\\(")[[1]][2]))
  dt$chrom=as.numeric(sapply(dt$SP2,function(SP2) strsplit(SP2,"_|\\(")[[1]][1]))
  dt$SP2=NULL
  dt$chromStart=pmin(dt$SNP,dt$SP2_new)
  dt$chromEnd=pmax(dt$SNP,dt$SP2_new)
  dt$chromStart=dt$chromStart-1 ##converting 1-based corrdinate to 0-based corrdinate as per bed format
  dt=dt[, -c(1,2)] ##delete first two columns named 'SNP' and SP2_new
  fwrite(dt,paste0(output_filename,".bed"),row.names = F,col.names = F,quote = F,sep="\t")
}