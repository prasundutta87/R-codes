library(readr)
library(stringr) #for str_sort
library(R.utils) #for listDirectory
library(data.table) ##for rbindlist
library(qqman)
options(scipen = 999)

#"setwd("~/Thesis_PhD/Graphs_tables_for_thesis/xpehh/")Banni_vs_Med/")
#setwd("~/Thesis_PhD/Graphs_tables_for_thesis/xpehh/")
#tiff("all.tiff", units="in", width=25, height=11, res=100)
#par(mfrow=c(2,2))
dir = listDirectory(pattern = "*_vs_*")
xpehh_tables = list()

for (i in 1:length(dir)) {
  filenames = str_sort(list.files(dir[i], pattern = "*_vs_*"), numeric = T)
  for (j in 1:length(filenames)) {
    dt = read_delim(
      paste0(dir[i], "/", filenames[j]),
      "\t",
      escape_double = FALSE,
      col_types = cols(
        iHH_A1 = col_skip(),
        iHH_B1 = col_skip(),
        iHH_P1 = col_skip()
      ),
      trim_ws = TRUE
    )
    dt$chrom = j
    dt$SNP_ID = paste(dt$chrom, dt$Location, sep = "_")
    dt$z_score = scale(dt$XPEHH, center = T, scale = T)
    dt$z_score_abs = abs(dt$z_score)
    dt$p_value = (1 - pnorm(dt$z_score_abs)) * 2
    #dt = dt[, c(3,1,4,2,5,6)]
    xpehh_tables[[j]] = dt
  }
  dt_final = rbindlist(xpehh_tables)
  png(
    paste0(dir[i], "_QQ_plot", ".png"),
    units = "in",
    width = 15,
    height = 11,
    res = 600
  )
  qq(dt_final$p_value, main = dir[i], col = "blue4")
  dev.off()
}

