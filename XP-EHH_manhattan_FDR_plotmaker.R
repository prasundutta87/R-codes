library(readr)
library(stringr) #for str_sort
library(R.utils) #for listDirectory
library(data.table) ##for rbindlist
library(qqman)
options(scipen = 999)

setwd("~/Thesis_PhD/Graphs_tables_for_thesis/xpehh/")

dir = listDirectory(pattern = "*_vs_*")
xpehh_tables = list()

for (i in 1:length(dir)) {
  filenames = str_sort(list.files(dir[i],pattern = "*_vs_*"), numeric = T)
  for (j in 1:length(filenames)) {
    dt = read_delim(
      paste0(dir[i],"/",filenames[j]),
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
    dt$FDR_adj_p_value=p.adjust(dt$p_value,method="BH",n=length(dt$p_value))
    xpehh_tables[[j]] = dt
  }
  dt_final= rbindlist(xpehh_tables)
  temp_dt=subset(dt_final,dt_final$FDR_adj_p_value <=0.05)
  z=temp_dt$z_score[which.max(temp_dt$FDR_adj_p_value)]
  png(paste0(dir[i],"_FDR_line",".png"), units="in", width=25, height=11, res=600)
  manhattan(
    dt_final,
    chr = "chrom",
    bp = "Location",
    p = "z_score",
    snp = "SNP_ID",
    logp = FALSE,
    ylim = c(min(dt_final$z_score) - 0.5, max(dt_final$z_score) + 0.5),
    ylab = "XPEHH Z-scores",
    genomewideline = z,
    suggestiveline = -z,
    main = dir[i]
  )
  dev.off()
}
