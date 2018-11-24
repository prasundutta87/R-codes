##This script converts a VCF file produced from plink using '--recode vcf' into a form suitable for Graphia with added metadata##
##Usage: Rscript script_name.R <VCF_file> <P-value_file_in_xls_or_xlsx_format> <metadata_file_in_xls_or_xlsx_format>##
library(readr) ##for read_delim
library(readxl) ##for read_excel
library(plyr) ##for join

args = commandArgs(trailingOnly = TRUE) ##accept filename from command prompt
system(paste("grep", " '^##' -v ", args[1], " > temporary_file.vcf")) ## make a temporary file without VCF header
graphia_input = read_delim("temporary_file.vcf",
                           "\t",
                           trim_ws = TRUE,
                           col_names = T) ##Read input file as data frame
system("rm temporary_file.vcf") ## remove the temporary file
p_value = read_excel(args[2], na = "NA", trim_ws = T) ##Read p value input file as data frame
graphia_input = join(graphia_input, p_value, by = 'ID') ##add p values to graphia_input data frame, 3 columns added to graphia_input data frame
graphia_input = graphia_input[, c(2:1, 3:9, 952:954, 10:951)] ##reorder data frame, 951 being total columns in graphia_input data frame

graphia_input[graphia_input == '0/0'] = 0 ##change 0/0 to 0
graphia_input[graphia_input == '0/1'] = 1 ##change 0/1 to 1
graphia_input[graphia_input == '1/1'] = 2 ##change 1/1 to 2

##Loop to change all column names from 0_<sample_number> to only <sample_number>
for (i in 13:ncol(graphia_input))
{
  colnames(graphia_input)[i] = strsplit(colnames(graphia_input)[i], "_")[[1]][2]
}

graphia_input[, 1:12] <-
  lapply(graphia_input[, 1:12], as.character) ##change graphia_input data frame columns to characters for
##easy binding to metadata data frame

################add metadata rows to graphia input file################
metadata = read_excel(args[3],
                      sheet = "USe me",
                      na = "NA",
                      trim_ws = T) ##read excel file containing metadata as data frame
temp = data.frame(colnames(graphia_input[13:ncol(graphia_input)])) ##create dataframe of sample names
colnames(temp)[1] = 'bird_id' ## change column name to 'bird_id' to match metadata excel file
metadata_new = join(temp, metadata, by = 'bird_id')
metadata_t = as.data.frame(t(metadata_new[, -1])) ##transpose data frame
colnames(metadata_t) = metadata_new$bird_id ##update columnnames to bird IDs
new_metadata_t = cbind(
  "POS" = "",
  "#CHROM" = "",
  "ID" = "",
  "REF" = "",
  "ALT" = "",
  "QUAL" = "",
  "FILTER" = "",
  "INFO" = "",
  "FORMAT" = "",
  "WG_p" = "",
  "CLS_p" = "",
  "IL10_p" = "",
  metadata_t
) ##add extra columns to match column numbers/names of graphia_input data frame
new_metadata_t[, 1] = rownames(new_metadata_t) ##make rownames as data for first column
rownames(new_metadata_t) = c() ##remove row names
##################################################################

graphia_input_with_metadata = rbind(new_metadata_t, graphia_input) ##Bind the two data frames
write.csv(
  graphia_input_with_metadata,
  "graphia_input_processed_for_Tom_with_metadata.csv",
  row.names = F,
  quote = F
) ##write file to disk
