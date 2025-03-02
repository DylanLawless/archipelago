library(dplyr)
library(data.table)
library(stringr)

path_mart <- "~/web/archipelago/validation/07_SKATO_pathway/mart_export_GRCh37p13.txt"
path_plink <- "~/web/archipelago_validation/04_Data_QC/sample_data.clean.bim"

process_bim_file <- function(file_name) {
  print("read bim")
  data_bim <- fread(file_name, header = FALSE, stringsAsFactors = FALSE)
  setnames(data_bim, c("chromosome", "SNPID", "misc", "position", "allele1", "allele2"))
  print(head(data_bim))
  
  print("Reorder allele columns")
  data_bim[, SNPID := paste(chromosome, position, allele2, allele1, sep = ":")]
  print(head(data_bim))
  
  print("read mart")
  data_mart <- fread(path_mart, header = TRUE, sep = "\t")
  setnames(data_mart, c("Gene_stable_ID", "Gene_start", "Gene_end", "chromosome", "Gene_name"))
  data_mart <- unique(data_mart[, .(Gene_stable_ID, Gene_start, Gene_end, chromosome, Gene_name)])
  print(head(data_mart))
  
  data_bim[, chromosome := as.character(chromosome)]
  data_mart[, chromosome := as.character(chromosome)]
  data_bim[, position := as.numeric(as.character(position))]
  print(head(data_mart))
  
  result <- data_mart[data_bim,
                      on = .(chromosome = chromosome,
                             Gene_start <= position,
                             Gene_end >= position),
                      .(chromosome, SNPID = i.SNPID, misc = i.misc, position = i.position, allele1 = i.allele1, allele2 = i.allele2, Gene_name = Gene_name)]
  
  result <- unique(result)
  print("head(result)")
  print(head(result))
  
  out_file_full <- "~/web/archipelago_validation/04_Data_QC/sample_data.clean.full.tsv"
  fwrite(result, out_file_full, sep = "\t")
  
  result_slim <- result %>% dplyr::select(Gene_name, SNPID)
  result_slim <- result_slim[result_slim$Gene_name != "" & !is.na(result_slim$Gene_name), ]
  out_file_slim <- "~/web/archipelago_validation/04_Data_QC/sample_data.clean.bim.tsv"
  fwrite(result_slim, out_file_slim, sep = "\t", col.names = FALSE)
  
  return(list(full = result, slim = result_slim))
}

results <- process_bim_file(path_plink)
head(results$full)
head(results$slim)
