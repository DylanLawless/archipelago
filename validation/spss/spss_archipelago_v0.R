# library(dplyr)

# Check data funciton ----
check_archipelago_input <- function(df1, df2) {
  # Required columns
  required_df1 <- c("set_ID", "P")
  required_df2 <- c("set_ID", "P", "CHR", "BP")
  
  # Check required columns
  missing_df1 <- setdiff(required_df1, names(df1))
  missing_df2 <- setdiff(required_df2, names(df2))
  if (length(missing_df1) > 0) stop("df1 is missing: ", paste(missing_df1, collapse = ", "))
  if (length(missing_df2) > 0) stop("df2 is missing: ", paste(missing_df2, collapse = ", "))
  
  # Check for NAs
  na_df1 <- required_df1[sapply(df1[required_df1], anyNA)]
  na_df2 <- required_df2[sapply(df2[required_df2], anyNA)]
  
  if (length(na_df1) > 0) {
    msg <- paste("❌ [df1 contains NA in:", paste(na_df1, collapse = ", "))
    if ("set_ID" %in% na_df1) msg <- paste0(msg, "\nSuggestion: df1$set_ID[is.na(df1$set_ID)] <- 0")
    stop(msg)
  }
  
  if (length(na_df2) > 0) {
    msg <- paste("❌ [df2 contains NA in:", paste(na_df2, collapse = ", "))
    if ("set_ID" %in% na_df2) msg <- paste0(msg, "\nSuggestion: df2$set_ID[is.na(df2$set_ID)] <- 0")
    stop(msg)
  }
  
  # Type checks
  if (!is.numeric(df2$CHR) || !is.numeric(df2$BP)) stop("df2$CHR and df2$BP must be numeric")
  if (!is.numeric(df1$P) || !is.numeric(df2$P)) stop("P columns must be numeric")
  
  # Transformation warnings
  if (max(df1$P, na.rm = TRUE) > 1.5) {
    message("❌ [warning] df1$P appears to be -log10 transformed. To revert: df1$P <- 10^(-df1$P)")
  }
  if (any(df1$P <= 0, na.rm = TRUE)) {
    message("❌ [warning]  df1 contains P values <= 0. Consider: df1 <- dplyr::filter(df1, P > 0)")
  }
  
  # ID consistency checks
  missing_in_df1 <- setdiff(unique(df2$set_ID), unique(df1$set_ID))
  missing_in_df2 <- setdiff(unique(df1$set_ID), unique(df2$set_ID))
  if (length(missing_in_df2) > 0) {
    message("❌ [warning] Some set_IDs in df1 have no matching rows in df2.")
    message("          Suggested: df1 <- dplyr::filter(df1, set_ID %in% df2$set_ID)")
  }
  if (length(missing_in_df1) > 0) {
    message("[info] df2 contains set_IDs not in df1 (ignored during plotting).")
  }
  
  message("✅ Input checks passed")
}

# real data ----
# Define output path
# data_dir <- "../data/"
data_dir <- "~/mnt/atlas_data_big/data/archipelago_usage/"

df1 <- read.table(file = file.path(data_dir, "spss/df1_SPSS_archipelago.tsv"), sep = "\t", header = TRUE)

df2 <- read.table(file = file.path(data_dir, "spss/df2_SPSS_archipelago.tsv"), sep = "\t", header = TRUE)

df1 <- df1 |> dplyr::select(set_ID, P)
df2 <- df2 |> dplyr::select(set_ID, BP, P, CHR, SNP)

check_archipelago_input(df1, df2)
df2$set_ID[is.na(df2$set_ID)] <- 0
df1 <- df1 |> dplyr::mutate(P = 10^(-P))
df1 <- df1 |> dplyr::filter(P > 0) # any Inf
df1 <- dplyr::filter(df1, set_ID %in% df2$set_ID)
check_archipelago_input(df1, df2)

# Install and run ----
# remove.packages("archipelago")
# detach("package:archipelago", unload = TRUE)
# .rs.restartR() 
# install.packages("/Users/dylanlawless/web/archipelago_0.1.0.tar.gz", repos = NULL, type = "source")
# "archipelago" %in% rownames(installed.packages())
# find.package("archipelago")
library(archipelago)
# archipelago_plot(df1, df2) 

plot_subtitle <- paste0("Dataset: ", nrow(df1), " VSAT and ", format(nrow(df2), big.mark = ","), " SNPs")
plot_subtitle

plot_title <- "SPSS Trait 1"
plot_subtitle <- plot_subtitle
output_path = paste0(data_dir, "output/archipelago_plot_spss")
output_raw = paste0(data_dir, "output/archipelago_plot_vsat_raw_spss")

archipelago_plot(df1 = df1,
                 d = df2,
                 add_title = TRUE,
                 plot_title = plot_title,
                 add_subtitle = TRUE,
                 plot_subtitle = plot_subtitle,
                 color_theme = 'retro',
                 output_path = output_path,
                 output_raw = output_raw
)
