# ------------------------------------------------------------------------------
# UK Biobank validation of archipelago using GWAS and DeepRVAT gene-level tests
# ------------------------------------------------------------------------------

# This script demonstrates how the archipelago method integrates SNP-level GWAS 
# and gene-level rare variant association test (RVAT) results using data from 
# the UK Biobank. We focus on the trait "platelet distribution width" (phenocode 30110).

# Summary:
# - GWAS summary statistics (29 million SNPs) were obtained from the Pan-UK Biobank 
#   project (Broad Institute). We downsampled every 200th row (~145k SNPs) for plotting.
# - Gene-level p-values were sourced from DeepRVAT, a method that integrates rare 
#   variant annotations using deep set neural networks to improve power in RVAT.
# - DeepRVAT tested 3.49 million gene–trait pairs across 97 traits in 470,000 WES samples.
# - For phenocode 30110, 35,968 gene-level tests were available.
# - We used Ensembl GRCh38 gene coordinates to link SNPs to genes and mapped unmatched 
#   genes to their nearest SNP proxy using genomic midpoints.
# - The resulting SNP and gene-level p-values were plotted with archipelago, 
#   highlighting overlap and signal proximity between rare and common variants.

# Observation:
# - Only one gene (DOCK5) showed significant evidence in both layers.
# - Archipelago identified 26 additional RVAT-significant genes with nearby GWAS signal.

# Source studies and references:
# - Pan-UK Biobank GWAS manuscript: 
#     Pan-UK Biobank GWAS improves discovery, analysis of genetic architecture, and resolution into ancestry-enriched effects. Konrad J. Karczewski, Rahul Gupta, Masahiro Kanai, Wenhan Lu, Kristin Tsuo, Ying Wang, Raymond K. Walters, Patrick Turley, Shawneequa Callier, Nirav N. Shah, Nikolas Baya, Duncan S. Palmer, Jacqueline I. Goldstein, Gopal Sarma, Matthew Solomonson, Nathan Cheng, Sam Bryant, Claire Churchhouse, Caroline M. Cusick, Timothy Poterba, John Compitello, Daniel King, Wei Zhou, Cotton Seed, Hilary K. Finucane, Mark J. Daly, Benjamin M. Neale, Elizabeth G. Atkinson, Alicia R. Martin. medRxiv 2024.03.13.24303864; doi: https://doi.org/10.1101/2024.03.13.24303864
# - Pan-UK Biobank GWAS dataset:
#     Pan-UK Biobank GWAS summary statistics: https://pan.ukbb.broadinstitute.org. GWAS file (continuous-30110-both_sexes-irnt.tsv.bgz) (size 2.1G)
# 
# - DeepRVAT manuscript:
#     Clarke, B., Holtkamp, E., Öztürk, H. et al. 
#     Integration of variant annotations using deep set networks boosts rare variant association testing.
#     Nat Genet 56, 2271–2280 (2024). https://doi.org/10.1038/s41588-024-01919-z
# - DeepRVAT dataset:
#     Clarke, B., & Holtkamp, E. (2024). DeepRVAT gene-trait association testing results on the 470k UK Biobank WES dataset [Data set]. Zenodo. https://doi.org/10.5281/zenodo.12736824  (UKBB_470k_deeprvat_results.csv) (size  387.7 MB)

# This script downloads DeepRVAT gene-level summary statistics and selects a specific 
# phenotype ("platelet distribution width", phenocode 30110) matched to the Pan-UK Biobank study manifest 
# (<https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz>).
# It then retrieves the corresponding GWAS summary statistics 
# (<https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/continuous-30110-both_sexes-irnt.tsv.bgz>),
# and prepares both datasets for integrated visualisation with Archipelago.

# Description of DeepRVAT:
# DeepRVAT is a rare variant association framework that uses deep set neural networks 
# to learn gene impairment scores from variant annotations and phenotypes. It improves 
# gene discovery and trait prediction by enabling scalable, well-calibrated RVAT at 
# biobank scale. In the UK Biobank, DeepRVAT identified more associated genes and 
# enhanced detection of high-risk individuals compared to conventional RVAT.

# ------------------------------------------------------------------------------

# TODO: Add downsamplign as package feature for massive studies. It is included manually below. Our final study method includes the following note: "We downsampled the 35,968 gene-level RVAT results by selecting the top 0.1\% most significant associations and sampling an additional 1\% from the remaining results using a ramped log-scale weighting based on $-\log_{10}(p)$. This prioritised moderately significant genes while reducing visual clutter from uniformly null signals, yielding a focused subset of 396 genes for plotting."

# Install and run ----
# remove.packages("archipelago")
# detach("package:archipelago", unload = TRUE)
# .rs.restartR()
# install.packages("/Users/dylanlawless/web/archipelago_0.1.0.tar.gz", repos = NULL, type = "source")
# "archipelago" %in% rownames(installed.packages())
# find.package("archipelago")
# library(archipelago)
# archipelago_plot(df1, df2, annotate_thresholds = FALSE)

# Load required packages
library(ggplot2)
library(readr)
library(dplyr)
library(data.table)
library(R.utils)
library(biomaRt)
library(readr)

print("Large data stored on local NAS")
print("~/so/src/data_big/mount_atlas.sh")

data_dir <- "~/mnt/atlas_data_big/data/archipelago_usage/"
# data_dir <- "~/Desktop/"

# Results from DeepRVAT ----
deeprvat_url <- "https://zenodo.org/records/12736824/files/UKBB_470k_deeprvat_results.csv?download=1"
deeprvat_file <- paste0(data_dir, "ukbb/UKBB_470k_deeprvat_results.csv")

if (!file.exists(deeprvat_file)) {
  message("⬇️ Downloading DeepRVAT results (this may take a few minutes)...")
  tryCatch({
    curl::curl_download(url = deeprvat_url, destfile = deeprvat_file, mode = "wb")
    message("✅ Download completed: ", basename(deeprvat_file))
  }, error = function(e) {
    message("❌ Download failed: ", e$message)
    if (file.exists(deeprvat_file)) file.remove(deeprvat_file)
    stop("Download of DeepRVAT results failed. Try again with a stable connection.")
  })
} else {
  message("✅ DeepRVAT results already exist: ", basename(deeprvat_file))
}

deeprvat_file <- paste0(data_dir, "ukbb/UKBB_470k_deeprvat_results.csv")

df_rvat <- read_csv(deeprvat_file)
dim(df_rvat)

# Examples - "Albumin", "Heart failure", "Mean reticulocyte volume".
df_rvat <- df_rvat |> filter(Trait == "Platelet distribution width")

dim(df_rvat)

# Prepare the data
df_rvat <- df_rvat |>
  filter(!is.na(pval) & pval > 0) |>
  mutate(logp = -log10(pval),
         gene = factor(`Gene symbol`, levels = unique(`Gene symbol`)),
         index = row_number())

# # Pan UKBB ----
print("The manifest of the Pan UKBB phenotypes manifest contaning updated study download links are stated here (<https://pan.ukbb.broadinstitute.org/downloads>): 
> `All heritability estimates (see here for more information on our approach) are available for download in the following formats:
Flat files: the manifest flat file is available on AWS (tarball here) or on Google Sheets. Our topline results can be found as part of the main phenotype manifest (Amazon AWS or Google Sheets)`
(<https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz>)
      ")

# Define output path
dest_gz <- paste0(data_dir, "ukbb/Pan_UK_Biobank_phenotype_manifest_manifest.tsv.bgz")
dest_tsv <- paste0(data_dir, "ukbb/Pan_UK_Biobank_phenotype_manifest_manifest.tsv")

# Download from AWS ----
if (!file.exists(dest_tsv)) {
  if (!file.exists(dest_gz)) {
    download.file(
      url = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz",
      destfile = dest_gz,
      mode = "wb"
    )
  }
  R.utils::gunzip(dest_gz, destname = dest_tsv, overwrite = TRUE)
} else {
  message("✅ Skipping download and decompression: ", basename(dest_tsv))
}

# Decompress .bgz to .tsv (R.utils::gunzip handles .bgz)
# R.utils::gunzip(dest_gz, destname = dest_tsv, overwrite = TRUE)

# Import using data.table::fread (faster) or readr::read_tsv
manifest <- fread(dest_tsv)
dim(manifest)
# Filter for "Albumin"
 
# Examples - phenocode == "411.8", phenocode == "30260" # Desc:  Mean reticulocyte volume.
filtered <- manifest |> filter(phenocode == "30110") # description == "Platelet distribution width"
t(filtered)

url <- filtered$aws_path

url <- gsub(
  "^s3://pan-ukb-us-east-1/", 
  "https://pan-ukb-us-east-1.s3.amazonaws.com/", 
  filtered$aws_path
)

outdir <- path.expand(paste0(data_dir, "ukbb"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

outfile <- file.path(outdir, basename(url))  # keeps full filename

if (file.exists(outfile)) {
  message("✅ Skipping, already downloaded: ", basename(outfile))
} else {
  message("⬇️ Downloading: ", basename(outfile))
  cmd <- sprintf(
    "curl -L -o %s --retry 3 --speed-limit 1000000 --speed-time 10 %s",
    shQuote(outfile), shQuote(url)
  )
  system(cmd)
}

# Read UKBB GWAS ----
library(readr)
library(dplyr)

# Raw SNPs: 28,987,535
filename <- basename(outfile)  # "continuous-30110-both_sexes-irnt.tsv.bgz"
subset_file <- file.path(dirname(outfile), sub("\\.tsv\\.bgz$", "_subset.tsv", filename))

# Downsample every 200th row + header.
if (!file.exists(subset_file)) {
  message("📂 Creating subset: ", basename(subset_file))
  cmd_subset <- sprintf(
    "gzcat %s | awk 'NR == 1 || NR %% 200 == 0' > %s",
    shQuote(outfile), shQuote(subset_file)
  )
  system(cmd_subset)
} else {
  message("✅ Subset already exists: ", basename(subset_file))
}

# Import subset using fread
library(data.table)

# entire GWAS ----
# Run this for a very slow complete test 
df_gwas <- fread(outfile,
                 select = c("chr", "pos", "ref", "alt", "beta_meta_hq", "neglog10_pval_meta"))

# downsampled GWAS ----
# Run this for an much faster test
# df_gwas <- fread(subset_file,
#                  select = c("chr", "pos", "ref", "alt", "beta_meta_hq", "neglog10_pval_meta"))

dim(df_gwas)

rm(filtered)
rm(manifest)

# Map SNPs to genes (variant to set_ID) ----

outfile <- paste0(data_dir, "ukbb/gene_coords_grch38.tsv")

if (file.exists(outfile)) {
  message("✓ gene_coords already downloaded. Reading from file.")
  gene_coords <- read_tsv(outfile, show_col_types = FALSE)
} else {
  message("Downloading gene coordinates from Ensembl...")
  
  # GRCh38 is the default; no need to specify GRCh
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  gene_coords <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position"),
    filters = "biotype",
    values = "protein_coding",
    mart = ensembl
  )
  
  # Clean and filter chromosomes before renaming
  gene_coords <- gene_coords |>
    filter(chromosome_name %in% c(as.character(1:22), "X", "Y")) |>
    mutate(chr = chromosome_name) |>
    dplyr::select(
      ensembl_gene_id,
      external_gene_name,
      chr,
      start = start_position,
      end = end_position
    )
  
  write_tsv(gene_coords, outfile)
  message("✓ gene coordinates saved to: ", outfile)
}

head(gene_coords)



# Prepare gene coordinates
dt_genes <- as.data.table(gene_coords)
dt_genes[, chr := paste0("chr", chr)]
setkey(dt_genes, chr, start, end)

# Prepare SNP data
dt_gwas <- as.data.table(df_gwas)
dt_gwas[, chr := paste0("chr", chr)]  # match 'chr' format with gene table
dt_gwas[, start := pos]
dt_gwas[, end := pos]
setkey(dt_gwas, chr, start, end)

# Perform overlap join: returns SNPs falling within gene ranges
df_gwas_mapped <- foverlaps(dt_gwas, dt_genes, nomatch = 0L)

df_gwas_mapped <- df_gwas_mapped |>
  dplyr::select(chr, pos, ref, alt, ensembl_gene_id, external_gene_name,
                beta_meta_hq, neglog10_pval_meta)
# View result
head(df_gwas_mapped)
head(df_rvat)

# Clean names ----
# Rename GWAS-mapped SNP-level data
df_gwas_mapped <- df_gwas_mapped |>
  dplyr::rename(
    CHR = chr,
    BP = pos,
    LOGP = neglog10_pval_meta,
    ENSEMBL = ensembl_gene_id,
    GENE = external_gene_name
  )

# Rename RVAT gene-level data
df_rvat <- df_rvat |>
  dplyr::rename(
    TRAIT = Trait,
    GENE = `Gene symbol`,
    ENSEMBL = `Ensembl ID`,
    P = pval
  )

df_gwas_mapped <- df_gwas_mapped |>
  dplyr::mutate(
    CHR = gsub("^chr", "", CHR),
    P = 10^(-LOGP)
  ) |>
  dplyr::select(-LOGP)

df_gwas_mapped$set_ID <- df_gwas_mapped$ENSEMBL
df_rvat$set_ID <- df_rvat$ENSEMBL

df_gwas_mapped <- df_gwas_mapped |> 
  dplyr::select(CHR, BP, P, set_ID, GENE) 

df_rvat <- df_rvat |> 
  dplyr::select(TRAIT, GENE, P, set_ID) 

# View result
head(df_gwas_mapped)
head(df_rvat)

# save archipelago-ready ----

write_tsv(df_rvat, paste0(data_dir, "ukbb/df1_UKBB_TRAIT1.tsv"))
write_tsv(df_gwas_mapped, paste0(data_dir, "ukbb/df2_UKBB_TRAIT1.tsv"))

rm(df_rvat, df_gwas, df_gwas_mapped)
rm(dt_genes, dt_gwas)

# Archipelago -----
## Check data function ----
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

# Real data ----
df1 <- fread(file.path(data_dir, "ukbb/df1_UKBB_TRAIT1.tsv"))
df2 <- fread(file.path(data_dir, "ukbb/df2_UKBB_TRAIT1.tsv"))

df1 <- df1 |> dplyr::select(set_ID, P)
df2 <- df2 |> dplyr::select(set_ID, BP, P, CHR)

# Down weight benign P-vaules to better viz ----

# Since the Pan UKBB gene-level study has so many VSAT tests
# Step 1: weight definition using your existing `P` column
df1_weighted <- df1 |>
  filter(!is.na(P)) |>
  mutate(
    logP = -log10(P + 1e-300),
    weight = 1 / (1 + exp(-(logP - 2)))
  )

# Step 2: split into top 10% and the rest
n_top <- round(nrow(df1_weighted) * 0.001)
df1_sorted <- df1_weighted |> arrange(P)
df1_top10 <- df1_sorted |> slice(1:n_top)
df1_rest <- df1_sorted |> slice((n_top + 1):n())

# Step 3: sample additional rows from the rest using weights
set.seed(42)
n_extra <- round(nrow(df1_weighted) * 0.01)  # or whatever % you'd like additionally

df1_rest_sampled <- df1_rest |>
  slice_sample(n = n_extra, weight_by = weight)

# Step 4: combine top and sampled rest, and sort
df1_sampled <- bind_rows(df1_top10, df1_rest_sampled) |>
  arrange(P)

df1 <- df1_sampled

# Special case -----
# Some RVAT do not have a SNP overlapping in GWAS since they are on rare and common variants, respectively.
# df1 (gene-level) contains all genes tested via RVAT
# df2 (SNP-level) only contains SNPs with GWAS p-values, and not all genes are tagged by nearby SNPs
# So: we can't always match on set_ID, but we can use genomic coordinates to find the nearest gene that has a SNP in df2, and link the gene-level p-value (df1) to that proxy set_ID.

# Some RVAT genes (df1) do not have overlapping SNPs in GWAS (df2)
# Use matched_snps to find a nearby gene ID (set_ID) from df2 to use instead

# matched_snps must contain: proxy_for = original gene ID from df1
# and set_ID = real gene ID from df2

# Prepare input data tables
dt_rvat <- as.data.table(df1)
dt_gwas <- as.data.table(df2)
dt_coords <- as.data.table(gene_coords)

# Standardise chromosome format
dt_coords[, chr := gsub("^chr", "", chr)]
dt_coords[, mid := (start + end) %/% 2]
dt_coords <- dt_coords[, .(set_ID = ensembl_gene_id, mid, chr)]

# Merge coordinates into RVAT
dt_rvat <- merge(dt_rvat, dt_coords, by = "set_ID", all.x = TRUE)
dt_rvat <- dt_rvat[!is.na(mid)]
dt_rvat <- dt_rvat[, .(set_ID, P, mid, chr)]

# Standardise chr format in GWAS
setnames(dt_gwas, "CHR", "chr")
dt_gwas[, chr := as.character(chr)]
dt_rvat[, chr := as.character(chr)]

# Set keys
setkey(dt_gwas, chr, BP)
setkey(dt_rvat, chr, mid)

# Identify missing genes and keep mid
missing_genes <- setdiff(dt_rvat$set_ID, dt_gwas$set_ID)
dt_missing <- dt_rvat[set_ID %in% missing_genes]
dim(dt_missing)

# Add a row index to preserve original identity
dt_missing[, row_id := .I]
dim(dt_missing)

# Nearest join
nearest_matches <- dt_gwas[dt_missing, on = .(chr, BP = mid), roll = "nearest", nomatch = 0L]
dim(nearest_matches)

# Rebuild mapping with row_id
df_proxymap <- nearest_matches[, .(
  row_id = row_id,
  original_set_ID = i.set_ID,
  proxy_set_ID = set_ID,
  rvat_P = i.P,
  proxy_chr = chr
)]

# Keep only closest match per row
df_proxymap <- df_proxymap[!duplicated(row_id)]

# Keep only closest match per gene
df_proxymap <- df_proxymap[!duplicated(original_set_ID)]

head(df_proxymap)

# Now you can join this back to the original `df_rvat` or `df2` to enable downstream plotting
# Filter out genes in df1 not already in df2
df1_missing <- df1 |> 
  filter(!set_ID %in% df2$set_ID)

# Join proxy mappings
df1_proxy_snps <- df1_missing |> 
  left_join(df_proxymap, by = c("set_ID" = "original_set_ID")) |> 
  filter(!is.na(proxy_set_ID))

# Pull in corresponding SNP rows from df2 using proxy_set_ID
df2_proxies <- df2 |> 
  filter(set_ID %in% df1_proxy_snps$proxy_set_ID) |> 
  left_join(df1_proxy_snps |> dplyr::select(set_ID, proxy_set_ID), by = c("set_ID" = "proxy_set_ID")) |> 
  mutate(set_ID = set_ID.y) |> 
  dplyr::select(names(df2))  # make sure column order matches

# Combine original df2 with proxy-augmented SNPs
df2_extended <- bind_rows(df2, df2_proxies) |> distinct()

# Replace original df2
df2 <- df2_extended

# Final filtering to keep only set_IDs present in df2
df1 <- dplyr::filter(df1, set_ID %in% df2$set_ID)

# Special case -----
df1 <- as.data.frame(df1)
df2 <- as.data.frame(df2)

# check_archipelago_input(df1, df2)
df2$CHR <- as.numeric(df2$CHR)
df2$BP <- as.numeric(df2$BP)
df2 <- na.omit(df2)
df1 <- df1 |> dplyr::filter(P > 0) # any Inf
df2 <- df2 |> dplyr::filter(P > 0) # any Inf

# cant have duplicate VSAT ?
df1 <- df1 |>
  dplyr::group_by(set_ID) |>
  dplyr::slice_min(P, n = 1, with_ties = FALSE) |>
  dplyr::ungroup()

df1 <- dplyr::filter(df1, set_ID %in% df2$set_ID)

check_archipelago_input(df1, df2)

print("Note that some top RVAT hit genes are not present as SNPs in GWAS directly")

# Install and run ----
library(archipelago)
# archipelago_plot(df1, df2) 

# Set the true count threshold for Bonferroni since we down-sample
count_snp <- 28987535
count_vsat <- 35968

# plot_subtitle <- paste0("Dataset: ", nrow(df1), " VSAT and ", format(nrow(df2), big.mark = ","), " SNPs")
plot_subtitle <- paste0(
  "Dataset: ", count_vsat, " VSAT and ",
  format(count_snp, big.mark = ","), " SNPs"
  )

plot_subtitle

plot_title <- "UKBB Trait 1"
plot_subtitle <- plot_subtitle

output_path = paste0(data_dir, "/output/archipelago_plot_ukbb")
output_raw = paste0(data_dir, "/output/vsat_raw_ukbb")

# output_path = paste0(data_dir, "/output/archipelago_plot_ukbb_full")
# output_raw = paste0(data_dir, "/output/vsat_raw_ukbb_full")

# remove the sporadit outliers that are distractingly high
df1_filt <- df1 |> filter(P > 1e-75)
df2_filt <- df2 |> filter(P > 1e-75)

archipelago_plot(df1 = df1_filt,
                 d = df2_filt,
                 add_title = TRUE,
                 plot_title = plot_title,
                 add_subtitle = TRUE,
                 plot_subtitle = plot_subtitle,
                 color_theme = 'retro',
                 output_path = output_path,
                 crit_val_VSAT = (0.05/count_vsat),
                 crit_val_single_variant = (0.05/count_snp),
                 # fig_width = 12,
                 # fig_height = 4,
                 annotate_thresholds = FALSE
)

df1_sig_named <- df1 |> 
  filter(P < (0.05 / count_vsat)) |>
  left_join(
    gene_coords |> 
      dplyr::select(
        ensembl_gene_id,
        external_gene_name,
        chr,
        start,
        end
      ),
    by = c("set_ID" = "ensembl_gene_id")
  ) |>
  arrange(P)

df2_sig_named <- df2 |> 
  # filter(P < (0.05 / count_snp)) |> # test
  left_join(
    gene_coords |> 
      dplyr::select(
        ensembl_gene_id,
        external_gene_name,
        chr,
        start,
        end
      ),
    by = c("set_ID" = "ensembl_gene_id")
  )

df2_sig_named <- df2_sig_named |> 
  group_by(set_ID) |>
  slice_min(P, with_ties = FALSE) |>
  ungroup() |>
  arrange(P)

write_tsv(df1_sig_named, paste0(data_dir, "ukbb/df1_UKBB_TRAIT1_hit_genes.tsv"))
write_tsv(df2_sig_named,  paste0(data_dir, "ukbb/df2_UKBB_TRAIT1_hit_genes.tsv"))

# Print overlaps for review ----
# Identify overlapping gene IDs
# overlap_ids <- intersect(df1_sig_named$external_gene_name, df2_sig_named$external_gene_name)
# 
# # Add source label for plotting
# df1_sig_named <- df1_sig_named |> mutate(source = "df1")
# df2_sig_named <- df2_sig_named |> mutate(source = "df2") |>
#   group_by(set_ID) |>
#   slice_min(P, with_ties = FALSE) |>
#   ungroup() |>
#   arrange(P)
# 
# # Combine and order by genomic position
# merged_df <- bind_rows(df1_sig_named, df2_sig_named) |>
#   mutate(chr = as.character(chr)) |>
#   arrange(as.numeric(gsub("chr", "", chr)), start) |>
#   mutate(index = row_number())
# 
# # Identify overlaps
# overlap_ids <- intersect(df1_sig_named$set_ID, df2_sig_named$set_ID)
# 
# # Plot
# ggplot(merged_df, aes(x = index, y = -log10(P), colour = source)) +
#   geom_point(alpha = 0.7) +
#   ggrepel::geom_text_repel(
#     data = merged_df |> filter(set_ID %in% overlap_ids),
#     aes(label = external_gene_name),
#     size = 3,
#     max.overlaps = Inf
#   ) +
#   theme_bw() +
#   labs(x = "Gene (ordered by genome)", y = expression(-log[10](P)), colour = "Source")

# recheck: ----
# Significance thresholds
threshold_df1 <- 0.05 / 35968
threshold_df2 <- 0.05 / 28987535

# Reattach gene names for df1
df1_sig_named <- df1 |>
  filter(P < threshold_df1) |>
  left_join(
    gene_coords |> dplyr::select(ensembl_gene_id, external_gene_name),
    by = c("set_ID" = "ensembl_gene_id")
  ) |>
  mutate(source = "df1")

# Reattach gene names for df2
df2_sig_named <- df2 |>
  filter(P < threshold_df2) |>
  group_by(set_ID) |>
  slice_min(P, with_ties = FALSE) |>
  ungroup() |>
  left_join(
    gene_coords |> dplyr::select(ensembl_gene_id, external_gene_name),
    by = c("set_ID" = "ensembl_gene_id")
  ) |>
  mutate(source = "df2")

# Identify overlaps using set_ID
overlap_genes <- inner_join(
  df1_sig_named |> dplyr::select(set_ID, external_gene_name),
  df2_sig_named |> dplyr::select(set_ID, external_gene_name),
  by = "set_ID",
  suffix = c("_rvat", "_gwas")
)

# Output unique overlapping gene names
unique_overlap_gene_names <- unique(overlap_genes$external_gene_name_rvat)
print(unique_overlap_gene_names)
