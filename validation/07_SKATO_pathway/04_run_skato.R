library(SKAT)
library(data.table)
library(dplyr)

setwd("~/web/archipelago_validation/07_SKATO_pathway")

# This requires ~ 160 GB memory and 1 hour.
# You can also split the VSAT by set_ID first and run each separately in parallel
# We have done approach in the past with SkatRbrain
# I think that once SSD is generated the memory drops to ~10GB so is worth saving.
# There is a comment by author in the sparse.R code noting this I think: Get Sparse Matrix From SSD file.


# 1 genotypeFile="../04_Data_QC/sample_data.clean"¬
# 2 phenotypeFile="../01_Dataset/1kgeas_binary.txt"¬
# 3 covariateFile="../05_PCA/plink_results_projected.sscore"¬

# setID_prefix <- "1KG.EAS.auto.snp.norm.nodup.split.rare002.common015.missing.bim.tsv"
# file_prefix <- "1KG.EAS.auto.snp.norm.nodup.split.rare002.common015.missing"

setID_path <- "../04_Data_QC/"
setID_prefix <- "sample_data.clean.bim.tsv"
base_path <- "../04_Data_QC/"
file_prefix <- "sample_data.clean"
pheno_path <- "../01_Dataset/"
pheno_name <- "1kgeas_binary.txt"
result_file <- paste0("SKAT0_results_", file_prefix, pheno_name)
covar <- "../05_PCA/plink_results_projected.sscore"

# Construct file paths
path_setID <- paste0(base_path, setID_prefix)
path_fam   <- paste0(base_path, file_prefix, ".fam")
path_bim   <- paste0(base_path, file_prefix, ".bim")
path_bed   <- paste0(base_path, file_prefix, ".bed")
path_ssd   <- paste0(base_path, file_prefix, ".SSD")
path_info  <- paste0(base_path, file_prefix, "SSD.info")

print(paste("path:",  path_setID))
print(paste("path:",  path_fam))
print(paste("path:",  path_bim))
print(paste("path:",  path_bed))
print(paste("path:",  path_ssd))
print(paste("path:",  path_info)) 

# temp <- read.table(file = path_setID)
# print(head(temp))
# temp_bim <- read.table(file = path_bim)
# print(head(temp_bim))

# merege SNP, gene, pathway ----
# genes
setID <- read.table(file = path_setID, header = FALSE)
colnames(setID)[colnames(setID) == 'V1'] <- 'Gene'

# pathways
pathwayID_path <- "~/web/archipelago/validation/07_SKATO_pathway/mcl_clusters_recode_df_whole_genome_v1_c7_700.csv"
pathwayID <- read.csv(pathwayID_path)
colnames(pathwayID)[colnames(pathwayID) == 'Items'] <- 'Gene'
pathwayID <- pathwayID |> select(ID, Gene)

bim <- read.table(file = path_bim)
bim <- bim |> select(V2)
bim_setID <- merge(bim, setID, all.x = TRUE)

# drop new duplicate variants (this is a problem with an example dataset)
print("this is slow but correct")
bim_setID <- bim_setID |>
  dplyr::group_by(V2) |>
  dplyr::slice(if(any(!is.na(Gene))) which(!is.na(Gene))[1] else 1) |>
  dplyr::ungroup()

head(bim_setID)
head(setID)
head(pathwayID)

bim_set_pathway <- merge(bim_setID, pathwayID, all.x = TRUE)
bim_set_pathway$ID[is.na(bim_set_pathway$ID)] <- "0"
bim_set_pathway <- bim_set_pathway |> select(ID, V2)
head(bim_set_pathway)
bim_set_pathway_path <- paste0(base_path, "mcl_clusters_recode_df_whole_genome_v1_c7_700_1KG.EAS.tsv")
write.table(bim_set_pathway, file = bim_set_pathway_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

read.table(path_setID) |> head()
read.table(bim_set_pathway_path) |> head()
read.table(path_setID) |> dim()
read.table(bim_set_pathway_path) |> dim()

# skat ----
# Load data and generate SSD and Info files if they don't exist
# Generate_SSD_SetID(path_bed, path_bim, path_fam, path_setID, path_ssd, path_info) # this would give single-variant tests (slowly)
Generate_SSD_SetID(path_bed, path_bim, path_fam, bim_set_pathway_path, path_ssd, path_info)

# Load external pheno - must be exact order and size of fam
fam_data <- read.table(path_fam)
names(fam_data) <- c("FID", "IID", "V3", "V4", "V5", "V6")
phenotype_data <- read.table(file = paste0(pheno_path, pheno_name), header = TRUE)
phenotype_data <- merge(fam_data, phenotype_data, by = c("FID", "IID"), all.x = TRUE) # ensure match and ordered
phenotype_data <- phenotype_data[[7]]
y.b <- ifelse(phenotype_data == 2, 1, 0)
print(y.b)

covar_data <- read.table(file = covar, header = FALSE) # not full match?
names(covar_data) <- c("FID", "IID",  "V3",  "V4",  "V5",  "V6",  "V7",  "V8",  "V9","V10", "V11", "V12", "V13", "V14")
fam_data_names <- fam_data |> select("FID", "IID")
covar_data <- merge(fam_data_names, covar_data, by = c("FID", "IID"), all.x = TRUE) # ensure match and ordered
print(head(covar_data))
X <- covar_data |> select(-"FID", -"IID") |> as.matrix()

# Open SSD
SSD.INFO <- Open_SSD(path_ssd, path_info)

# Define and run null model
obj_null <- SKAT_Null_Model(y.b ~ X, out_type = "D")
results <- SKAT.SSD.All(SSD.INFO, obj_null, method="SKATO")

print("Note for this dataset size we can get error : cannot allocate vector of size 67.7 Gb")

# Close SSD file
Close_SSD()

# Prepare output
output_df <- as.data.frame(results$results)

# Optionally, save the final results to a file
write.table(output_df, paste0(result_file, ".tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
output_df <- read.csv(file= paste0(result_file, ".tsv"), sep = "\t")

# Plot ----
library(ggplot2)
library(dplyr)
library(ggrepel)

# install.packages("ggrepel")

test_count <- nrow(output_df)
p_sig <- 0.05/test_count
p_sig <- round(p_sig, digits = 6)
p_sig_desc <- paste0("Psig threshold =\n", p_sig)
p_sig_desc
output_df$gene_count <- rownames(output_df) |> as.numeric()

names(output_df)
output_df$Chromosome <- output_df$SetID

# Assuming output_df is already loaded and contains the P.value and Gene_name
top_genes <- output_df %>%
  arrange(P.value) %>%
  top_n(-10, P.value)  # Selects the top 10 genes with the smallest P.value

p1 <- output_df %>%
  ggplot(aes(x = SetID, y = -log10(P.value), color = ifelse((Chromosome %% 2) == 0, "Even", "Odd"))) +
  geom_point() +
  geom_text_repel(data = top_genes, aes(label = SetID), box.padding = 1, point.padding = 0.3) +
  geom_text(label = p_sig_desc, x= (nrow(output_df)/22), y= -log10(p_sig), color = "black") +
  geom_hline(yintercept = -log10(p_sig), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("Even" = "black", "Odd" = "black")) +
  theme_minimal() +
  guides(color = FALSE) 
  # labs(title = "VSAT gene-level for rare variants in SwissPedHealth",
       # subtitle = paste0("WGS data, settings exp 2, variants with affect on coding genes, n=133.\n",
       # pheno_name))

p1
ggsave( paste0(result_file, ".pdf"), p1, width = 12, height = 8)

names(output_df)


print("simulate a strong associated pathway for testing")
top_genes_sim <- output_df %>%
  arrange(P.value) %>%
  top_n(-1, P.value)  # Selects the top 10 genes with the smallest P.value

top_genes_sim$P.value <- top_genes_sim$P.value/100000

p2 <- output_df %>%
  ggplot(aes(x = SetID, y = -log10(P.value), color = ifelse((Chromosome %% 2) == 0, "Even", "Odd"))) +
  geom_point() +
  geom_point(data = top_genes_sim) +
  geom_text_repel(data = top_genes_sim, aes(label = SetID),  box.padding = 1, point.padding = 0.3) +
  geom_text(label = p_sig_desc, x= (nrow(output_df)/22), y= -log10(p_sig), color = "black") +
  geom_hline(yintercept = -log10(p_sig), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("Even" = "black", "Odd" = "black")) +
  theme_minimal() +
  guides(color = FALSE) 
  # labs(title = "VSAT gene-level for rare variants in SwissPedHealth",
       # subtitle = paste0("WGS data, settings exp 2, variants with affect on coding genes, n=133.\n",
                         # pheno_name))

p2
ggsave( paste0(result_file, "_simulated_sig.pdf"), p1, width = 12, height = 8)

