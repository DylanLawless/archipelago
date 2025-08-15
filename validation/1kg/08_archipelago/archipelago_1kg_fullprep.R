library(dplyr)
data_dir <- "~/mnt/atlas_data_big/data/archipelago_usage/"

# In this script we first import the results of GWAS and SKAT-O pathway.
# We merge GWAS positions on the file which matches SNP id to pathway ID.
# Then we can merge with SKAT-O pathway results using their common pathway ID. 

# Run test build
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::load_all("~/web/archipelago")

# get setIDs for SNPs ----
bim_set_pathway_path <- paste0(data_dir, "1kg/04_Data_QC/mcl_clusters_recode_df_whole_genome_v1_c7_700_1KG.EAS.tsv")

bim_set_pathway <- read.table(bim_set_pathway_path)
head(bim_set_pathway)
names(bim_set_pathway) <- c("set_ID", "ID")

# Import data from plink2. We used linear and logistic/Firth regression with covariates: 1kgeas.B1.glm.firth
df_snp <- read.table("~/web/archipelago_validation/06_Association_tests/1kgeas.B1.glm.firth", header=TRUE, comment.char="", check.names=FALSE)
head(df_snp)

df_snp_id <- merge(df_snp, bim_set_pathway, by = "ID")
head(df_snp_id)
# df_snp_id$set_ID <- as.numeric(df_snp_id$set_ID)

colnames(df_snp_id)[colnames(df_snp_id) == '#CHROM'] <- 'CHR'
colnames(df_snp_id)[colnames(df_snp_id) == 'POS'] <- 'BP'
df_snp_id <- df_snp_id |> dplyr::select(set_ID, BP, P, CHR)
df_snp_id$SNP <- df_snp_id$BP

# Import SKAT-O pathway analysis ----
df_pathway <- read.table("~/web/archipelago_validation/07_SKATO_pathway/SKAT0_results_sample_data.clean1kgeas_binary.txt.tsv", header=TRUE, comment.char="", check.names=FALSE)
df_pathway <- df_pathway[1:2]
head(df_pathway)
names(df_pathway) <- c("set_ID", "P")
df_pathway$set_ID <- as.numeric(df_pathway$set_ID)

df_pathway <- na.omit(df_pathway)
df_snp_id <- na.omit(df_snp_id)
df_pathway <- df_pathway |> filter(set_ID > 0)
df_snp_id <- df_snp_id |> filter(set_ID > 0)

print("Simulate a strong associated pathway for testing")
df_pathway_sim <- df_pathway
idx <- which.min(df_pathway_sim$P)
df_pathway_sim$P[idx] <- df_pathway_sim$P[idx] / 100000

# run achipelago on df1, df2 ----
output_path = "./1kgeas_archipelago_plot"
output_raw = "./1kgeas_vsat_raw_plot"

write.csv(df_snp_id, file = "./df_snp_id.csv", row.names = FALSE, quote = FALSE)
write.csv(df_pathway_sim, file = "./df_pathway_sim.csv", row.names = FALSE, quote = FALSE)

head(df_pathway_sim)
head(df_snp_id)

class(df_pathway_sim$set_ID)
class(df_snp_id$set_ID)

plot_subtitle <- paste0("Dataset: ", nrow(df_pathway_sim), " VSAT and ", format(nrow(df_snp_id), big.mark = ","), " SNPs")
plot_subtitle

plot_title <- "1KG Trait 1"
plot_subtitle <- plot_subtitle

archipelago_plot(df_pathway_sim,
                 df_snp_id,
                 add_title = TRUE,
                 plot_title = plot_title,
                 add_subtitle = TRUE,
                 plot_subtitle = plot_subtitle,
                 # color_theme = 'retro',
                 output_path = output_path,
                 output_raw = output_raw,
                 better_space = TRUE,
                 alpha_point = 0.8
)

# devtools::load_all("~/web/archipelago")
# archipelago_plot(df1, df2, better_space = TRUE)

# Get summary stats ----
head(bim_set_pathway)
head(df_snp_id)
head(df_pathway_sim)


# print(paste("For the GWAS, the allele frequency of enriched variants was", ))
library(dplyr)
library(ggplot2)

enriched <- df_snp %>% filter(P < 0.05/nrow(df_snp))
mean_freq <- mean(enriched$A1_FREQ, na.rm = TRUE)
median_freq <- median(enriched$A1_FREQ, na.rm = TRUE)
min_freq <- min(enriched$A1_FREQ, na.rm = TRUE)
max_freq <- max(enriched$A1_FREQ, na.rm = TRUE)

print(paste("For the GWAS, the allele frequency of enriched variants (n =", nrow(enriched), 
            ") was: mean =", round(mean_freq, 3),
            ", median =", round(median_freq, 3),
            ", min =", round(min_freq, 3),
            ", max =", round(max_freq, 3),
            "which indicates that these individual SNP associations were due to highly common variants."
            ))

# Get the set_ID corresponding to the minimum P in df_pathway_sim
min_set_id <- df_pathway_sim$set_ID[which.min(df_pathway_sim$P)]
subset_snp <- subset(df_snp_id, set_ID == min_set_id)
n_rows <- nrow(subset_snp)
print(paste("The top hit in VSAT was pathway ID:", min_set_id))
print(paste("The number of variants for the enriched VSAT was", n_rows))

# For subset_snp, match variants in df_snp by CHR and POS (df_snp_id$BP corresponds to df_snp$POS)
subset_enriched <- inner_join(df_snp, subset_snp, by = c("#CHROM" = "CHR", "POS" = "BP"))
n_enriched <- nrow(subset_enriched)
# print(paste("For the top VSAT hit (set_ID", min_set_id, "), there were", n_enriched, "variants with available allele frequency data."))

# Compute allele frequency metrics for the enriched VSAT variants
mean_freq_vsat <- mean(subset_enriched$A1_FREQ, na.rm = TRUE)
median_freq_vsat <- median(subset_enriched$A1_FREQ, na.rm = TRUE)
min_freq_vsat <- min(subset_enriched$A1_FREQ, na.rm = TRUE)
max_freq_vsat <- max(subset_enriched$A1_FREQ, na.rm = TRUE)

print(paste("For the top VSAT hit (set_ID", min_set_id, "), there were", n_enriched, "variants.",
"The allele frequency of variants was: mean =", round(mean_freq_vsat, 3),
            ", median =", round(median_freq_vsat, 3),
            ", min =", round(min_freq_vsat, 3),
            ", max =", round(max_freq_vsat, 3),
"This demonstrates the scenarios where rare variants within a collapsed set can have a different association than the GWAS region, although this is not necessary"))

p1 <- df_snp %>%
  filter(P < 0.05/nrow(df_snp)) %>%
  ggplot(aes(x = A1_FREQ)) +
  geom_histogram(binwidth = 0.05, fill = "#f6a192", color = "black") +
  labs(title = "Enriched variants in GWAS", x = "A1 Frequency", y = "Count") +
  theme_minimal()


p2 <- subset_enriched %>%
  # filter(P < 0.05/nrow(df_snp)) %>%
  ggplot(aes(x = A1_FREQ)) +
  geom_histogram(binwidth = 0.05, fill = "#27afea", color = "black") +
  labs(title = "Enriched variants in VSAT", x = "A1 Frequency", y = "Count")+
  theme_minimal()


library(patchwork)
patch1 <- p1 + p2

ggsave(filename = "plot_enriched_variant_frequencies.pdf", width = 8, height =  4)
