# Archipelago plot

<img src="man/figures/logo.png" width="200px" align="right"/>

## Summary and illustration of variant set association test statistics

Variant set association tests (VSAT), particularly those incorporating minor allele frequency variants, have become invaluable in genetic association studies by allowing robust statistical analysis with variant collapse. Unlike single variant tests, VSAT statistics cannot be assigned to a genomic coordinate for visual interpretation by default. To address these challenges, we introduce the Archipelago plot, a graphical method for interpreting both VSAT p-values and individual variant contributions. The Archipelago method assigns a meaningful genomic coordinate to the VSAT p-value, enabling its simultaneous visualization alongside individual variant p-values. This results in an intuitive and rich illustration akin to an archipelago of clustered islands, enhancing the understanding of both collective and individual impacts of variants. The Archipelago plot is applicable in any genetic association study that uses variant collapse to evaluate both individual variants and variant sets, and its customizability facilitates clear communication of complex genetic data. By integrating two dimensions of genetic data into a single visualization, VSAT results can be easily read and aid in identification of potential causal variants in variant sets such as protein pathways.

## How to

```r
# install the package
install.packages("~/archipelago_0.0.0.9000.tar.gz", repos = NULL, type = "source")

# load
library(archipelago)

# load test data
data("vsat_pval")
data("variant_pval")

# basic usage with defaults
archipelago_plot(
  df1 = vsat_pval,
  df2 = variant_pval
)

# specify a built-in colour theme
# available themes: 'retro', 'metro', 'summer', 'messenger', 'sunset', 'alice', 'buckley',
# 'romance', 'meme', 'saiko', 'pagliacci', 'ambush', 'sunra', 'caliber', 'yawn', 'lawless'
archipelago_plot(
  df1 = vsat_pval,
  df2 = variant_pval,
  color_theme = "alice"
)

# customised usage
output_path <- "./archipelago_plot_custom_color"  # file_type controls extension
output_raw  <- "./vsat_raw_plot"

color_labels  <- c("Label 1", "Label 2", "Label 3", "Label 4")
custom_colors <- c("#9abfd8", "#cac1f3", "#371c4b", "#2a5b7f")  # buckley palette

archipelago_plot(
  df1 = vsat_pval,                   # Input required VSAT level
  df2 = variant_pval,                # Input required variant level
  add_title = TRUE,
  plot_title = "Archipelago plot",
  add_subtitle = TRUE,
  plot_subtitle = "VSAT and single variant",
  show_legend = TRUE,                # Legend
  legend_position = "right",         # bottom top left right
  chr_ticks = FALSE,                 # chromosome ticks
  better_space = TRUE,               # prevent cramped chromosome labels
  color_theme = NULL,                # ignore built-in theme
  custom_colors = custom_colors,     # replace four default colors
  color_labels = color_labels,       # replace four default text in legend
  crit_val_VSAT = 0.05 / 300,        # VSAT threshold line
  crit_val_single_variant = 5e-8,    # typical single-variant threshold
  point_size = 0.5,
  point_size_large = 1.2,            # highlight larger points of interest
  alpha_point = 0.9,                 # point transparency
  alpha_seg = 0.5,                   # segment transparency
  fig_width = 8,                     # dimension of main figure
  fig_height = 4,                    # dimension of main figure
  raw_fig_width = 8,                 # dimension of raw VSAT figure
  raw_fig_height = 4,                # dimension of raw VSAT figure
  output_path = output_path,         # saved as output_path.file_type
  output_raw = output_raw,           # raw VSAT plot path
  file_type = "pdf"                  # png jpg pdf
)
```

## Significant threshold

The package expects `crit_val_VSAT` and `crit_val_single_variant` to be numeric P value cut-offs for drawing threshold lines. By default, if `crit_val_VSAT` is not provided, it is set to 0.05 / length(unique(df1$set_ID)), applying a Bonferroni correction based on the number of variant sets in df1. Likewise, if `crit_val_single_variant` is not provided, it can be set to 0.05 / length(unique(df2$set_ID)) to apply the same principle to single-variant results. This ensures that significant VSAT and single-variant results are highlighted consistently.

## Output images

The default output format is PNG. You can change `file_type` to `"pdf"`, `"jpg"`, or `"png"`. For large datasets, PDF files may be slow to render in a PDF viewer, whereas PNG or JPG formats usually display more smoothly.

## Validation datasets

Zenodo record: [https://doi.org/10.5281/zenodo.16880622](https://doi.org/10.5281/zenodo.16880622)

This repository contains validation studies for the Archipelago method, a framework for integrating single-variant GWAS and variant set association test (VSAT) statistics. Three datasets were used for validation:

* **1KG** — 1000 Genomes Project (East Asian samples), real GWAS trait and simulated pathway-level trait.
* **Pan-UKBB** — UK Biobank (platelet distribution width, DeepRVAT), real GWAS and gene-level trait.
* **UKBB WGS** — UK Biobank WGS (UTR collapsing PheWAS), real GWAS and gene-level trait.

Each dataset directory contains its own README, source data references, and output plots. The `src/` directory provides reproducible scripts for preparing inputs and generating Archipelago visualisations.

# Dev notes

## Set up
```
Rscript -e "usethis::create_package(getwd())"
```

## Docs with 
```
Rscript -e "devtools::document()"
```

## Build and install
```
Rscript -e "devtools::build()"
Rscript -e "devtools::install()"
```

## Install from zip in R
```
install.packages("/path/to/your/package/MyPackage.tar.gz", repos = NULL, type = "source")
```

## Example of set ID

An example of the protein pathway set IDs can be found at <https://github.com/DylanLawless/ProteoMCLustR/tree/main/data/ppi_examples> where we have a numeric pathway ID for each ensemble gene name, using MCL clustering on STRINGdb.



## Citation

Archipelago method for variant set association test statistics
Dylan Lawless, et al.
*medRxiv* 2025.03.17.25324111
[https://doi.org/10.1101/2025.03.17.25324111](https://doi.org/10.1101/2025.03.17.25324111)

