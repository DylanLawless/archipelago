# archipelago

[![CRAN status](https://www.r-pkg.org/badges/version/archipelago)](https://cran.r-project.org/package=archipelago)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/archipelago)](https://cran.r-project.org/package=archipelago)

<img src="man/figures/logo.png" width="200px" align="right"/>

Manhattan plots are for GWAS.  
Archipelago plots are for complex variant association studies.

`archipelago` provides an R implementation of the Archipelago plot, a visualisation method for integrating variant set association test statistics with single-variant association statistics in a shared genomic view.

The method is published in *Genetic Epidemiology*:

Lawless D, Saadat A, Ait Oumelloul M, Schlapbach LJ, Fellay J.  
**Archipelago Method for Variant Set Association Test Statistics.**  
*Genetic Epidemiology*. 2026;50(1):e70025.  
doi: [10.1002/gepi.70025](https://doi.org/10.1002/gepi.70025)

## Overview

Variant set association tests, including rare variant association tests and burden-based approaches, are widely used to test aggregated genetic effects across genes, pathways, regulatory regions, and other biologically defined sets.

Unlike single-variant GWAS results, variant set association test results do not naturally have a single genomic coordinate. This makes them difficult to visualise alongside single-variant association statistics.

The Archipelago method addresses this by assigning each variant set a representative genomic coordinate derived from its constituent variants. This enables set-level and variant-level association signals to be visualised together in a familiar genome-wide plot.

Archipelago does not perform association testing. It visualises results produced by existing GWAS, RVAT, VSAT, or related workflows.

## When to use Archipelago

Archipelago is useful when:

* You have variant set association results, such as gene-level, pathway-level, region-level, or burden test results.
* You also have single-variant association results.
* Both result tables share a common set identifier.
* You want to interpret aggregated association signals in genomic context.
* You want to identify whether a set-level signal is driven by one or more individual variants.

The method is applicable to focused disease cohorts, sequencing studies, population-scale biobanks, and any study design where variant collapse is used.

## Installation

Install the released version from CRAN:

```r
install.packages("archipelago")
````

Load the package:

```r
library(archipelago)
```

## Quick start

```r
library(archipelago)

data("vsat_pval", package = "archipelago")
data("variant_pval", package = "archipelago")

p_basic <- archipelago_plot(
  df1 = vsat_pval,
  df2 = variant_pval,
  output_path = tempfile(),
  output_raw  = tempfile()
)

p_basic
```

## Input data

`archipelago_plot()` requires two input data frames.

### Variant set association results

The first input, `df1`, contains set-level association results.

Required columns:

```r
set_ID
P
```

Example:

```r
head(vsat_pval)
```

```text
  set_ID           P
1      1 0.002039443
2      2 0.003459603
3      3 0.060544051
```

### Single-variant association results

The second input, `df2`, contains single-variant association results linked to the same set identifiers.

Required columns:

```r
set_ID
CHR
BP
P
SNP
```

Example:

```r
head(variant_pval)
```

```text
  set_ID     BP         P CHR    SNP
1      1 351696 0.9211610   6 351696
2      2 988282 0.8652950   9 988282
3      3 929171 0.6916336  12 929171
```

The `set_ID` column links the set-level and variant-level results.

## Basic usage

```r
p_basic <- archipelago_plot(
  df1 = vsat_pval,
  df2 = variant_pval,
  output_path = tempfile(),
  output_raw  = tempfile()
)

p_basic
```

## Colour themes

Built-in colour themes allow rapid visual changes.

```r
p_theme <- archipelago_plot(
  df1 = vsat_pval,
  df2 = variant_pval,
  color_theme = "alice",
  output_path = tempfile(),
  output_raw  = tempfile()
)

p_theme
```

Several predefined themes are available.

```r
"retro", "metro", "summer", "messenger", "sunset", "alice",
"buckley", "romance", "meme", "saiko", "pagliacci", "ambush",
"sunra", "caliber", "yawn", "lawless"
```

## Customised example

```r
custom_colors <- c("#9abfd8", "#cac1f3", "#371c4b", "#2a5b7f")

p_custom <- archipelago_plot(
  df1 = vsat_pval,
  df2 = variant_pval,
  add_title = TRUE,
  plot_title = "Custom Archipelago Plot",
  add_subtitle = TRUE,
  plot_subtitle = "Variant set and variant signals",
  show_legend = TRUE,
  legend_position = "bottom",
  chr_ticks = TRUE,
  point_size = 0.6,
  point_size_large = 1.2,
  custom_colors = custom_colors,
  color_labels = c(
    "Chromosome A",
    "Chromosome B",
    "Highlighted variants",
    "Variant set result"
  ),
  crit_val_VSAT = 0.05 / 300,
  crit_val_single_variant = 5e-8,
  annotate_thresholds = TRUE,
  fig_width = 10,
  fig_height = 5,
  output_path = tempfile(),
  output_raw  = tempfile(),
  file_type = "pdf"
)

p_custom
```

## Significance thresholds

`crit_val_VSAT` and `crit_val_single_variant` define the P value thresholds used to draw significance lines.

For set-level results, a common choice is a Bonferroni-style threshold based on the number of tested variant sets:

```r
0.05 / length(unique(df1$set_ID))
```

For single-variant results, a common choice is a genome-wide or variant-level threshold. For example:

```r
5e-8
```

or, for a Bonferroni-style correction over variants in the input table:

```r
0.05 / length(unique(df2$SNP))
```

Study-specific thresholds can be provided directly:

```r
crit_val_VSAT = 0.05 / 300
crit_val_single_variant = 5e-8
```

## Output files

Plots are saved automatically using the specified output paths and formats.

In examples and vignettes, temporary paths are used so files are not persisted:

```r
output_path = tempfile()
output_raw  = tempfile()
```

For saved output, provide explicit paths:

```r
output_path = "./archipelago_plot"
output_raw  = "./archipelago_raw_plot"
file_type   = "png"
```

Supported output formats are:

```r
"png"
"jpg"
"pdf"
```

For large genome-wide datasets, PDF files may be slow to render in some viewers. PNG or JPG output is usually preferable for rapid inspection, while PDF output may be preferred for publication workflows.

## Validation datasets

Validation datasets and reproducible scripts are available through Zenodo:

[https://doi.org/10.5281/zenodo.16880622](https://doi.org/10.5281/zenodo.16880622)

The validation repository contains three complementary examples:

* **1000 Genomes Project**: East Asian samples with a real GWAS trait and simulated pathway-level trait.
* **Pan-UK Biobank with DeepRVAT**: platelet distribution width using WES GWAS and gene-level rare variant association results.
* **UK Biobank WGS UTR PheWAS**: whole-genome sequencing with rare non-coding burden tests.

Each dataset directory contains its own README, source data references, processing scripts, and output plots.

## Related resources

* CRAN package: [https://cran.r-project.org/package=archipelago](https://cran.r-project.org/package=archipelago)
* Source code: [https://github.com/DylanLawless/archipelago](https://github.com/DylanLawless/archipelago)
* Method page: [https://switzerlandomics.ch/technologies/archipelago/](https://switzerlandomics.ch/technologies/archipelago/)
* Validation datasets: [https://doi.org/10.5281/zenodo.16880622](https://doi.org/10.5281/zenodo.16880622)

## Citation

If you use `archipelago`, please cite:

```bibtex
@article{2025lawlessArchipelagoMethodVariant,
  author  = {Lawless, Dylan and Saadat, Ali and Oumelloul, Mariam Ait and Schlapbach, Luregn J. and Fellay, Jacques},
  title   = {Archipelago Method for Variant Set Association Test Statistics},
  journal = {Genetic Epidemiology},
  volume  = {50},
  number  = {1},
  pages   = {e70025},
  year    = {2026},
  doi     = {10.1002/gepi.70025},
  url     = {https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.70025}
}
```

## Licence

`archipelago` is released under the MIT Licence.

It may be used, modified, and embedded in academic, clinical research, and commercial analysis pipelines, subject to the terms of the licence.

## Development

Generate documentation:

```r
Rscript -e "devtools::document()"
```

Build the package:

```r
Rscript -e "devtools::build()"
```

Install locally:

```r
Rscript -e "devtools::install()"
```

Install from a local source archive:

```r
install.packages(
  "/path/to/archipelago_0.0.0.9000.tar.gz",
  repos = NULL,
  type = "source"
)
```

## Example set identifiers

`set_ID` can represent any shared grouping used to connect variants and variant set statistics.

Examples include:

* Gene identifiers.
* Pathway identifiers.
* Protein interaction clusters.
* Regulatory regions.
* Sliding windows.
* User-defined variant groups.

An example of pathway-style set identifiers is available here:

[https://github.com/DylanLawless/ProteoMCLustR/tree/main/data/ppi_examples](https://github.com/DylanLawless/ProteoMCLustR/tree/main/data/ppi_examples)


