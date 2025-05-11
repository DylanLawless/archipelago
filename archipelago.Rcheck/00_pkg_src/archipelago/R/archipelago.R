#' Archipelago Plot
#'
#' Plot p-values from variant and variant-set tests.
#'
#' @name archipelago_plot
#' @param df1 A dataframe of variant-set data.
#' @param df2 A dataframe of variant data.
#' @param plot_title Title for the plot.
#' @param add_title Logical; add title if TRUE.
#' @param plot_subtitle Subtitle for the plot.
#' @param add_subtitle Logical; add subtitle if TRUE.
#' @param chr_ticks Logical; show chromosome ticks if TRUE.
#' @param show_legend Logical; display the legend if TRUE.
#' @param color_theme Name of the colour theme.
#' @param custom_colors Vector of custom colours.
#' @param color_labels Labels for the colour groups.
#' @param crit_val Critical p-value threshold.
#' @param point_size Size of the points.
#' @param point_size_large Size of the large points of interest such as VSAT.
#' @param fig_width Width of the archipelago plot.
#' @param fig_height Height of the archipelago plot.
#' @param raw_fig_width Width of the raw plot.
#' @param raw_fig_height Height of the raw plot.
#' @param output_path File path for the plot.
#' @param output_raw File path for the raw plot.
#' @param file_type Use .png, .jpg, or .pdf. Defaults to png. pdf is slow and large for many SNPs.
#' @param alpha_point Use the alpha_point value to set point transparency. 
#' @param alpha_seg Use the alpha_segvalue to set line segment transparency. 
#' @param better_space Use better_space = TRUE to make sure that x-axis chr do not squash. 
#' @param legend_position Default right, allows bottom top left right.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'   # Load example data for df1 (VSAT) and df2 (individual variant)
#'   data("vsat_pval")         # same structure as df1
#'   data("variant_pval")      # same structure as df2
#'
#'   # Basic usage with defaults
#'   archipelago_plot(
#'     df1 = vsat_pval,
#'     df2 = variant_pval
#'   )
#'
#'   # Specify a built-in colour theme
#'   archipelago_plot(
#'     df1 = vsat_pval,
#'     df2 = variant_pval,
#'     color_theme = "alice"
#'   )
#'
#'   # More customised usage
#'   output_path <- "./archipelago_plot_custom_color.pdf"
#'   output_raw <- "./vsat_raw_plot.pdf"
#'   color_labels <- c("Label_1", "Label_2", "Label_3", "Label_4")
#'   custom_colors <- c("#9abfd8", "#cac1f3", "#371c4b", "#2a5b7f")
#'
#'   archipelago_plot(
#'     df1 = vsat_pval,
#'     df2 = variant_pval,
#'     add_title = TRUE,
#'     plot_title = "My Archipelago Plot",
#'     add_subtitle = TRUE,
#'     plot_subtitle = "VSAT vs Single Variant",
#'     show_legend = TRUE,
#'     chr_ticks = FALSE,
#'     point_size = 0.5,
#'     color_theme = NULL,         # ignore built-in theme
#'     custom_colors = custom_colors,
#'     color_labels = color_labels,
#'     crit_val_VSAT = 0.05 / 300, # for highlighting VSAT p-values
#'     crit_val_single_variant = 5e-8, # typical single-variant threshold
#'     output_path = output_path, 
#'     output_raw = output_raw,   
#'     file_type = "pdf"          # save as PDF instead of PNG
#'   )
#' }
utils::globalVariables(c(
  "aes", "color_group", "metric", "CHR", "BP", "width",
  "new_pos_in_chr", "CHR_order", "pos", "set_ID", "P",
  "color_condition", "pos_variant_set", "pos_variant",
  "P_variant_set", "P_variant", "quantile"
))
#'
#' @export

archipelago_plot <- function(df1, df2, 
                             plot_title = "Archipelago Plot", 
                             add_title = FALSE,
                             plot_subtitle = "Variant Set Association Test\nwith individual variant contributions\nand contributions for significant VSAT",
                             add_subtitle = FALSE,
                             chr_ticks = TRUE, 
                             show_legend = TRUE,
                             color_theme = NULL,
                             custom_colors = NULL,
                             color_labels = c("Chr: Individual p-val", "Chr: Individual p-val" ,"Individual variant\nfrom enriched VSAT","VSAT p-val"),
                             crit_val_VSAT = NULL,
                             crit_val_single_variant = NULL,
                             point_size = 1,
                             point_size_large = 1,
                             fig_width = 8,
                             fig_height = 4,
                             raw_fig_width = 8,
                             raw_fig_height = 4,
                             output_path = "archipelago_plot",
                             output_raw = "archipelago_vsat_raw_plot",
                             file_type = "png",
                             alpha_point = 1,
                             alpha_seg = 0.3,
                             better_space = FALSE,
                             legend_position = "right"
                             ) {

# library(ggplot2)
# library(dplyr)

  # Output file naming
  if (!grepl(paste0("\\.", file_type, "$"), output_raw)) {
    output_raw <- paste0(output_raw, ".", file_type)
  }
  
  if (!grepl(paste0("\\.", file_type, "$"), output_path)) {
    output_path <- paste0(output_path, ".", file_type)
  }
  
  
print("Input df1 is for VSAT: set_ID and P")
print("Input df2 is for SNP: set_ID, BP, P, CHR, SNP")


if (!"P" %in% colnames(df1) || nrow(df1) < 1) {
  stop("df1 does not have a 'P' column or is empty.")
}
if (!"P" %in% colnames(df2) || nrow(df2) < 1) {
  stop("df2 does not have a 'P' column or is empty.")
}

# Ensure P columns are numeric
df1$P <- as.numeric(df1$P)
df2$P <- as.numeric(df2$P)


# The example data was from plink and skat. It was then saved to Rdata with
# write.csv(df1, file = "vsat_pval.csv", row.names = FALSE)
# save(vsat_pval, file = "vsat_pval.RData")
# alternative: usethis::use_data(vsat_pval, overwrite = TRUE)

# if(is.null(crit_val)) {
#   set_ID_max <- df1$set_ID |> unique() |> length()
#   crit_val <- .05/set_ID_max
# }

if(is.null(crit_val_VSAT)) {
    set_ID_max <- df1$set_ID |> unique() |> length()
    crit_val_VSAT <- .05/set_ID_max
  }

if(is.null(crit_val_single_variant)) {
  set_ID_max <- df2$BP |> unique() |> length()
  crit_val_single_variant <- .05/set_ID_max
}

df1$P <- -log10(df1$P)
df2$P <- -log10(df2$P)
df2 <- df2[order(df2$CHR, df2$BP), ]
df2$index=NA
ind = 0
for (i in unique(df2$CHR)){
  ind = ind + 1
  df2[df2$CHR==i,]$index = ind
}
df2$index = rep.int(seq_along(unique(df2$CHR)), times = tapply(df2$BP,df2$CHR,length))  

nchr = length(unique(df2$CHR))
df2$pos=NA

# The following code block replaces my old method - to more closely match qqman - to make it easier to integrate in future.
lastbase=0
ticks=NULL
for (i in unique(df2$index)) {
  if (i==1) {
    df2[df2$index==i, ]$pos=df2[df2$index==i, ]$BP
  } else {
    lastbase = lastbase +max(df2[df2$index==(i-1),"BP"])  
    df2[df2$index == i,"BP"] = df2[df2$index == i,"BP"]-min(df2[df2$index==i,"BP"]) +1
    df2[df2$index == i, "pos"] = df2[df2$index == i,"BP"] + lastbase 
    
  }
}
ticks <-tapply(df2$pos,df2$index,quantile,probs=0.5) 
xlabel = 'Chromosome'
labs <- unique(df2$CHR)

# This function prevents the chromosome ticks from squashing and overlapping. It uses the largest chr to space all others while keeping the structure. 
if (better_space) {
  # Find the maximum BP range across all chromosomes
  max_width <- df2 |>
    dplyr::group_by(CHR) |>
    dplyr::summarise(width = max(BP) - min(BP)) |>
    dplyr::ungroup() |>
    dplyr::summarise(max_width = max(width)) |>
    dplyr::pull(max_width)
  
  # Scale each chromosome's positions so that they span max_width and add a cumulative offset
  df2 <- df2 |>
    dplyr::group_by(CHR) |>
    dplyr::arrange(BP) |>
    dplyr::mutate(new_pos_in_chr = (BP - min(BP)) / (max(BP) - min(BP)) * max_width) |>
    dplyr::ungroup() |>
    dplyr::mutate(CHR_order = as.numeric(factor(CHR, levels = sort(unique(CHR)))),
           new_pos = new_pos_in_chr + (CHR_order - 1) * max_width)
  
  df2$pos <- df2$new_pos
}

# Merging df1 and df2
merged_df <- dplyr::bind_rows(df1, df2)

# Add a new column to specify variant_set points
merged_df$metric <- ifelse(is.na(merged_df$BP), "variant_set", "variant")

# Calculate the total of locations per set_ID
pos_sum <- stats::aggregate(merged_df$pos, by=list(merged_df$set_ID), FUN=sum, na.rm=TRUE)

# Rank set_ID by the sum of locations
pos_sum$rank <- rank(pos_sum$x)

# Calculate the average location within each set_ID group from the original dataframe
average_pos <- stats::aggregate(merged_df$pos, by=list(merged_df$set_ID), FUN=mean, na.rm=TRUE)

# Replace NA locations with the average location of their set_ID group
merged_df$pos[is.na(merged_df$pos)] <- average_pos$x[match(merged_df$set_ID[is.na(merged_df$pos)], average_pos$Group.1)]

# Create a new grouping variable
merged_df$grouping <- with(merged_df, ave(pos, set_ID, FUN = function(x) cumsum(!is.na(x))))

# Get a simple ordered distribution of variant_set p value positions for the x-axis
# sort the "variant_set" group by "pos"
variant_set_sorted <- merged_df |>
  dplyr::filter(metric == "variant_set") |>
  dplyr::arrange(pos)

# create a sequence of evenly spaced numbers across the range of all "pos" values
even_pos <- seq(from = min(merged_df$pos),
                to = max(merged_df$pos),
                length.out = nrow(variant_set_sorted))

# assign these evenly spaced numbers to the "pos" column of the "variant_set" group
variant_set_sorted$pos <- even_pos

# replace the "pos" values in the original dataframe with these new evenly spaced numbers
merged_df <- merged_df |>
  dplyr::mutate(pos = ifelse(metric == "variant_set",
                      variant_set_sorted$pos[match(set_ID, variant_set_sorted$set_ID)],
                      pos))

# The previous 4 code chunks constructs a natural spread of positions for "variant_set". However, their positions are distributioned slightly atrificially to give a more even distribution in their true odeder. This is ideal for a dense network. If you have [1] a less dense network, [2] don't mind high-density clustiner at the genome center, [3] want to see the less dispersed distribution, just uncomment the next chuck to override the previous distribution version.

# merged_df$pos <- ifelse(merged_df$metric == "variant_set",
# 								(merged_df$pos - min(merged_df$pos[merged_df$metric == "variant_set"])) /
# 									(max(merged_df$pos[merged_df$metric == "variant_set"]) - min(merged_df$pos[merged_df$metric == "variant_set"])) * max(df2$pos),
# 								merged_df$pos)

# Add a new column for color groups
merged_df$color_group <- ifelse(merged_df$metric == "variant_set", "variant_set", 
                                ifelse(merged_df$CHR %% 2 == 0, "Chr_B", "Chr_A"))

# Create a dataframe with the mid-point of each chromosome
chromosome_ticks <- merged_df |>
  dplyr::filter(metric == "variant") |>
  dplyr::group_by(CHR) |>
  dplyr::summarise(mid_point = mean(pos))

# # Instead of using mean(pos) for chromosome_ticks, use the midpoint of the range of positions:
# chromosome_ticks <- merged_df |>
#   dplyr::filter(metric == "variant") |>
#   dplyr::group_by(CHR) |>
#   dplyr::summarise(mid_point = (min(pos) + max(pos)) / 2)

# Split the data into variant_set and variant datasets
variant_set_data <- merged_df |> dplyr::filter(metric == "variant_set")
variant_data <- merged_df |> dplyr::filter(metric == "variant")

# Join the datasets to create a line dataset
line_data <- dplyr::left_join(variant_data, variant_set_data, by = "set_ID", suffix = c("_variant", "_variant_set"))

# Raw plot ----
raw <- merged_df |> 
  dplyr::filter(metric == "variant_set")
raw$rank <- rank(raw$P)

p_raw <- 
  raw |>
  ggplot2::ggplot() +
  ggplot2::geom_point(
    ggplot2::aes(x = rank, y = P), color='#27afea', size = point_size) +
  ggplot2::ylab("-log10 (p-value)") +
  ggplot2::geom_hline(linetype="dotted", 
             yintercept=-log10(crit_val_VSAT)) +
  ggplot2::theme_bw() +
  ggplot2::ggtitle("Raw VSAT p-value") + 
  ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "black"),
        plot.margin = ggplot2::margin(10, 10, 10, 10, "pt")) +
  ggplot2::guides(color = "none") 
p_raw
ggplot2::ggsave(p_raw, filename = output_raw, width = raw_fig_width, height = raw_fig_height, device = file_type)

# Plot 2 ----
# Highlight individual variant contributions 
# Create a new variable 'color_condition' that checks the condition

# Step 1: Identify all set_ID groups where the metric is "variant_set" and P > -log10(crit_val_VSAT)
set_IDs_to_color <- merged_df |>
  dplyr::filter(metric == "variant_set", P > -log10(crit_val_VSAT)) |>
  dplyr::pull(set_ID)

# Step 2: Assign "condition_met" to "variant" points within those groups
merged_df <- 
  merged_df |>
  dplyr::mutate(color_condition = 
           ifelse(metric == "variant" & set_ID %in% set_IDs_to_color, 
                  "condition_met", 
                  color_group))

# Change the alpha variable accordingly
merged_df$alpha <- ifelse(merged_df$color_condition == "condition_met", 1, ifelse(merged_df$metric == "variant_set", 1, 0.3))


# Plot 4 ----
# Clearer condition_met layer
# Separate the points and lines that meet the condition
condition_met_points <- merged_df[merged_df$color_condition == "condition_met", ]
condition_met_lines <- line_data[line_data$set_ID %in% condition_met_points$set_ID, ]

# Figure  ----
# Define color themes
color_themes <- list(
  retro = c("#ffe28a", "#6fcb9f", "#666547", "#fb2e01"),
  metro = c("#f37735", "#ffc425", "#d11141", "#00aedb"),
  summer = c("#ddf098", "#f9d62e", "#ff4e50", "#fc913a"),
  messenger = c("#44bec7", "#ffc300", "#0084ff", "#fa3c4c"),
  sunset = c("#eeaf61", "#fb9062", "#6a0d83", "#ee5d6c"),
  alice = c("#f7c297", "#ffecb8", "#90d2d8", "#f6a6b2"),
  buckley = c("#9abfd8", "#cac1f3", "#371c4b", "#2a5b7f"), 
  romance = c("#ffcad4", "#c08497", "#3a4440", "#cc8562"),
  meme = c("#f4d4b8", "#a9dada", "#474747", "#f17255"),
  saiko = c("#a2b3d8", "#6289b3", "#d0444a", "#005e90"),
  pagliacci = c("#ffdd75", "#59a4d2", "#f7955d", "#49518a"),
  ambush = c("#f5e9be", "#ffe184", "#174c4f", "#207178"),
  sunra = c("#34b1ff", "#2d93d6", "#eaa221", "#1a406f"),
  caliber = c("#cccccc", "#bbbbbb", "#900303", "#000000"),
  yawn = c("#bcbcbc", "#999999", "#000000", "#5b5b5b"),
  lawless = c('#f6d992', '#f6a192', 'black', '#27afea')
)

# Check if custom_colors is provided, otherwise use color_theme or default
if (!missing(custom_colors)) {
  colors <- custom_colors
} else if (!missing(color_theme) && color_theme %in% names(color_themes)) {
  colors <- color_themes[[color_theme]]
} else {
  colors <- color_themes[['lawless']] # default
}


# Create a data frame with only the "variant_set" points
variant_set_points <- merged_df[merged_df$metric == "variant_set", ]

p_arch_leg <- 
  ggplot2::ggplot() +
  # Adding points where color_condition is not "condition_met"
  ggplot2::geom_point(data = merged_df[merged_df$color_condition != "condition_met", ], 
             ggplot2::aes(x = pos, y = P, color = color_condition), size = point_size, alpha = alpha_point) +
  # Adding lines
  ggplot2::geom_segment(data = condition_met_lines, 
               ggplot2::aes(x = pos_variant_set, xend = pos_variant, y = P_variant_set, yend = P_variant), alpha = alpha_seg) +
  # Adding points where metric is "variant_set", always on top of SNP
  ggplot2::geom_point(data = variant_set_points, 
             ggplot2::aes(x = pos, y = P, color = color_condition), size = point_size_large, alpha = alpha_point) +
  # Adding points where color_condition is "condition_met"
  ggplot2::geom_point(data = condition_met_points, 
             ggplot2::aes(x = pos, y = P, color = color_condition), size = point_size, alpha = alpha_point) +
   ggplot2::scale_color_manual(values = colors, 
                     labels= color_labels) +
  ggplot2::ylab("-log10 (p-value)") +
  ggplot2::theme_bw() +
  ggplot2::labs(colour = "P value") +
  ggplot2::theme(legend.position = legend_position,
    panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "black"),
        plot.margin = ggplot2::margin(10, 10, 10, 10, "pt")) 

# Add the title only if add_title is TRUE
if (add_title) {
  p_arch_leg <- p_arch_leg + ggplot2::ggtitle(plot_title)
}

# Add the subtitle only if add_subtitle is TRUE
if (add_subtitle) {
  p_arch_leg <- p_arch_leg + ggplot2::labs(subtitle = plot_subtitle)
}

# Only show chromosome ticks if chr_ticks is TRUE
if (chr_ticks) {
  p_arch_leg <- p_arch_leg + ggplot2::scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR) + ggplot2::xlab("Chromosome") 
}

# Only show the legend if show_legend is TRUE
if (!show_legend) {
  p_arch_leg <- p_arch_leg + ggplot2::guides(color = "none") 
}

# Set critical p-val line
if (!is.null(crit_val_VSAT)) {
  p_arch_leg <- p_arch_leg + 
    ggplot2::geom_hline(linetype="dotted", color = "#135775",
               yintercept=-log10(crit_val_VSAT)) +
    ggplot2::geom_hline(linetype="dotted", color = "#f6a192",
               yintercept=-log10(crit_val_single_variant))
  
}

# label each p-val line
if (!is.null(crit_val_VSAT)) {
  p_arch_leg <- p_arch_leg + 
    ggplot2::geom_text(ggplot2::aes(x = max(merged_df$pos), y = (-log10(crit_val_VSAT)+0.3), 
                  label = "VSAT\nthreshold"), hjust = 1.2) +
    ggplot2::geom_text(ggplot2::aes(x = max(merged_df$pos), y = (-log10(crit_val_single_variant)+0.3), 
                  label = "single-variant\nthreshold"), hjust = 1.2)
}

p_arch_leg
ggplot2::ggsave(p_arch_leg, filename = output_path, width = fig_width, height = fig_height, device = file_type)


# return the final plot
return(p_arch_leg)
}


# library(patchwork)
# # patch1 <- (p_arch_leg / p_raw)
# patch1 <- (p_arch_leg | p_raw) + plot_annotation(tag_levels = 'A')
# 
# ggplot2::ggsave(patch1, filename = output_patch_jpg, width = 14, height = 9)
# 
# # return the final plot
# return(p_arch_leg)
# }
# 
# # archipelago_plot(df1, df2)
# archipelago_plot(df1, df2_ready)
