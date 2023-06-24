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
                             crit_val = NULL,
                             point_size = 1,
                             output_path = "archipelago_plot.pdf",
                             output_raw = "vsat_raw_plot.pdf"
                             ) {


    
    
library(ggplot2)
library(dplyr)

# Set data ----
# # Set seed for reproducibility
# set.seed(123)
# 
# mulitpier <- 5000
# first <- 1
# last <- mulitpier # Number of variants
# set_ID_max <- (mulitpier/20)
# crit_val <- .05/set_ID_max
# pval_max <- 1
# pval_min1 <- .05/10
# pval_min2 <- .05/100
# 
# # Sample df1: random numbers from a log-normal distribution
# df1 <- data.frame(
#   set_ID = seq(1, set_ID_max, by = 1),
#   P = rlnorm(set_ID_max, meanlog = log(pval_min1), sdlog = 1.6)
# )
# 
# df1$P <- -log10(df1$P)
# 
# # Sample df2: random from range
# df2 <- data.frame(
#   set_ID = rep(seq(1, set_ID_max, by = 1), each = 1),
#   BP = sample(1:1e6, last, replace = TRUE),
#   P = runif(last, pval_min2, pval_max),
#   CHR = sample(1:22, last, replace = TRUE)  # Chromosome numbers 1-22
# )
# 
# df2$SNP <- df2$BP # this would be required for multiallelic sites
# 
# write.csv(df1, file="../data/vsat_pval.txt", row.names=FALSE )
# write.csv(df2, file="../data/variant_pval.txt", row.names=FALSE)

# Bonferroni correction
# crit_val <- .05/set_ID_max
# Calculate crit_val from the dataset if it is NULL
if(is.null(crit_val)) {
  set_ID_max <- df1$set_ID %>% unique() %>% length()
  crit_val <- .05/set_ID_max
}

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

# Merging df1 and df2
merged_df <- bind_rows(df1, df2)

# Add a new column to specify variant_set points
merged_df$metric <- ifelse(is.na(merged_df$BP), "variant_set", "variant")

# Calculate the total of locations per set_ID
pos_sum <- aggregate(merged_df$pos, by=list(merged_df$set_ID), FUN=sum, na.rm=TRUE)

# Rank set_ID by the sum of locations
pos_sum$rank <- rank(pos_sum$x) 

# Calculate the average location within each set_ID group from the original dataframe
average_pos <- aggregate(merged_df$pos, by=list(merged_df$set_ID), FUN=mean, na.rm=TRUE)

# Replace NA locations with the average location of their set_ID group
merged_df$pos[is.na(merged_df$pos)] <- average_pos$x[match(merged_df$set_ID[is.na(merged_df$pos)], average_pos$Group.1)]

# Create a new grouping variable
merged_df$grouping <- with(merged_df, ave(pos, set_ID, FUN = function(x) cumsum(!is.na(x))))

# Get a simple ordered distribution of variant_set p value positions for the x-axis
# sort the "variant_set" group by "pos"
variant_set_sorted <- merged_df %>%
  filter(metric == "variant_set") %>%
  arrange(pos)

# create a sequence of evenly spaced numbers across the range of all "pos" values
even_pos <- seq(from = min(merged_df$pos),
                to = max(merged_df$pos),
                length.out = nrow(variant_set_sorted))

# assign these evenly spaced numbers to the "pos" column of the "variant_set" group
variant_set_sorted$pos <- even_pos

# replace the "pos" values in the original dataframe with these new evenly spaced numbers
merged_df <- merged_df %>%
  mutate(pos = ifelse(metric == "variant_set",
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
chromosome_ticks <- merged_df %>%
  filter(metric == "variant") %>%
  group_by(CHR) %>%
  summarise(mid_point = mean(pos))

# Split the data into variant_set and variant datasets
variant_set_data <- merged_df %>% filter(metric == "variant_set")
variant_data <- merged_df %>% filter(metric == "variant")

# Join the datasets to create a line dataset
line_data <- left_join(variant_data, variant_set_data, by = "set_ID", suffix = c("_variant", "_variant_set"))

# Raw plot ----
raw <- merged_df %>% 
  filter(metric == "variant_set")
raw$rank <- rank(raw$P)

p_raw <- 
  raw %>%
  ggplot() +
  geom_point(aes(x = rank, y = P), color='#27afea', size = point_size) +
  ylab("-log10 (p-value)") +
  geom_hline(linetype="dotted", 
             yintercept=-log10(crit_val)) +
  theme_bw() +
  ggtitle("Raw VSAT p-value") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.margin = margin(10, 10, 10, 10, "pt")) +
  guides(color = "none") 
p_raw
ggsave(p_raw, filename = output_raw, width = 6, height = 4)

# Plot 2 ----
# Highlight individual variant contributions 
# Create a new variable 'color_condition' that checks the condition
merged_df <- 
  merged_df %>%
  group_by(set_ID) %>%
  mutate(color_condition = 
           ifelse(any(P > -log10(crit_val)) & metric == "variant", 
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


p_arch_leg <- 
  ggplot() +
  geom_point(data = merged_df[merged_df$color_condition != "condition_met", ], 
             aes(x = pos, y = P, color = color_condition), size = point_size) +
  geom_segment(data = condition_met_lines, 
               aes(x = pos_variant_set, xend = pos_variant, y = P_variant_set, yend = P_variant), alpha = 0.3) +
  geom_point(data = condition_met_points, 
             aes(x = pos, y = P, color = color_condition), size = point_size, alpha = 1) +
  scale_color_manual(values = colors, 
                     labels= color_labels) +
  ylab("-log10 (p-value)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.margin = margin(10, 10, 10, 10, "pt")) 

# Add the title only if add_title is TRUE
if (add_title) {
  p_arch_leg <- p_arch_leg + ggtitle(plot_title)
}

# Add the subtitle only if add_subtitle is TRUE
if (add_subtitle) {
  p_arch_leg <- p_arch_leg + labs(subtitle = plot_subtitle)
}

# Only show chromosome ticks if chr_ticks is TRUE
if (chr_ticks) {
  p_arch_leg <- p_arch_leg + scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR) + xlab("Chromosome") 
}

# Only show the legend if show_legend is TRUE
if (!show_legend) {
  p_arch_leg <- p_arch_leg + guides(color = "none") 
}

# Set critical p-val line
if (!is.null(crit_val)) {
  p_arch_leg <- p_arch_leg + geom_hline(linetype="dotted", yintercept=-log10(crit_val))
}


p_arch_leg
ggsave(p_arch_leg, filename = output_path, width = 8, height = 4)

# return the final plot
return(p_arch_leg)
}

# # Import user data
# df1 <- read.csv(file="../data/vsat_pval.txt")
# df2 <- read.csv(file="../data/variant_pval.txt")
# save(df1, file = "../data/vsat_pval.RData")
# save(df2, file = "../data/variant_pval.RData")
# 
# # Use default settings
# archipelago_plot(df1, df2)
# 
# # 16 color themes:
# # 'retro', 'metro', 'summer', 'messenger', 'sunset', 'alice', 'buckley', 'romance', 'meme', 'saiko', 'pagliacci', 'ambush', 'sunra', 'caliber', 'yawn', 'lawless'
# output_path = "./archipelago_plot_custom_color.pdf"
# output_raw = "./vsat_raw_plot.pdf"
# archipelago_plot(df1,
#                  df2,
#                  color_theme = 'alice',
#                  output_path = output_path)
# 
# # Custom everything
# color_labels <- c("Label_1", "Label_2" ,"Label_3","Label_4")
# custom_colors = c("#9abfd8", "#cac1f3", "#371c4b", "#2a5b7f") # buckley theme colors
# plot_title <- "Title"
# plot_subtitle <- "Subtitle"
# crit_val = .05/300 # P-value threshold line
# point_size = .5 # geom_point size
# output_path = "../output/archipelago_plot_custom_everything.pdf"
# archipelago_plot(df1 = df1,
#                  d = df2,
#                  add_title = TRUE,
#                  plot_title = plot_title,
#                  add_subtitle = TRUE,
#                  plot_subtitle = plot_subtitle,
#                  show_legend = TRUE,
#                  crit_val = crit_val,
#                  point_size = point_size,
#                  chr_ticks <- FALSE, # TRUE / FALSE
#                  show_legend <- TRUE, # TRUE / FALSE
#                  custom_colors = custom_colors,
#                  color_labels = color_labels,
#                  output_path = output_path
#                  )
# 
# 
# # Print every color theme plot for manual
# # List of color themes
# color_themes <- c('retro', 'metro', 'summer', 'messenger', 'sunset', 'alice', 'buckley', 'romance', 'meme', 'saiko', 'pagliacci', 'ambush', 'sunra', 'caliber', 'yawn', 'lawless')
# 
# # Loop over each color theme
# for (color_theme in color_themes) {
#   # Define output path
#   output_path <- paste0("../output/archipelago_plot_", color_theme, ".pdf")
#   # Generate and save plot
#   archipelago_plot(df1,
#                    df2,
#                    color_theme = color_theme,
#                    output_path = output_path,
#                    show_legend = FALSE)
# }
