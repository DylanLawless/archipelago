# Archipelago Plot ----
library(ggplot2)
# Make synthetic data. 
# Pathway P-value is likely lower than variant P-values. 
# Set seed for reproducibility
set.seed(123)

# To do - replace coordinates with genomic Chromosome positions. 
# Assign rank per variant for simple plotting x-axis.

mulitpier <- 1
first <- 1
last <- 5000*mulitpier
MCLID_max <- 20
crit_val <- 0.0002
pval_max <- 0.001

alternating_colors_max <- 500*mulitpier
alternating_colors_min <- 250*mulitpier

# Sample df1: random numbers from a log-normal distribution
df1 <- data.frame(
  MCLID = seq(1, MCLID_max, by = 1),
  p_value = rlnorm(MCLID_max, meanlog = log(0.001), sdlog = 1)
)

df1$location <- NA

# Sample df2: random from range
df2 <- data.frame(
  MCLID = rep(seq(1, MCLID_max, by = 1), each = 1),
  location = first:last,
  p_value = runif(last, 0.001, 0.01)
)

# Merging df1 and df2
merged_df <- merge(df1, df2, all = T)

# Add a new column to specify pathway points
merged_df$metric <- ifelse(is.na(merged_df$location), "pathway", "variant")

# Calculate the total of locations per MCLID
location_sum <- aggregate(merged_df$location, by=list(merged_df$MCLID), FUN=sum, na.rm=TRUE)

# Rank MCLID by the sum of locations
location_sum$rank <- rank(location_sum$x) 

# Calculate the average location within each MCLID group from the original dataframe
average_location <- aggregate(merged_df$location, by=list(merged_df$MCLID), FUN=mean, na.rm=TRUE)

# Replace NA locations with the average location of their MCLID group
merged_df$location[is.na(merged_df$location)] <- average_location$x[match(merged_df$MCLID[is.na(merged_df$location)], average_location$Group.1)]

# Create a new grouping variable
merged_df$grouping <- with(merged_df, ave(location, MCLID, FUN = function(x) cumsum(!is.na(x))))

# Transform the p_value column
merged_df$p_value <- -log10(merged_df$p_value)


# Narmalise the MCLID location genome wide

merged_df$location <- ifelse(merged_df$metric == "pathway", 
                             (merged_df$location - min(merged_df$location[merged_df$metric == "pathway"])) / 
                               (max(merged_df$location[merged_df$metric == "pathway"]) - min(merged_df$location[merged_df$metric == "pathway"])) * last,
                             merged_df$location)

library(dplyr)

# Separate the data
pathway_points <- merged_df[merged_df$grouping == 1, ]
variant_points <- merged_df[merged_df$grouping > 1, ]

# Join the data
line_data <- variant_points %>%
  left_join(pathway_points, by = "MCLID", suffix = c("", "_pathway"))



# Add a new column for color groups
merged_df$color_group <- ifelse(merged_df$metric == "pathway", "Pathway", 
                                ifelse(((merged_df$location-1) %% alternating_colors_max) < alternating_colors_min,"Chr_A", "Chr_B"))

# Plotting
ggplot() +
  geom_segment(data = line_data, 
               aes(x = location_pathway, xend = location, y = p_value_pathway, yend = p_value), alpha = 0.01) +
  geom_point(data = merged_df, aes(x = location, y = p_value, color = color_group), size = 3, alpha = ifelse(merged_df$metric == "pathway", 1, 0.5)) +
  scale_color_manual(values = c('#f6d992', '#f6a192', '#71c7ec')) +
  ylab("-log10 (p-value)") +
  xlab("Genome position") +
  geom_hline(linetype="dotted", 
             yintercept=-log10(crit_val)) +
  theme_bw() +
  ggtitle("Archipelago Plot") + 
  labs(subtitle = "Variant Set Association Test (blue)\nwith variant level contributions (alternating yellow/orange).") +
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      plot.margin = margin(10, 10, 10, 10, "pt")) +
  guides(color = "none")

# Archipelago_plot_1_20_5k.pdf 6x4


# Version 2 - highlight variant individual contributions ----

# Create a new variable 'color_condition' that checks the condition
merged_df <- 
  merged_df %>%
  group_by(MCLID) %>%
  mutate(color_condition = 
           ifelse(any(p_value > -log10(crit_val)) & metric == "variant", 
                  "condition_met", 
                  color_group))

# Change the alpha variable accordingly
merged_df$alpha <- ifelse(merged_df$color_condition == "condition_met", 1, ifelse(merged_df$metric == "pathway", 1, 0.3))

# Plotting
ggplot() +
  geom_segment(data = line_data, 
               aes(x = location_pathway, xend = location, y = p_value_pathway, yend = p_value), alpha = 0.01) +
  geom_point(data = merged_df, aes(x = location, y = p_value, color = color_condition), size = 3, alpha = merged_df$alpha) +
  scale_color_manual(values = c('#f6d992', '#f6a192', 'black', '#71c7ec')) +
  ylab("-log10 (p-value)") +
  xlab("Genome position") +
  geom_hline(linetype="dotted", 
             yintercept=-log10(crit_val)) +
  theme_bw() +
  ggtitle("Archipelago Plot") + 
  labs(subtitle = "Variant Set Association Test (blue)\nwith variant level contributions (alternating yellow/orange)\nand variant contribution highlights.") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.margin = margin(10, 10, 10, 10, "pt")) +
  guides(color = "none")

# Archipelago_plot_2_20_5k.pdf 6x4


# Version 3 - very large set -----
# We will plot points as layers for higherarchy
# Subset the data
variant_df <- merged_df[merged_df$metric == "variant",]
pathway_df <- merged_df[merged_df$metric == "pathway",]

ggplot() +
  geom_segment(data = line_data, 
               aes(x = location_pathway, xend = location, y = p_value_pathway, yend = p_value), alpha = 0.005) +
  # Layer 1: variant points
  geom_point(data = variant_df, aes(x = location, y = p_value, color = color_condition), size = 3, alpha = variant_df$alpha) +
  # Layer 2: pathway points
  geom_point(data = pathway_df, aes(x = location, y = p_value, color = color_condition), size = 3, alpha = pathway_df$alpha) +
  scale_color_manual(values = c('#f6d992', '#f6a192', 'black', '#71c7ec')) +
  ylab("-log10 (p-value)") +
  xlab("Genome position") +
  geom_hline(linetype="dotted", 
             yintercept=-log10(crit_val)) +
  theme_bw() +
  ggtitle("Archipelago Plot") + 
  labs(subtitle = "Variant Set Association Test (blue)\nwith variant level contributions (alternating yellow/orange)\nand variant contribution highlights.") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        plot.margin = margin(10, 10, 10, 10, "pt")) +
  guides(color = "none")

# Archipelago_plot_3_20_50k.pdf 6x4