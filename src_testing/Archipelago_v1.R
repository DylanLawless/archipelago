
# Archipelago Plot ----
library(ggplot2)
library(dplyr)
# Make synthetic data. 
# Pathway P-value is likely lower than variant P-values. 
# Set seed for reproducibility
set.seed(123)

# To do - replace coordinates with genomic Chromosome positions. 
# Assign rank per variant for simple plotting x-axis.

mulitpier <- 1
first <- 1
last <- 5000*mulitpier # Number of variants
MCLID_max <- 20
crit_val <- 0.0002
pval_max <- 0.001

# alternating_colors_max <- 500*mulitpier
# alternating_colors_min <- 250*mulitpier

# Sample df1: random numbers from a log-normal distribution
df1 <- data.frame(
  MCLID = seq(1, MCLID_max, by = 1),
  p_value = rlnorm(MCLID_max, meanlog = log(0.001), sdlog = 1)
)

df1$position <- NA
df1$location <- NA

# Sample df2: random from range
df2 <- data.frame(
  MCLID = rep(seq(1, MCLID_max, by = 1), each = 1),
  position = sample(1:1e6, last, replace = TRUE), 
  p_value = runif(last, 0.001, 0.01),
  chromosome = sample(1:22, last, replace = TRUE)  # Chromosome numbers 1-22
)

df2$location <- df2$position # location will be modified for the plot


# Calculate chromosome lengths based on min-max values in the data
chromosome_lengths_df <- df2 %>% 
  group_by(chromosome) %>%
  summarise(chromosome_length = max(location) + 1) %>%
  arrange(chromosome)

# Calculate cumulative chromosome lengths, 
# where we add +1 to each chromosome length except the first one (as it already starts from 1)
chromosome_lengths_df <- chromosome_lengths_df %>%
  mutate(cumulative_start = c(1, head(cumsum(chromosome_length) + 1, -1)))

# Join df2 with chromosome_lengths_df
df2 <- df2 %>%
  inner_join(chromosome_lengths_df, by = "chromosome")

# Calculate the new genome-wide location
df2 <- df2 %>%
  mutate(location = cumulative_start + location - 1)





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



# Normalize the genomic location
merged_df$location <- ifelse(merged_df$metric == "pathway", 
                             (merged_df$location - min(merged_df$location[merged_df$metric == "pathway"])) / 
                               (max(merged_df$location[merged_df$metric == "pathway"]) - min(merged_df$location[merged_df$metric == "pathway"])) * max(df2$location),
                             merged_df$location)


# Separate the data
pathway_points <- merged_df[merged_df$grouping == 1, ]
variant_points <- merged_df[merged_df$grouping > 1, ]

# Join the data
line_data <- variant_points %>%
  left_join(pathway_points, by = "MCLID", suffix = c("", "_pathway"))


# Add a new column for color groups
merged_df$color_group <- ifelse(merged_df$metric == "pathway", "Pathway", 
                                ifelse(merged_df$chromosome %% 2 == 0, "Chr_B", "Chr_A"))



# Create a dataframe with the mid-point of each chromosome
chromosome_ticks <- merged_df %>%
  filter(metric == "variant") %>%
  group_by(chromosome) %>%
  summarise(mid_point = mean(location))

# Plotting
ggplot() +
  geom_segment(data = line_data, 
               aes(x = location_pathway, xend = location, y = p_value_pathway, yend = p_value), alpha = 0.01) +
  geom_point(data = merged_df, aes(x = location, y = p_value, color = color_group), size = 3, alpha = ifelse(merged_df$metric == "pathway", 1, 0.5)) +
  scale_color_manual(values = c('#f6d992', '#f6a192', '#71c7ec')) +
  ylab("-log10 (p-value)") +
  xlab("Chromosome") +
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
  guides(color = "none")  + 
  scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$chromosome)

# Archipelago_plot_1_20_5k_v2.pdf 6x4



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
  xlab("Chromosome") +
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
  guides(color = "none") + 
  scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$chromosome)

# Archipelago_plot_2_20_5k_v2.pdf 6x4

