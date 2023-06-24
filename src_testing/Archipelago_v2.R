# To do - replace coordinates with genomic Chromosome positions. 
# Assign rank per variant for simple plotting x-axis.
# 
# mulitpier <- 1
# first <- 1
# last <- 5000*mulitpier # Number of variants
# MCLID_max <- 20
# crit_val <- 0.0002
# pval_max <- 0.001

# Archipelago Plot ----
library(ggplot2)
library(dplyr)

# Set seed for reproducibility
set.seed(123)

mulitpier <- 500
first <- 1
last <- mulitpier # Number of variants
MCLID_max <- (mulitpier/20)
crit_val <- .05/MCLID_max
pval_max <- 1
pval_min1 <- .05/8
pval_min2 <- .05/100

# alternating_colors_max <- 500*mulitpier
# alternating_colors_min <- 250*mulitpier

# Sample df1: random numbers from a log-normal distribution
df1 <- data.frame(
  MCLID = seq(1, MCLID_max, by = 1),
  p_value = rlnorm(MCLID_max, meanlog = log(pval_min1), sdlog = 1.4)
)

df1$position <- NA
df1$location <- NA

# Sample df2: random from range
df2 <- data.frame(
  MCLID = rep(seq(1, MCLID_max, by = 1), each = 1),
  position = sample(1:1e6, last, replace = TRUE), 
  p_value = runif(last, pval_min2, pval_max),
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

# Plott,ing
ggplot() +
  geom_segment(data = line_data, 
               aes(x = location_pathway, xend = location, y = p_value_pathway, yend = p_value), alpha = 0.01) +
  geom_point(data = merged_df, aes(x = location, y = p_value, color = color_group), size = 1, alpha = ifelse(merged_df$metric == "pathway", 1, 0.5)) +
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
  geom_point(data = merged_df, aes(x = location, y = p_value, color = color_condition), size = 1, alpha = merged_df$alpha) +
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




# Clearer condition_met layer ----

# Separate the points and lines that meet the condition
condition_met_points <- merged_df[merged_df$color_condition == "condition_met", ]
condition_met_lines <- line_data[line_data$MCLID %in% condition_met_points$MCLID, ]

# Plotting
ggplot() +
  # First layer: all other points and lines
  # geom_segment(data = line_data[!line_data$MCLID %in% condition_met_lines$MCLID, ], 
  #              aes(x = location_pathway, xend = location, y = p_value_pathway, yend = p_value), alpha = 0.01) +
  geom_point(data = merged_df[merged_df$color_condition != "condition_met", ], 
             aes(x = location, y = p_value, color = color_condition), size = 1#, alpha = merged_df$alpha
             ) +
  # Second layer: 'condition_met' points and lines
  geom_segment(data = condition_met_lines, 
               aes(x = location_pathway, xend = location, y = p_value_pathway, yend = p_value), alpha = 0.3) +
  geom_point(data = condition_met_points, 
             aes(x = location, y = p_value, color = color_condition), size = 1, alpha = 1) +
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

# Archipelago_plot_6_250_5k_v2.pdf






# ACAT plot ----
# Load necessary libraries
library(tidyverse)
library(ggplot2)

pval_min3 <- .05/10
pval_min4 <- .05/3
pval_min5 <- .05/2
pval_min6 <- .05/1
citial_individual <- 50

df3 <- data.frame(
	MCLID = rep(seq(1, 1, by = 1), each = 1),
	position = sample(1:1e6, citial_individual, replace = TRUE), 
	chromosome = sample(1:22, citial_individual, replace = TRUE), # Chromosome numbers 1-22
	gnomad_AF =  rlnorm(citial_individual, meanlog = log(pval_min3), sdlog = 1),
	CADD_phred = rlnorm(citial_individual, meanlog = log(pval_min4), sdlog = 1),
	REVEL =      rlnorm(citial_individual, meanlog = log(pval_min5), sdlog = 1),
	FATHMM =     rlnorm(citial_individual, meanlog = log(pval_min6), sdlog = 1)
)

df3$location <- df3$position # location will be modified for the plot

# Calculate chromosome lengths based on min-max values in the data
chromosome_lengths_df <- df3 %>% 
	group_by(chromosome) %>%
	summarise(chromosome_length = max(location) + 1) %>%
	arrange(chromosome)

# Calculate cumulative chromosome lengths, 
# where we add +1 to each chromosome length except the first one (as it already starts from 1)
chromosome_lengths_df <- chromosome_lengths_df %>%
	mutate(cumulative_start = c(1, head(cumsum(chromosome_length) + 1, -1)))

# Join df3 with chromosome_lengths_df
df3 <- df3 %>%
	inner_join(chromosome_lengths_df, by = "chromosome")

# Calculate the new genome-wide location
df3 <- df3 %>%
	mutate(location = cumulative_start + location - 1)


# Pivot data to long format
df_long <- df3 %>%
	pivot_longer(cols = 4:7,
					 names_to = "method",
					 values_to = "p_value")

# Create color groups based on method
df_long$P_annotation <- as.factor(df_long$method)

# Define critical p-value for significance
crit_val <- 0.05

# Create Archipelago plot
ggplot(df_long) +
	geom_point(aes(x = location, y = -log10(p_value), color = P_annotation), size = 1, alpha = 1) +
	# geom_hline(linetype="dotted", yintercept=-log10(crit_val)) +
	scale_color_manual(values = c('#f6d992', '#f6a192', '#71c7ec', '#92c6f6')) +
	ylab("-log10 (p-value)") +
	xlab("Position") +
	theme_bw() +
	ggtitle("Archipelago Plot") + 
	labs(subtitle = "ACAT variant level contributions.") +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) # +
	# guides(color = "none")



#

# Load necessary libraries
library(tidyverse)
library(ggplot2)

pval_min3 <- .05/10
pval_min4 <- .05/3
pval_min5 <- .05/2
pval_min6 <- .05/1
citial_individual <- 50
expand_axis <- 100000000

df3 <- data.frame(
	MCLID = rep(seq(1, 1, by = 1), each = 1),
	position = sample(1:1e6, citial_individual, replace = TRUE), 
	chromosome = sample(1:22, citial_individual, replace = TRUE), # Chromosome numbers 1-22
	gnomad_AF =  rlnorm(citial_individual, meanlog = log(pval_min3), sdlog = 1),
	CADD_phred = rlnorm(citial_individual, meanlog = log(pval_min4), sdlog = 1),
	REVEL =      rlnorm(citial_individual, meanlog = log(pval_min5), sdlog = 1),
	FATHMM =     rlnorm(citial_individual, meanlog = log(pval_min6), sdlog = 1)
)

df3$location <- df3$position # location will be modified for the plot

# Calculate chromosome lengths based on min-max values in the data
chromosome_lengths_df <- df3 %>% 
	group_by(chromosome) %>%
	summarise(chromosome_length = max(location) + 1) %>%
	arrange(chromosome)

# Calculate cumulative chromosome lengths, 
# where we add +1 to each chromosome length except the first one (as it already starts from 1)
chromosome_lengths_df <- chromosome_lengths_df %>%
	mutate(cumulative_start = c(1, head(cumsum(chromosome_length) + 1, -1)))

# Join df3 with chromosome_lengths_df
df3 <- df3 %>%
	inner_join(chromosome_lengths_df, by = "chromosome")

# Calculate the new genome-wide location
df3 <- df3 %>%
	mutate(location = cumulative_start + location - 1)



total_genome_length <- sum(chromosome_lengths_df$chromosome_length)

padding <-  total_genome_length - (chromosome_lengths_df$chromosome_length)

df3 <- df3 %>%
	mutate(padded_position = location + (chromosome - 1) * padding)

chromosome_ticks <- df_long %>%
	group_by(chromosome) %>%
	summarise(mid_point = mean(padded_position))


# Pivot data to long format
df_long <- df3 %>%
	pivot_longer(cols = 4:7,
					 names_to = "method",
					 values_to = "p_value")

# Create color groups based on method
df_long$P_annotation <- as.factor(df_long$method)

# Define critical p-value for significance
crit_val <- 0.05

# Calculate averages
avg_p_values <- df_long %>%
	group_by(P_annotation) %>%
	summarise(avg_p_value = mean(p_value))


# Create a dataframe with the mid-point of each chromosome
chromosome_ticks <- df_long %>%
	group_by(chromosome) %>%
	summarise(mid_point = mean(padded_position))

# # Create Archipelago plot
# ggplot(df_long) +
# 	geom_point(aes(x = location, y = -log10(p_value), color = P_annotation), size = 1, alpha = 1) +
# 	geom_hline(data = avg_p_values, aes(yintercept = -log10(avg_p_value), color = P_annotation), linetype="dashed") +
# 	scale_color_manual(values = c('#f6d992', '#f6a192', '#71c7ec', '#92c6f6')) +
# 	ylab("-log10 (p-value)") +
# 	xlab("Position") +
# 	labs(subtitle = "ACAT variant level contributions.") +
# 	# geom_text(data = avg_p_values, aes(x = Inf, y = -log10(avg_p_value), label = P_annotation, color = P_annotation), hjust = "inward", size = 5) +
# 	# geom_text(data = avg_p_values, aes(x = max(df_long$location) + expand_axis, y = -log10(avg_p_value), label = P_annotation, color = P_annotation), hjust = "inward", size = 5) +
# 	# coord_cartesian(xlim = c(min(df_long$location), max(df_long$location) + expand_axis)) +
# 	theme_bw() +
# 	# ggtitle("Archipelago Plot") + 
# 	labs(subtitle = "ACAT variant level contributions.") +
# 	theme(panel.grid.major = element_blank(),
# 			panel.grid.minor = element_blank(),
# 			panel.border = element_blank(),
# 			axis.line = element_line(color = "black"),
# 			plot.margin = margin(10, 10, 10, 10, "pt")) +
# 	guides(color = "none") +
# 	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$chromosome)
# 

# Create ACAT plot
ggplot(df_long) +
	geom_point(aes(x = padded_position, y = -log10(p_value), color = P_annotation), size = 1, alpha = 1) +
	scale_color_manual(values = c('#f6d992', '#f6a192', '#71c7ec', '#92c6f6')) +
	ylab("-log10 (p-value)") +
	xlab("Position") +
	labs(subtitle = "ACAT variant level contributions.") +
	theme_bw() +
	labs(subtitle = "ACAT variant level contributions.") +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) +
	guides(color = "none") +
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$chromosome) +
	geom_text(data = avg_p_values, aes(x = max(df_long$padded_position) + expand_axis*1.1, y = -log10(avg_p_value), label = P_annotation, color = P_annotation), hjust = "inward", size = 5) +
	geom_hline(data = avg_p_values, aes(yintercept = -log10(avg_p_value), color = P_annotation), linetype="dashed", alpha = 0.5) +
		coord_cartesian(xlim = c(min(df_long$padded_position), max(df_long$padded_position) + expand_axis))
# ACAT_50var_4anno_v1.pdf 6x4

ggplot(df_long) +
	geom_point(aes(x = padded_position, y = -log10(p_value), color = P_annotation), size = 1, alpha = 1) +
	geom_violin(scale = "area", aes(x = padded_position, y = -log10(p_value), group = padded_position), alpha = 0.1) +
	scale_color_manual(values = c('#f6d992', '#f6a192', '#71c7ec', '#92c6f6')) +
	ylab("-log10 (p-value)") +
	xlab("Position") +
	labs(subtitle = "ACAT variant level contributions.") +
	theme_bw() +
	labs(subtitle = "ACAT variant level contributions.") +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) +
	guides(color = "none") +
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$chromosome) +
	geom_text(data = avg_p_values, aes(x = max(df_long$padded_position) + expand_axis*1.1, y = -log10(avg_p_value), label = P_annotation, color = P_annotation), hjust = "inward", size = 5) +
	geom_hline(data = avg_p_values, aes(yintercept = -log10(avg_p_value), color = P_annotation), linetype="dashed", alpha = 0.5) +
	coord_cartesian(xlim = c(min(df_long$padded_position), max(df_long$padded_position) + expand_axis))
# ACAT_50var_4anno_v2.pdf 6x4