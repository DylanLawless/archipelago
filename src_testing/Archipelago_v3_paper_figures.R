library(ggplot2)
library(dplyr)

# Set data ----
# Set seed for reproducibility
set.seed(123)

mulitpier <- 5000
first <- 1
last <- mulitpier # Number of variants
MCLID_max <- (mulitpier/20)
crit_val <- .05/MCLID_max
pval_max <- 1
pval_min1 <- .05/10
pval_min2 <- .05/100

# Sample df1: random numbers from a log-normal distribution
df1 <- data.frame(
	MCLID = seq(1, MCLID_max, by = 1),
	P = rlnorm(MCLID_max, meanlog = log(pval_min1), sdlog = 1.6)
)

head(df1$P)
df1$P <- -log10(df1$P)
head(df1$P)

# Sample df2: random from range
d <- data.frame(
	MCLID = rep(seq(1, MCLID_max, by = 1), each = 1),
	BP = sample(1:1e6, last, replace = TRUE), 
	P = runif(last, pval_min2, pval_max),
	CHR = sample(1:22, last, replace = TRUE)  # Chromosome numbers 1-22
)

d$SNP <- d$BP # this can be modified for multiallelic sites

d$P <- -log10(d$P)
d <- d[order(d$CHR, d$BP), ]
d$index=NA
ind = 0
for (i in unique(d$CHR)){
	ind = ind + 1
	d[d$CHR==i,]$index = ind
}
d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$BP,d$CHR,length))  # replcace the for loop of line 92-96 to improve efficiency

nchr = length(unique(d$CHR))
d$pos=NA

# The following code block replaces my old method - to more closely match qqman - to make it easier to integrate in future.
# if (nchr==1) { ## For a single chromosome
# 	## Uncomment the next two linex to plot single chr results in Mb
# 	#options(scipen=999)
# 	#d$pos=d$BP/1e6
# 	# d$pos=d$BP
# 	# #  ticks=floor(length(d$pos))/2+1          ## unused, from code line: 169
# 	# xlabel = paste('Chromosome',unique(d$CHR),'position')
# 	#  labs = ticks          ## unused, from code line: 169
# } else { ## For multiple chromosomes
	lastbase=0
	ticks=NULL
	for (i in unique(d$index)) {
		if (i==1) {
			d[d$index==i, ]$pos=d[d$index==i, ]$BP
		} else {
			## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced. 
			lastbase = lastbase +max(d[d$index==(i-1),"BP"])   # replace line 128
			d[d$index == i,"BP"] = d[d$index == i,"BP"]-min(d[d$index==i,"BP"]) +1
			d[d$index == i, "pos"] = d[d$index == i,"BP"] + lastbase    # replace line 129
			# lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
			# d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
			
		}
		# Old way: assumes SNPs evenly distributed
		# ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
		# New way: doesn't make that assumption
		# ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)  # see line 136, to reduce the burden of for loop 
	}
	ticks <-tapply(d$pos,d$index,quantile,probs=0.5)   # replace line 135
	xlabel = 'Chromosome'
	labs <- unique(d$CHR)
# }

# Merging df1 and df2
merged_df <- bind_rows(df1, d)

# Add a new column to specify variant_set points
merged_df$metric <- ifelse(is.na(merged_df$BP), "variant_set", "variant")

# Calculate the total of locations per MCLID
pos_sum <- aggregate(merged_df$pos, by=list(merged_df$MCLID), FUN=sum, na.rm=TRUE)

# Rank MCLID by the sum of locations
pos_sum$rank <- rank(pos_sum$x) 

# Calculate the average location within each MCLID group from the original dataframe
average_pos <- aggregate(merged_df$pos, by=list(merged_df$MCLID), FUN=mean, na.rm=TRUE)

# Replace NA locations with the average location of their MCLID group
merged_df$pos[is.na(merged_df$pos)] <- average_pos$x[match(merged_df$MCLID[is.na(merged_df$pos)], average_pos$Group.1)]

# Create a new grouping variable
merged_df$grouping <- with(merged_df, ave(pos, MCLID, FUN = function(x) cumsum(!is.na(x))))



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
							  variant_set_sorted$pos[match(MCLID, variant_set_sorted$MCLID)],
							  pos))

# The previous 4 code chunks constructs a natural spread of positions for "variant_set". However, their positions are distributioned slightly atrificially to give a more even distribution in their true odeder. This is ideal for a dense network. If you have [1] a less dense network, [2] don't mind high-density clustiner at the genome center, [3] want to see the less dispersed distribution, just uncomment the next chuck to override the previous distribution version. 

# merged_df$pos <- ifelse(merged_df$metric == "variant_set",
# 								(merged_df$pos - min(merged_df$pos[merged_df$metric == "variant_set"])) /
# 									(max(merged_df$pos[merged_df$metric == "variant_set"]) - min(merged_df$pos[merged_df$metric == "variant_set"])) * max(d$pos),
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
line_data <- left_join(variant_data, variant_set_data, by = "MCLID", suffix = c("_variant", "_variant_set"))

# Raw plot ----
# 
raw <- merged_df %>% 
  filter(metric == "variant_set")
raw$rank <- rank(raw$P)

raw %>%
  ggplot() +
  geom_point(aes(x = rank, y = P), color='#27afea', size = 1) +
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
ggsave("../images/raw_vsat_plot_250set_v1.pdf", width = 6, height = 4)

# Plot 1 ----
ggplot() +
	geom_segment(data = line_data, 
					 aes(x = pos_variant_set, xend = pos_variant, y = P_variant_set, yend = P_variant), alpha = 0.02) +
	geom_point(data = merged_df, aes(x = pos, y = P, color = color_group), size = 1, alpha = ifelse(merged_df$metric == "variant_set", 1, 0.5)) +
	scale_color_manual(values = c('#f6d992', '#f6a192', '#27afea')) +
	ylab("-log10 (p-value)") +
	xlab("Chromosome") +
	geom_hline(linetype="dotted", 
				  yintercept=-log10(crit_val)) +
	theme_bw() +
	ggtitle("Archipelago Plot") + 
	labs(subtitle = "Variant Set Association Test (blue)\nwith individual variant contributions (alternating yellow/orange).") +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) +
	guides(color = "none")  + 
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR)

ggsave("../images/Archipelago_plot_5kvar_250set_v3.pdf", width = 6, height = 4)

# Plot 2 ----
# Highlight individual variant contributions 
# Create a new variable 'color_condition' that checks the condition
merged_df <- 
	merged_df %>%
	group_by(MCLID) %>%
	mutate(color_condition = 
			 	ifelse(any(P > -log10(crit_val)) & metric == "variant", 
			 			 "condition_met", 
			 			 color_group))

# Change the alpha variable accordingly
merged_df$alpha <- ifelse(merged_df$color_condition == "condition_met", 1, ifelse(merged_df$metric == "variant_set", 1, 0.3))

# Plotting
ggplot() +
	geom_segment(data = line_data, 
					 aes(x = pos_variant_set, xend = pos_variant, y = P_variant_set, yend = P_variant), alpha = 0.01) +
	geom_point(data = merged_df, aes(x = pos, y = P, color = color_condition), size = 1, alpha = merged_df$alpha) +
	scale_color_manual(values = c('#f6d992', '#f6a192', 'black', '#27afea')) +
	ylab("-log10 (p-value)") +
	xlab("Chromosome") +
	geom_hline(linetype="dotted", 
				  yintercept=-log10(crit_val)) +
	theme_bw() +
	ggtitle("Archipelago Plot") + 
	labs(subtitle = "Variant Set Association Test (blue)\nwith individual variant contributions (alternating yellow/orange)\nand contributions for significant VSAT (black)") +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) +
	guides(color = "none") + 
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR)

ggsave("../images/Archipelago_plot_5kvar_250set_v2.pdf", width = 6, height = 4)

# Plot 4 ----
# Clearer condition_met layer
# Separate the points and lines that meet the condition
condition_met_points <- merged_df[merged_df$color_condition == "condition_met", ]
condition_met_lines <- line_data[line_data$MCLID %in% condition_met_points$MCLID, ]

# Plotting
p_arch <- 
	ggplot() +
	geom_point(data = merged_df[merged_df$color_condition != "condition_met", ], 
				  aes(x = pos, y = P, color = color_condition), size = 1) +
	geom_segment(data = condition_met_lines, 
					 aes(x = pos_variant_set, xend = pos_variant, y = P_variant_set, yend = P_variant), alpha = 0.3) +
	geom_point(data = condition_met_points, 
				  aes(x = pos, y = P, color = color_condition), size = 1, alpha = 1) +
	scale_color_manual(values = c('#f6d992', '#f6a192', 'black', '#27afea')) +
	ylab("-log10 (p-value)") +
	xlab("Chromosome") +
	geom_hline(linetype="dotted", 
				  yintercept=-log10(crit_val)) +
	theme_bw() +
	ggtitle("Archipelago Plot") + 
	labs(subtitle = "Variant Set Association Test (blue)\nwith individual variant contributions (alternating yellow/orange)\nand contributions for significant VSAT (black)") +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) +
	guides(color = "none") +
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR)

p_arch
ggsave("../images/Archipelago_plot_5kvar_250set_v1.pdf", width = 6, height = 4)



# Figure legend ----
p_arch_leg <- 
	ggplot() +
	geom_point(data = merged_df[merged_df$color_condition != "condition_met", ], 
				  aes(x = pos, y = P, color = color_condition), size = 1) +
	geom_segment(data = condition_met_lines, 
					 aes(x = pos_variant_set, xend = pos_variant, y = P_variant_set, yend = P_variant), alpha = 0.3) +
	geom_point(data = condition_met_points, 
				  aes(x = pos, y = P, color = color_condition), size = 1, alpha = 1) +
	scale_color_manual(values = c('#f6d992', '#f6a192', 'black', '#27afea'), labels= c("Chr: Individual p-val", "Chr: Individual p-val" ,"Individual variant\nfrom enriched VSAT","VSAT p-val")) +
	ylab("-log10 (p-value)") +
	xlab("Chromosome") +
	geom_hline(linetype="dotted", 
				  yintercept=-log10(crit_val)) +
	theme_bw() +
	ggtitle("Archipelago Plot") + 
	labs(subtitle = "Variant Set Association Test (blue)\nwith individual variant contributions (alternating yellow/orange)\nand contributions for significant VSAT (black)") +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) +
	# guides(color = "none") + 
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR)

p_arch_leg
ggsave("../images/Archipelago_plot_5kvar_250set_legend_v1.pdf", width = 8, height = 4)


# For more complex geom_segments we may need multi-color scale

# Plot 5 ----
ggplot() +
	geom_point(data = merged_df[merged_df$color_condition != "condition_met", ], 
				  aes(x = pos, y = P, color = color_condition), size = 1) +
	geom_segment(data = condition_met_lines, 
					 aes(x = pos_variant_set, xend = pos_variant, y = P_variant_set, yend = P_variant, color = as.factor(MCLID)), alpha = 0.6) +
	geom_point(data = condition_met_points, 
				  aes(x = pos, y = P, color = color_condition), size = 1, alpha = 1) +
	scale_color_manual(values = c('#0224a9', '#a602a9', '#f6d992', '#f6a192', 'black', '#27afea' )) +
	ylab("-log10 (p-value)") +
	xlab("Chromosome") +
	geom_hline(linetype="dotted", 
				  yintercept=-log10(crit_val)) +
	theme_bw() +
	ggtitle("Archipelago Plot") + 
	labs(subtitle = "Variant Set Association Test (blue)\nwith individual variant contributions (alternating yellow/orange)\nand contributions for significant VSAT (black)") +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) +
	guides(color = "none") + 
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR)
ggsave("../images/Archipelago_plot_5kvar_250set_seg_v1.pdf", width = 6, height = 4)


# ACAT plot ----
# Load necessary libraries
library(tidyverse)
set.seed(123)
# We will take our top hits from the previous example
df3 <- merged_df %>% filter(color_condition == "condition_met")
names(df3)

# And now generate some example p-value for each annotation layer
pval_min3 <- .05/10
pval_min4 <- .05/3
pval_min5 <- .05/2
pval_min6 <- .05/1

# Get the number of rows in the data frame
num_rows <- nrow(df3)

# Add the new p-values for each row
df3 <- df3 %>% 
	ungroup() %>%
	mutate(
		gnomad_AF =  rlnorm(num_rows, meanlog = log(pval_min3), sdlog = 1),
		CADD_phred = rlnorm(num_rows, meanlog = log(pval_min4), sdlog = 1),
		REVEL =      rlnorm(num_rows, meanlog = log(pval_min5), sdlog = 1),
		FATHMM =     rlnorm(num_rows, meanlog = log(pval_min6), sdlog = 1)
	)

# Pivot data to long format
df_long <- df3 %>%
	pivot_longer(cols = gnomad_AF:FATHMM,
					 names_to = "method",
					 values_to = "p_value")

# Create color groups based on method
df_long$P_annotation <- as.factor(df_long$method)

# Calculate averages
avg_p_values <- df_long %>%
	group_by(P_annotation) %>%
	summarise(avg_p_value = mean(p_value))

# Set the levels of P_annotation according to the order of avg_p_value
df_long$P_annotation <- factor(df_long$P_annotation, levels = avg_p_values$P_annotation[order(avg_p_values$avg_p_value)])

# Create ACAT plot
p_ACAT <- ggplot(df_long) +
	geom_point(aes(x = pos, y = -log10(p_value), color = P_annotation), size = 1, alpha = 1) +
	scale_color_manual(values =  c('#f6d992', '#f6a192', '#d192f6', '#6469ff'), name = "Average p-value\nannotation") +
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR) +
	geom_hline(data = avg_p_values, aes(yintercept = -log10(avg_p_value), color = P_annotation), linetype="dashed", alpha = 1) +
	labs(subtitle = "ACAT variant level contributions.") +
	ylab("-log10 (p-value)") +
	xlab("Position") +
	theme_bw() +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) 

p_ACAT
ggsave("../images/ACAT_5kvar_1set_40var_4anno_v1.pdf", width = 6, height = 4)


# Add a single-variant box/line
p_ACAT2 <- ggplot(df_long) +
	geom_point(aes(x = pos, y = -log10(p_value), color = P_annotation), size = 1, alpha = 1) +
	geom_boxplot(scale = "area", aes(x = pos, y = -log10(p_value), group = pos), alpha = 0.5) +
	scale_color_manual(values = c('#f6d992', '#f6a192', '#d192f6', '#6469ff'), name = "Average p-value\nannotation") +
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR) +
	geom_hline(data = avg_p_values, aes(yintercept = -log10(avg_p_value), color = P_annotation), linetype="dashed", alpha = 1) +
	labs(subtitle = "ACAT variant level contributions.") +
	ylab("-log10 (p-value)") +
	xlab("Position") +
	theme_bw() +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) 

p_ACAT2
ggsave("../images/ACAT_5kvar_1set_40var_4anno_v2.pdf", width = 6, height = 4)


# Save the arranged plots to a PDF
# ggsave("../images/ACAT_5kvar_1set_40var_4anno_v3.pdf")

# install.packages("patchwork")  # Uncomment if package is not installed yet
library(patchwork)

p_arch / p_ACAT / p_ACAT2 + plot_layout(guides = "collect")
ggsave("../images/stack_v1.pdf", width = 6, height = 9)



# Small version ----


library(ggplot2)
library(dplyr)

# Set data ----
# Set seed for reproducibility
set.seed(123)

mulitpier <- 500
first <- 1
last <- mulitpier # Number of variants
MCLID_max <- (mulitpier/20)
crit_val <- .05/MCLID_max
pval_max <- 1
pval_min1 <- .05/7
pval_min2 <- .05/10

# Sample df1: random numbers from a log-normal distribution
df1 <- data.frame(
	MCLID = seq(1, MCLID_max, by = 1),
	P = rlnorm(MCLID_max, meanlog = log(pval_min1), sdlog = 1)
)

head(df1$P)
df1$P <- -log10(df1$P)
head(df1$P)

# Sample df2: random from range
d <- data.frame(
	MCLID = rep(seq(1, MCLID_max, by = 1), each = 1),
	BP = sample(1:1e6, last, replace = TRUE), 
	P = runif(last, pval_min2, pval_max),
	CHR = sample(1:22, last, replace = TRUE)  # Chromosome numbers 1-22
)


d$SNP <- d$BP # this can be modified for multiallelic sites

d$P <- -log10(d$P)
d <- d[order(d$CHR, d$BP), ]
d$index=NA
ind = 0
for (i in unique(d$CHR)){
	ind = ind + 1
	d[d$CHR==i,]$index = ind
}
d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$BP,d$CHR,length))  # replcace the for loop of line 92-96 to improve efficiency

nchr = length(unique(d$CHR))
d$pos=NA

# The following code block replaces my old method - to more closely match qqman - to make it easier to integrate in future.
# if (nchr==1) { ## For a single chromosome
# 	## Uncomment the next two linex to plot single chr results in Mb
# 	#options(scipen=999)
# 	#d$pos=d$BP/1e6
# 	# d$pos=d$BP
# 	# #  ticks=floor(length(d$pos))/2+1          ## unused, from code line: 169
# 	# xlabel = paste('Chromosome',unique(d$CHR),'position')
# 	#  labs = ticks          ## unused, from code line: 169
# } else { ## For multiple chromosomes
lastbase=0
ticks=NULL
for (i in unique(d$index)) {
	if (i==1) {
		d[d$index==i, ]$pos=d[d$index==i, ]$BP
	} else {
		## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced. 
		lastbase = lastbase +max(d[d$index==(i-1),"BP"])   # replace line 128
		d[d$index == i,"BP"] = d[d$index == i,"BP"]-min(d[d$index==i,"BP"]) +1
		d[d$index == i, "pos"] = d[d$index == i,"BP"] + lastbase    # replace line 129
		# lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
		# d[d$index==i, ]$pos=d[d$index==i, ]$BP+lastbase
		
	}
	# Old way: assumes SNPs evenly distributed
	# ticks=c(ticks, d[d$index==i, ]$pos[floor(length(d[d$index==i, ]$pos)/2)+1])
	# New way: doesn't make that assumption
	# ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)  # see line 136, to reduce the burden of for loop 
}
ticks <-tapply(d$pos,d$index,quantile,probs=0.5)   # replace line 135
xlabel = 'Chromosome'
labs <- unique(d$CHR)
# }

# Merging df1 and df2
merged_df <- bind_rows(df1, d)

# Add a new column to specify variant_set points
merged_df$metric <- ifelse(is.na(merged_df$BP), "variant_set", "variant")

# Calculate the total of locations per MCLID
pos_sum <- aggregate(merged_df$pos, by=list(merged_df$MCLID), FUN=sum, na.rm=TRUE)

# Rank MCLID by the sum of locations
pos_sum$rank <- rank(pos_sum$x) 

# Calculate the average location within each MCLID group from the original dataframe
average_pos <- aggregate(merged_df$pos, by=list(merged_df$MCLID), FUN=mean, na.rm=TRUE)

# Replace NA locations with the average location of their MCLID group
merged_df$pos[is.na(merged_df$pos)] <- average_pos$x[match(merged_df$MCLID[is.na(merged_df$pos)], average_pos$Group.1)]

# Create a new grouping variable
merged_df$grouping <- with(merged_df, ave(pos, MCLID, FUN = function(x) cumsum(!is.na(x))))



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
							  variant_set_sorted$pos[match(MCLID, variant_set_sorted$MCLID)],
							  pos))

# The previous 4 code chunks constructs a natural spread of positions for "variant_set". However, their positions are distributioned slightly atrificially to give a more even distribution in their true odeder. This is ideal for a dense network. If you have [1] a less dense network, [2] don't mind high-density clustiner at the genome center, [3] want to see the less dispersed distribution, just uncomment the next chuck to override the previous distribution version. 

# merged_df$pos <- ifelse(merged_df$metric == "variant_set",
# 								(merged_df$pos - min(merged_df$pos[merged_df$metric == "variant_set"])) /
# 									(max(merged_df$pos[merged_df$metric == "variant_set"]) - min(merged_df$pos[merged_df$metric == "variant_set"])) * max(d$pos),
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
line_data <- left_join(variant_data, variant_set_data, by = "MCLID", suffix = c("_variant", "_variant_set"))

# Plot 1 ----
ggplot() +
	geom_segment(data = line_data, 
					 aes(x = pos_variant_set, xend = pos_variant, y = P_variant_set, yend = P_variant), alpha = 0.02) +
	geom_point(data = merged_df, aes(x = pos, y = P, color = color_group), size = 1, alpha = ifelse(merged_df$metric == "variant_set", 1, 0.5)) +
	scale_color_manual(values = c('#f6d992', '#f6a192', '#27afea')) +
	ylab("-log10 (p-value)") +
	xlab("Chromosome") +
	geom_hline(linetype="dotted", 
				  yintercept=-log10(crit_val)) +
	theme_bw() +
	ggtitle("Archipelago Plot") + 
	labs(subtitle = "Variant Set Association Test (blue)\nwith individual variant contributions (alternating yellow/orange).") +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) +
	guides(color = "none")  + 
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR)

ggsave("../images/Archipelago_plot_500var_25set_small_v3.pdf", width = 6, height = 4)

# Plot 2 ----
# Highlight individual variant contributions 
# Create a new variable 'color_condition' that checks the condition
merged_df <- 
	merged_df %>%
	group_by(MCLID) %>%
	mutate(color_condition = 
			 	ifelse(any(P > -log10(crit_val)) & metric == "variant", 
			 			 "condition_met", 
			 			 color_group))

# Change the alpha variable accordingly
merged_df$alpha <- ifelse(merged_df$color_condition == "condition_met", 1, ifelse(merged_df$metric == "variant_set", 1, 0.3))

# Plotting
ggplot() +
	geom_segment(data = line_data, 
					 aes(x = pos_variant_set, xend = pos_variant, y = P_variant_set, yend = P_variant), alpha = 0.01) +
	geom_point(data = merged_df, aes(x = pos, y = P, color = color_condition), size = 1, alpha = merged_df$alpha) +
	scale_color_manual(values = c('#f6d992', '#f6a192', 'black', '#27afea')) +
	ylab("-log10 (p-value)") +
	xlab("Chromosome") +
	geom_hline(linetype="dotted", 
				  yintercept=-log10(crit_val)) +
	theme_bw() +
	ggtitle("Archipelago Plot") + 
	labs(subtitle = "Variant Set Association Test (blue)\nwith individual variant contributions (alternating yellow/orange)\nand contributions for significant VSAT (black)") +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) +
	guides(color = "none") + 
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR)

ggsave("../images/Archipelago_plot_500var_25set_small_v2.pdf", width = 6, height = 4)

# Plot 4 ----
# Clearer condition_met layer
# Separate the points and lines that meet the condition
condition_met_points <- merged_df[merged_df$color_condition == "condition_met", ]
condition_met_lines <- line_data[line_data$MCLID %in% condition_met_points$MCLID, ]

# Plotting
p_arch <- 
	ggplot() +
	geom_point(data = merged_df[merged_df$color_condition != "condition_met", ], 
				  aes(x = pos, y = P, color = color_condition), size = 1) +
	geom_segment(data = condition_met_lines, 
					 aes(x = pos_variant_set, xend = pos_variant, y = P_variant_set, yend = P_variant), alpha = 0.3) +
	geom_point(data = condition_met_points, 
				  aes(x = pos, y = P, color = color_condition), size = 1, alpha = 1) +
	scale_color_manual(values = c('#f6d992', '#f6a192', 'black', '#27afea')) +
	ylab("-log10 (p-value)") +
	xlab("Chromosome") +
	geom_hline(linetype="dotted", 
				  yintercept=-log10(crit_val)) +
	theme_bw() +
	ggtitle("Archipelago Plot") + 
	labs(subtitle = "Variant Set Association Test (blue)\nwith individual variant contributions (alternating yellow/orange)\nand contributions for significant VSAT (black)") +
	theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			axis.line = element_line(color = "black"),
			plot.margin = margin(10, 10, 10, 10, "pt")) +
	guides(color = "none") + 
	scale_x_continuous(breaks = chromosome_ticks$mid_point, labels = chromosome_ticks$CHR)

p_arch
ggsave("../images/Archipelago_plot_500var_25set_small_v1.pdf", width = 6, height = 4)


