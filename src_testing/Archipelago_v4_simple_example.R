library(ggplot2)
library(dplyr)

# Set data ----
# Set seed for reproducibility
set.seed(123)


mulitpier <- 5000
first <- 1
last <- mulitpier # Number of variants
set_ID_max <- (mulitpier/20)
crit_val <- .05/set_ID_max
pval_max <- 1
pval_min1 <- .05/10
pval_min2 <- .05/100

# Sample df1: random numbers from a log-normal distribution
df1 <- data.frame(
	set_ID = seq(1, set_ID_max, by = 1),
	P = rlnorm(set_ID_max, meanlog = log(pval_min1), sdlog = 1.6)
)

df1$P <- -log10(df1$P)

# Sample df2: random from range
d <- data.frame(
	set_ID = rep(seq(1, set_ID_max, by = 1), each = 1),
	BP = sample(1:1e6, last, replace = TRUE), 
	P = runif(last, pval_min2, pval_max),
	CHR = sample(1:22, last, replace = TRUE)  # Chromosome numbers 1-22
)

d$SNP <- d$BP # this would be required for multiallelic sites

# mulitpier <- number_of_variants
# set_ID_max <- number_of_sets

d$P <- -log10(d$P)
d <- d[order(d$CHR, d$BP), ]
d$index=NA
ind = 0
for (i in unique(d$CHR)){
	ind = ind + 1
	d[d$CHR==i,]$index = ind
}
d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$BP,d$CHR,length))  

nchr = length(unique(d$CHR))
d$pos=NA

# The following code block replaces my old method - to more closely match qqman - to make it easier to integrate in future.
lastbase=0
ticks=NULL
for (i in unique(d$index)) {
	if (i==1) {
		d[d$index==i, ]$pos=d[d$index==i, ]$BP
	} else {
		lastbase = lastbase +max(d[d$index==(i-1),"BP"])  
		d[d$index == i,"BP"] = d[d$index == i,"BP"]-min(d[d$index==i,"BP"]) +1
		d[d$index == i, "pos"] = d[d$index == i,"BP"] + lastbase 
		
	}
}
ticks <-tapply(d$pos,d$index,quantile,probs=0.5) 
xlabel = 'Chromosome'
labs <- unique(d$CHR)

# Merging df1 and df2
merged_df <- bind_rows(df1, d)

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
line_data <- left_join(variant_data, variant_set_data, by = "set_ID", suffix = c("_variant", "_variant_set"))

# Raw plot ----
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
# ggsave("../images/raw_vsat_plot_250set_v1.pdf", width = 6, height = 4)

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
# ggsave("../images/Archipelago_plot_5kvar_250set_legend_v1.pdf", width = 8, height = 4)
