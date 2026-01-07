# Load the packages to verify installation
setwd('/home/veo/pantanos_stats')
library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)
library(viridis)
library(dendextend)
library(tidyr)
library(dplyr)
library(patchwork)
library(ggpubr)
library(dunn.test)
library(reshape2)

tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam-v2.tsv", header=T, row.names=1, check.names=F, sep="\t"))
print("Successfully read taxonomy file")
sample_info_tab <- read.csv("pantano_metadata.csv", header=T, row.names=1, check.names=F)
print("Successfully read metadata file")
pre_count_tab <- read.table("ASVs_counts-no-contam-v2.tsv", header=T, row.names=1, check.names=F, sep="\t")[ , -c(7:9)]
print("Successfully read counts file")

colnames(pre_count_tab) <- gsub("reads/", "", colnames(pre_count_tab))
pre_count_tab <- pre_count_tab + 1 # Because of the zeros
deseq_counts <- DESeqDataSetFromMatrix(pre_count_tab, colData = sample_info_tab, design = ~mes) # Cambiarlo para que se introduzca desde param_file
deseq_counts <- estimateSizeFactors(deseq_counts)
count_tab <- counts(deseq_counts, normalized = TRUE)

# Process taxonomy table to replace NAs with "Unclassified"
tax_tab_processed <- tax_tab
tax_tab_processed[is.na(tax_tab_processed)] <- "Unclassified"

# Calculate relative abundance for each sample
rel_abund <- sweep(count_tab, 2, colSums(count_tab), "/") * 100

# Create a data frame with taxonomic and abundance information
taxa_abund_df <- data.frame(
  ASV = rownames(tax_tab_processed),
  Class = tax_tab_processed[, "class"],
  stringsAsFactors = FALSE
)

# Replace NA with "Unclassified" in the dataframe
taxa_abund_df$Class[is.na(taxa_abund_df$Class)] <- "Unclassified"

# Add abundance data
for (sample in colnames(rel_abund)) {
  taxa_abund_df[[sample]] <- rel_abund[taxa_abund_df$ASV, sample]
}

# Calculate total abundance per class across all samples
class_totals <- aggregate(. ~ Class, data = taxa_abund_df[, -1], sum)
class_totals$total <- rowSums(class_totals[, -1])
class_totals <- class_totals[order(class_totals$total, decreasing = TRUE), ]

# Get top 10 classes
top_classes <- class_totals$Class[1:min(10, nrow(class_totals))]

# Reshape data for plotting
plot_data <- taxa_abund_df
plot_data$Major_Taxa <- ifelse(plot_data$Class %in% top_classes, plot_data$Class, "Other")

# Aggregate by Major_Taxa
plot_data_agg <- aggregate(. ~ Major_Taxa, data = plot_data[, -c(1, 2)], sum)

# Reshape to long format for ggplot
plot_data_long <- melt(plot_data_agg, id.vars = "Major_Taxa", variable.name = "Sample", value.name = "Proportion")

# Create a custom pastel/dull color palette
pastel_colors <- c(
  "#8FB9AA", "#F2D096", "#F2A679", "#D95D39", "#35524A",
  "#7C9EB2", "#C9ADA7", "#9A8C98", "#797596", "#4A4E69",
  "#A9927D"  # For "Other" category
)

# Order Major_Taxa to have "Other" at the bottom of the stack
plot_data_long$Major_Taxa <- factor(plot_data_long$Major_Taxa, 
                                    levels = c(top_classes, "Other"))

# Create the stacked barplot
ggplot(plot_data_long, aes(x = Sample, y = Proportion, fill = Major_Taxa)) +
  geom_bar(width = 0.6, stat = "identity") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.4, hjust = 1),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_fill_manual(values = pastel_colors) +
  labs(
    x = "Sample",
    y = "% of 16S rRNA gene copies recovered",
    title = "Top 10 Bacterial and Archaeal Classes"
  )

# Define lake_order
lake_order <- sample_info_tab %>% 
  arrange(lago) %>% 
  rownames()

#----- CLASS ------
# 1. Prepare the data
# Create combined taxonomy-abundance dataframe
taxa_abund_df <- data.frame(
  ASV = rownames(tax_tab_processed),
  Domain = tax_tab_processed[, "domain"],
  Class = tax_tab_processed[, "class"],
  stringsAsFactors = FALSE
) %>% 
  mutate(Class = ifelse(is.na(Class), "Unclassified", Class),
         Domain = ifelse(is.na(Domain), "Unassigned", Domain))

# Add relative abundances
for(sample in colnames(rel_abund)) {
  taxa_abund_df[[sample]] <- rel_abund[taxa_abund_df$ASV, sample]
}

# 2. Split into bacterial and archaeal datasets
bacterial_df <- taxa_abund_df %>% filter(Domain == "Bacteria")
archaeal_df <- taxa_abund_df %>% filter(Domain == "Archaea")

# 3. Process bacterial data
bacterial_agg <- bacterial_df %>% 
  group_by(Class) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(total))

bacterial_top10 <- bacterial_agg$Class[1:min(10, nrow(bacterial_agg))]

bacterial_plot_data <- bacterial_df %>% 
  mutate(Major_Taxa = ifelse(Class %in% bacterial_top10, Class, "Other")) %>% 
  group_by(Major_Taxa) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  # Normalize to 100%
  mutate(across(where(is.numeric), ~ .x / sum(.x) * 100)) %>% 
  melt(id.vars = "Major_Taxa", variable.name = "Sample", value.name = "Proportion") %>% 
  merge(sample_info_tab, by.x = "Sample", by.y = "row.names")

# 4. Process archaeal data
archaeal_agg <- archaeal_df %>% 
  group_by(Class) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(total))

archaeal_top10 <- archaeal_agg$Class[1:min(10, nrow(archaeal_agg))]

archaeal_plot_data <- archaeal_df %>% 
  mutate(Major_Taxa = ifelse(Class %in% archaeal_top10, Class, "Other")) %>% 
  group_by(Major_Taxa) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  # Normalize to 100%
  mutate(across(where(is.numeric), ~ .x / sum(.x) * 100)) %>% 
  melt(id.vars = "Major_Taxa", variable.name = "Sample", value.name = "Proportion") %>% 
  merge(sample_info_tab, by.x = "Sample", by.y = "row.names")

# Create color palettes
bacterial_colors <- c(
  "#8FB9AA", "#F2D096", "#F2A679", "#D95D39", "#35524A",
  "#7C9EB2", "#C9ADA7", "#9A8C98", "#797596", "#4A4E69",
  "#A9927D"  # For "Other" category
)

archaeal_colors <- c(
  "#E6B89C", "#FECEA8", "#FF847C", "#E84A5F", "#2A363B",
  "#F7C59F", "#D6A2AD", "#B07BAC", "#774C60", "#CFBACE",
  "#7D7ABC"  # For "Other" category
)

# 6. Generate the plots
bacterial_plot <- ggplot(bacterial_plot_data, aes(x = Sample, y = Proportion, fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = bacterial_colors) +
  labs(title = "Top 10 Bacterial Classes",
       y = "Relative Abundance (%)",
       x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

archaeal_plot <- ggplot(archaeal_plot_data, aes(x = Sample, y = Proportion, fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = archaeal_colors) +
  labs(title = "Top 10 Archaeal Classes",
       y = "Relative Abundance (%)",
       x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

# 7. Arrange and save plots
combined_plots <- bacterial_plot / archaeal_plot + 
  plot_layout(heights = c(1, 1))

#ggsave("relative_abundance_plots.png", combined_plots, width = 12, height = 10, dpi = 300)

# 8. Show plots
combined_plots

#---- CLASS BY LAKES -------
# 1. Update the bacterial plot with boxes around lakes
bacterial_plot <- ggplot(bacterial_plot_data, 
                         aes(x = factor(Sample, levels = lake_order), 
                             y = Proportion, 
                             fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ lago, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = bacterial_colors) +
  labs(y = "Relative Abundance (%)",
       x = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.2, "lines"),
    panel.background = element_rect(color = "black", fill = NA),  # Add panel borders
    strip.background = element_rect(fill = "gray90", color = "black")  # Style facet labels
  )

# 2. Update the archaeal plot similarly
archaeal_plot <- ggplot(archaeal_plot_data, 
                        aes(x = factor(Sample, levels = lake_order), 
                            y = Proportion, 
                            fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ lago, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = archaeal_colors) +
  labs(y = "Relative Abundance (%)",
       x = "Sample") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.2, "lines"),
    panel.background = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "gray90", color = "black")
  )

# 3. Combine plots
combined_plots <- bacterial_plot / archaeal_plot + 
  plot_layout(heights = c(1, 1))

# 4. Save with wider aspect ratio
ggsave(file="Fig4_StackBars_By_Class.svg", plot=combined_plots, width=8, height=6)

# Show plots
combined_plots

#---- ORDER ------
# 1. Prepare the data
# Create combined taxonomy-abundance dataframe
taxa_abund_df <- data.frame(
  ASV = rownames(tax_tab_processed),
  Domain = tax_tab_processed[, "domain"],
  Order = tax_tab_processed[, "order"],
  stringsAsFactors = FALSE
) %>% 
  mutate(Order = ifelse(is.na(Order), "Unorderified", Order),
         Domain = ifelse(is.na(Domain), "Unassigned", Domain))

# Add relative abundances
for(sample in colnames(rel_abund)) {
  taxa_abund_df[[sample]] <- rel_abund[taxa_abund_df$ASV, sample]
}

# 2. Split into bacterial and archaeal datasets
bacterial_df <- taxa_abund_df %>% filter(Domain == "Bacteria")
archaeal_df <- taxa_abund_df %>% filter(Domain == "Archaea")

# 3. Process bacterial data
bacterial_agg <- bacterial_df %>% 
  group_by(Order) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(total))

bacterial_top10 <- bacterial_agg$Order[1:min(10, nrow(bacterial_agg))]

bacterial_plot_data <- bacterial_df %>% 
  mutate(Major_Taxa = ifelse(Order %in% bacterial_top10, Order, "Other")) %>% 
  group_by(Major_Taxa) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  # Normalize to 100%
  mutate(across(where(is.numeric), ~ .x / sum(.x) * 100)) %>% 
  melt(id.vars = "Major_Taxa", variable.name = "Sample", value.name = "Proportion") %>% 
  merge(sample_info_tab, by.x = "Sample", by.y = "row.names")

# 4. Process archaeal data
archaeal_agg <- archaeal_df %>% 
  group_by(Order) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(total))

archaeal_top10 <- archaeal_agg$Order[1:min(10, nrow(archaeal_agg))]

archaeal_plot_data <- archaeal_df %>% 
  mutate(Major_Taxa = ifelse(Order %in% archaeal_top10, Order, "Other")) %>% 
  group_by(Major_Taxa) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  # Normalize to 100%
  mutate(across(where(is.numeric), ~ .x / sum(.x) * 100)) %>% 
  melt(id.vars = "Major_Taxa", variable.name = "Sample", value.name = "Proportion") %>% 
  merge(sample_info_tab, by.x = "Sample", by.y = "row.names")

# Create color palettes
bacterial_colors <- c(
  "#8FB9AA", "#F2D096", "#F2A679", "#D95D39", "#35524A",
  "#7C9EB2", "#C9ADA7", "#9A8C98", "#797596", "#4A4E69",
  "#A9927D"  # For "Other" category
)

archaeal_colors <- c(
  "#E6B89C", "#FECEA8", "#FF847C", "#E84A5F", "#2A363B",
  "#F7C59F", "#D6A2AD", "#B07BAC", "#774C60", "#CFBACE",
  "#7D7ABC"  # For "Other" category
)

# 6. Generate the plots
bacterial_plot <- ggplot(bacterial_plot_data, aes(x = Sample, y = Proportion, fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = bacterial_colors) +
  labs(title = "Top 10 Bacterial Orderes",
       y = "Relative Abundance (%)",
       x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

archaeal_plot <- ggplot(archaeal_plot_data, aes(x = Sample, y = Proportion, fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = archaeal_colors) +
  labs(title = "Top 10 Archaeal Orderes",
       y = "Relative Abundance (%)",
       x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

# 7. Arrange and save plots
combined_plots <- bacterial_plot / archaeal_plot + 
  plot_layout(heights = c(1, 1))

#ggsave("relative_abundance_plots.png", combined_plots, width = 12, height = 10, dpi = 300)

# 8. Show plots
combined_plots

# ------- ORDER BY LAKES ----
# 1. Update the bacterial plot with boxes around lakes
bacterial_plot <- ggplot(bacterial_plot_data, 
                         aes(x = factor(Sample, levels = lake_order), 
                             y = Proportion, 
                             fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ lago, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = bacterial_colors) +
  labs(title = "Top 10 Bacterial Orderes by Lake",
       y = "Relative Abundance (%)",
       x = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.2, "lines"),
    panel.background = element_rect(color = "black", fill = NA),  # Add panel borders
    strip.background = element_rect(fill = "gray90", color = "black")  # Style facet labels
  )

# 2. Update the archaeal plot similarly
archaeal_plot <- ggplot(archaeal_plot_data, 
                        aes(x = factor(Sample, levels = lake_order), 
                            y = Proportion, 
                            fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ lago, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = archaeal_colors) +
  labs(title = "Top 10 Archaeal Orderes by Lake",
       y = "Relative Abundance (%)",
       x = "Sample") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.2, "lines"),
    panel.background = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "gray90", color = "black")
  )

# 3. Combine plots
combined_plots <- bacterial_plot / archaeal_plot + 
  plot_layout(heights = c(1, 1))

# 4. Save with wider aspect ratio
ggsave(file="Fig4_StackBars_By_Order.svg", plot=combined_plots, width=8, height=6)

# Show plots
combined_plots
#---- FAMILIA ------
# 1. Prepare the data
# Create combined taxonomy-abundance dataframe
taxa_abund_df <- data.frame(
  ASV = rownames(tax_tab_processed),
  Domain = tax_tab_processed[, "domain"],
  Family = tax_tab_processed[, "family"],
  stringsAsFactors = FALSE
) %>% 
  mutate(Family = ifelse(is.na(Family), "Unfamilyified", Family),
         Domain = ifelse(is.na(Domain), "Unassigned", Domain))

# Add relative abundances
for(sample in colnames(rel_abund)) {
  taxa_abund_df[[sample]] <- rel_abund[taxa_abund_df$ASV, sample]
}

# 2. Split into bacterial and archaeal datasets
bacterial_df <- taxa_abund_df %>% filter(Domain == "Bacteria")
archaeal_df <- taxa_abund_df %>% filter(Domain == "Archaea")

# 3. Process bacterial data
bacterial_agg <- bacterial_df %>% 
  group_by(Family) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(total))

bacterial_top10 <- bacterial_agg$Family[1:min(10, nrow(bacterial_agg))]

bacterial_plot_data <- bacterial_df %>% 
  mutate(Major_Taxa = ifelse(Family %in% bacterial_top10, Family, "Other")) %>% 
  group_by(Major_Taxa) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  # Normalize to 100%
  mutate(across(where(is.numeric), ~ .x / sum(.x) * 100)) %>% 
  melt(id.vars = "Major_Taxa", variable.name = "Sample", value.name = "Proportion") %>% 
  merge(sample_info_tab, by.x = "Sample", by.y = "row.names")

# 4. Process archaeal data
archaeal_agg <- archaeal_df %>% 
  group_by(Family) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(total))

archaeal_top10 <- archaeal_agg$Family[1:min(10, nrow(archaeal_agg))]

archaeal_plot_data <- archaeal_df %>% 
  mutate(Major_Taxa = ifelse(Family %in% archaeal_top10, Family, "Other")) %>% 
  group_by(Major_Taxa) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  # Normalize to 100%
  mutate(across(where(is.numeric), ~ .x / sum(.x) * 100)) %>% 
  melt(id.vars = "Major_Taxa", variable.name = "Sample", value.name = "Proportion") %>% 
  merge(sample_info_tab, by.x = "Sample", by.y = "row.names")

# Create color palettes
bacterial_colors <- c(
  "#8FB9AA", "#F2D096", "#F2A679", "#D95D39", "#35524A",
  "#7C9EB2", "#C9ADA7", "#9A8C98", "#797596", "#4A4E69",
  "#A9927D"  # For "Other" category
)

archaeal_colors <- c(
  "#E6B89C", "#FECEA8", "#FF847C", "#E84A5F", "#2A363B",
  "#F7C59F", "#D6A2AD", "#B07BAC", "#774C60", "#CFBACE",
  "#7D7ABC"  # For "Other" category
)

# 6. Generate the plots
bacterial_plot <- ggplot(bacterial_plot_data, aes(x = Sample, y = Proportion, fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = bacterial_colors) +
  labs(title = "Top 10 Bacterial Familyes",
       y = "Relative Abundance (%)",
       x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

archaeal_plot <- ggplot(archaeal_plot_data, aes(x = Sample, y = Proportion, fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = archaeal_colors) +
  labs(title = "Top 10 Archaeal Familyes",
       y = "Relative Abundance (%)",
       x = "Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right")

# 7. Arrange and save plots
combined_plots <- bacterial_plot / archaeal_plot + 
  plot_layout(heights = c(1, 1))

#ggsave("relative_abundance_plots.png", combined_plots, width = 12, height = 10, dpi = 300)

# 8. Show plots
combined_plots

# ------- FAMILIA BY LAKES ----
# 1. Update the bacterial plot with boxes around lakes
bacterial_plot <- ggplot(bacterial_plot_data, 
                         aes(x = factor(Sample, levels = lake_order), 
                             y = Proportion, 
                             fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ lago, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = bacterial_colors) +
  labs(title = "Top 10 Bacterial Familyes by Lake",
       y = "Relative Abundance (%)",
       x = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.2, "lines"),
    panel.background = element_rect(color = "black", fill = NA),  # Add panel bfamilys
    strip.background = element_rect(fill = "gray90", color = "black")  # Style facet labels
  )

# 2. Update the archaeal plot similarly
archaeal_plot <- ggplot(archaeal_plot_data, 
                        aes(x = factor(Sample, levels = lake_order), 
                            y = Proportion, 
                            fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ lago, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = archaeal_colors) +
  labs(title = "Top 10 Archaeal Familyes by Lake",
       y = "Relative Abundance (%)",
       x = "Sample") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.2, "lines"),
    panel.background = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "gray90", color = "black")
  )

# 3. Combine plots
combined_plots <- bacterial_plot / archaeal_plot + 
  plot_layout(heights = c(1, 1))

# 4. Save with wider aspect ratio
ggsave(file="Fig4_StackBars_By_Family.svg", plot=combined_plots, width=8, height=6)

# Show plots
combined_plots

#---- GENUS ------

# 1. Prepare the data
taxa_abund_df <- data.frame(
  ASV = rownames(tax_tab_processed),
  Domain = tax_tab_processed[, "domain"],
  Genus = tax_tab_processed[, "genus"],
  stringsAsFactors = FALSE
) %>% 
  mutate(Genus = ifelse(is.na(Genus), "Unclassified", Genus),
         Domain = ifelse(is.na(Domain), "Unassigned", Domain))

# Add relative abundances
for(sample in colnames(rel_abund)) {
  taxa_abund_df[[sample]] <- rel_abund[taxa_abund_df$ASV, sample]
}

# 2. Split into bacterial and archaeal datasets
bacterial_df <- taxa_abund_df %>% filter(Domain == "Bacteria")
archaeal_df <- taxa_abund_df %>% filter(Domain == "Archaea")

# 3. Process bacterial data
bacterial_agg <- bacterial_df %>% 
  group_by(Genus) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(total))

bacterial_top10 <- bacterial_agg$Genus[1:min(10, nrow(bacterial_agg))]

bacterial_plot_data <- bacterial_df %>% 
  mutate(Major_Taxa = ifelse(Genus %in% bacterial_top10, Genus, "Other")) %>% 
  group_by(Major_Taxa) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(across(where(is.numeric), ~ .x / sum(.x) * 100)) %>% 
  melt(id.vars = "Major_Taxa", variable.name = "Sample", value.name = "Proportion") %>% 
  merge(sample_info_tab, by.x = "Sample", by.y = "row.names")

# 4. Process archaeal data
archaeal_agg <- archaeal_df %>% 
  group_by(Genus) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  arrange(desc(total))

archaeal_top10 <- archaeal_agg$Genus[1:min(10, nrow(archaeal_agg))]

archaeal_plot_data <- archaeal_df %>% 
  mutate(Major_Taxa = ifelse(Genus %in% archaeal_top10, Genus, "Other")) %>% 
  group_by(Major_Taxa) %>% 
  summarise(across(where(is.numeric), sum)) %>% 
  mutate(across(where(is.numeric), ~ .x / sum(.x) * 100)) %>% 
  melt(id.vars = "Major_Taxa", variable.name = "Sample", value.name = "Proportion") %>% 
  merge(sample_info_tab, by.x = "Sample", by.y = "row.names")

# 5. Color palettes
bacterial_colors <- c(
  "#8FB9AA", "#F2D096", "#F2A679", "#D95D39", "#35524A",
  "#7C9EB2", "#C9ADA7", "#9A8C98", "#797596", "#4A4E69",
  "#A9927D"
)

archaeal_colors <- c(
  "#E6B89C", "#FECEA8", "#FF847C", "#E84A5F", "#2A363B",
  "#F7C59F", "#D6A2AD", "#B07BAC", "#774C60", "#CFBACE",
  "#7D7ABC"
)

# 6. Plots by lake
bacterial_plot <- ggplot(bacterial_plot_data, 
                         aes(x = factor(Sample, levels = lake_order), 
                             y = Proportion, 
                             fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ lago, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = bacterial_colors) +
  labs(y = "Relative Abundance (%)", x = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.2, "lines"),
    panel.background = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "gray90", color = "black")
  )

archaeal_plot <- ggplot(archaeal_plot_data, 
                        aes(x = factor(Sample, levels = lake_order), 
                            y = Proportion, 
                            fill = Major_Taxa)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ lago, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = archaeal_colors) +
  labs(y = "Relative Abundance (%)", x = "Sample") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.2, "lines"),
    panel.background = element_rect(color = "black", fill = NA),
    strip.background = element_rect(fill = "gray90", color = "black")
  )

combined_plots <- bacterial_plot / archaeal_plot + 
  plot_layout(heights = c(1, 1))

ggsave(file="Fig_StackBars_By_Genus.svg", plot=combined_plots, width=8, height=6)

combined_plots
