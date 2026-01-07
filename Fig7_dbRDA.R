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

tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam-v2.tsv", header=T, row.names=1, check.names=F, sep="\t"))
print("Successfully read taxonomy file")
sample_info_tab <- read.csv("pantano_metadata.csv", header=T, row.names=1, check.names=F)
print("Successfully read metadata file")
count_tab <- read.table("ASVs_counts-no-contam-v2.tsv", header=T, row.names=1, check.names=F, sep="\t")[ , -c(7:9)]
print("Successfully read counts file")

colnames(count_tab) <- gsub("reads/", "", colnames(count_tab))

# Crear objeto phyloseq (opcional, pero útil para integración)
ps <- phyloseq(
  otu_table(count_tab, taxa_are_rows = TRUE),
  tax_table(tax_tab),
  sample_data(sample_info_tab)
)

# Normalización Hellinger (para RDA y PCoA)
#otu_hell <- decostand(count_tab, method = "hellinger")

# Matriz de distancias Bray-Curtis
#dist_bray <- vegdist(otu_hell, method = "bray")

# Seleccionar variables ambientales (excluyendo NA)
#env_vars <- sample_info_tab[, c("temp", "conduct", "pH", "oxy", "solid_disuelto")]
#env_vars <- na.omit(env_vars)

#-----
ps.ord <- ordinate(ps, method = "NMDS", distance = "bray")
p_lagos_1 <- plot_ordination(
  physeq = ps,
  ordination = ps.ord,
  type = "samples",
  color = "lago",  # Color por grupos de lago
  title = "NMDS coloreado por lago (colores personalizados)"
) +
  geom_point(size = 3, aes(color = lago)) +  # Mapear color a la variable "lago"
  stat_ellipse(level = 0.95) +
  scale_color_manual(  # Asignar colores manualmente usando "color_lago"
    values = unique(sample_data(ps)$color_lago)[order(unique(sample_data(ps)$lago))],
    labels = unique(sample_data(ps)$lago)
  ) +
  theme_minimal()

p_lagos_1

p_lagos <- plot_ordination(
  physeq = ps,
  ordination = ps.ord,
  type = "samples",
  color = "lago"
) +
  geom_point(
    size = 3,
    aes(fill = lago),  # Relleno por lago
    shape = 21,        # Forma con borde
    color = "black",   # Borde negro
    stroke = 0.5       # Grosor del borde
  ) +
  geom_text(
    aes(label = rownames(sample_data(ps))),  # Etiquetas de muestras
    color = "black",
    size = 3,
    hjust = -0.2,
    vjust = 0.5
  ) +
  stat_ellipse(
    aes(color = lago),  # Elipses con colores por lago
    level = 0.95,
    linewidth = 0.5
  ) +
  scale_fill_manual(
    values = unique(sample_data(ps)$color_lago)[order(unique(sample_data(ps)$lago))]
  ) +
  scale_color_manual(
    values = unique(sample_data(ps)$color_lago)[order(unique(sample_data(ps)$lago))]
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Marco negro
    legend.position = "right"
  )

p_lagos

#ggsave(file="Fig6_NMDS_plot.svg", plot=parameters_plot, width=8, height=6)


#-----

# First, let's separate our phyloseq object into bacteria and archaea
# Note: Modify the Domain column name if your taxonomy has a different name

# Check taxonomy structure to find the Domain column
tax_colnames <- colnames(tax_tab)
print(tax_colnames)

# Assuming the first column is Domain (adjust if needed)
domain_col <- 1  # Change this if your Domain is in a different column

# Create bacteria-only phyloseq object
bacteria_indices <- which(tax_tab[,domain_col] == "Bacteria")
bacteria_ps <- prune_taxa(taxa_names(ps)[bacteria_indices], ps)

# Create archaea-only phyloseq object
archaea_indices <- which(tax_tab[,domain_col] == "Archaea")
archaea_ps <- prune_taxa(taxa_names(ps)[archaea_indices], ps)

# Print the number of taxa in each domain
cat("Total ASVs:", ntaxa(ps), "\n")
cat("Bacterial ASVs:", ntaxa(bacteria_ps), "\n")
cat("Archaeal ASVs:", ntaxa(archaea_ps), "\n")

# Extract environmental variables (excluding colors and categorical)
env_vars <- sample_data(ps)[, c("temp", "conduct", "pH", "oxy", "salinidad")]
env_vars <- as(env_vars, "data.frame")
# Handle NaN values
env_vars <- na.omit(env_vars)

# For each domain, we'll perform:
# 1. Hellinger transformation
# 2. db-RDA with environmental variables
# 3. Significance testing
# 4. Visualization

# Function to perform the complete analysis for a given phyloseq object
analyze_domain <- function(ps_obj, domain_name, env_data) {
  # Subset to samples that have environmental data
  ps_obj <- prune_samples(rownames(env_data), ps_obj)
  
  # Skip if we don't have enough data
  if(ntaxa(ps_obj) < 2 || nsamples(ps_obj) < 3) {
    cat("Not enough data for", domain_name, "analysis\n")
    return(NULL)
  }
  
  # Hellinger transformation
  otu_mat <- as(otu_table(ps_obj), "matrix")
  hell_transformed <- decostand(otu_mat, method="hellinger")
  
  # Calculate Bray-Curtis distances
  bray_dist <- vegdist(t(hell_transformed), method="bray")
  
  # Perform db-RDA
  dbrda_result <- try(dbrda(bray_dist ~ ., data=env_data), silent=TRUE)
  
  if(inherits(dbrda_result, "try-error")) {
    cat("Error in db-RDA for", domain_name, "\n")
    return(NULL)
  }
  
  # Calculate variance explained
  var_explained <- round(100 * eigenvals(dbrda_result)/sum(eigenvals(dbrda_result)), 1)
  
  # Test significance 
  global_test <- anova(dbrda_result)
  var_test <- anova(dbrda_result, by="terms")
  
  # Perform NMDS for comparison
  nmds_result <- metaMDS(bray_dist, trymax=100)
  
  # Return results as a list
  return(list(
    dbrda = dbrda_result,
    global_test = global_test,
    var_test = var_test,
    var_explained = var_explained,
    nmds = nmds_result,
    bray_dist = bray_dist,
    hell_transformed = hell_transformed
  ))
}

# 1. Crear un vector de colores con nombres según los lagos
color_mapping <- sample_data(ps)$color_lago
names(color_mapping) <- sample_data(ps)$lago
color_mapping <- color_mapping[!duplicated(names(color_mapping))]

# Analyze each domain
bacteria_results <- analyze_domain(bacteria_ps, "Bacteria", env_vars)
archaea_results <- analyze_domain(archaea_ps, "Archaea", env_vars)

# Plot the results for bacteria
if(!is.null(bacteria_results)) {
  # Plot db-RDA
  dbrda_bact <- bacteria_results$dbrda
  
  # Extract scores and check column names
  site_scores <- scores(dbrda_bact, display="sites")
  print("Column names of site scores:")
  print(colnames(site_scores))
  
  # Create plot data - use actual column names
  site_data <- data.frame(site_scores, 
                          sample_data(bacteria_ps)[rownames(site_scores),])
  
  # Extract arrow data
  arrow_data <- scores(dbrda_bact, display="bp")
  print("Column names of biplot scores:")
  print(colnames(arrow_data))
  arrow_df <- data.frame(arrow_data)
  arrow_df$label <- rownames(arrow_df)
  
  # Create the plot using the correct column names
  # Replace CAP1 and CAP2 with the actual column names from your output
  axis1 <- colnames(site_scores)[1]  # First axis name
  axis2 <- colnames(site_scores)[2]  # Second axis name
  
  bact_plot <- ggplot() +
    geom_point(data=site_data, aes_string(x=axis1, y=axis2, color="lago"), size=3) +
    geom_text(data=site_data, aes_string(x=axis1, y=axis2, label="rownames(site_data)"), 
              hjust=-0.2, vjust=0, size=3) +
    geom_segment(data=arrow_df, aes_string(x=0, y=0, xend=paste(axis1, "*2", sep=""), 
                                           yend=paste(axis2, "*2", sep="")),
                 arrow=arrow(length=unit(0.2, "cm")), color="red") +
    geom_text(data=arrow_df, aes_string(x=paste(axis1, "*2.2", sep=""), 
                                        y=paste(axis2, "*2.2", sep=""), label="label"),
              color="red", size=4) +
    theme_bw() +
    labs(title="db-RDA for Bacterial Communities",
         subtitle=paste(axis1, ":", bacteria_results$var_explained[1], "%, ", 
                        axis2, ":", bacteria_results$var_explained[2], "%"),
         caption=paste("Global test p-value:", 
                       round(bacteria_results$global_test$`Pr(>F)`[1], 3))) +
    scale_color_manual(
      values = color_mapping,  # Usa el mapeo de colores exacto
      breaks = names(color_mapping)  # Asegura el orden correcto
    ) +
    coord_equal()
  
  bact_plot
  
  # Print significance test results
  cat("\n========== Bacteria db-RDA Results ==========\n")
  print(bacteria_results$global_test)
  cat("\n---- Effects of individual variables ----\n")
  print(bacteria_results$var_test)
}

# Plot the results for archaea
if(!is.null(archaea_results)) {
  # Plot db-RDA
  dbrda_arch <- archaea_results$dbrda
  
  # Extract scores and check column names
  site_scores <- scores(dbrda_arch, display="sites")
  print("Column names of archaea site scores:")
  print(colnames(site_scores))
  
  # Create plot data - use actual column names
  site_data <- data.frame(site_scores, 
                          sample_data(archaea_ps)[rownames(site_scores),])
  
  # Extract arrow data
  arrow_data <- scores(dbrda_arch, display="bp")
  print("Column names of archaea biplot scores:")
  print(colnames(arrow_data))
  arrow_df <- data.frame(arrow_data)
  arrow_df$label <- rownames(arrow_df)
  
  # Create the plot using the correct column names
  # Use the actual column names from your output
  axis1 <- colnames(site_scores)[1]  # First axis name
  axis2 <- colnames(site_scores)[2]  # Second axis name
  
  arch_plot <- ggplot() +
    geom_point(data=site_data, aes_string(x=axis1, y=axis2, color="lago"), size=3) +
    geom_text(data=site_data, aes_string(x=axis1, y=axis2, label="rownames(site_data)"), 
              hjust=-0.2, vjust=0, size=3) +
    geom_segment(data=arrow_df, aes_string(x=0, y=0, xend=paste(axis1, "*2", sep=""), 
                                           yend=paste(axis2, "*2", sep="")),
                 arrow=arrow(length=unit(0.2, "cm")), color="red") +
    geom_text(data=arrow_df, aes_string(x=paste(axis1, "*2.2", sep=""), 
                                        y=paste(axis2, "*2.2", sep=""), label="label"),
              color="red", size=4) +
    theme_bw() +
    labs(title="db-RDA for Archaeal Communities",
         subtitle=paste(axis1, ":", archaea_results$var_explained[1], "%, ", 
                        axis2, ":", archaea_results$var_explained[2], "%"),
         caption=paste("Global test p-value:", 
                       round(archaea_results$global_test$`Pr(>F)`[1], 3))) +
    scale_color_manual(
      values = color_mapping,  # Usa el mapeo de colores exacto
      breaks = names(color_mapping)  # Asegura el orden correcto
    ) + theme(legend.position = "none") +
    coord_equal()
  
  arch_plot
  
  # Print significance test results
  cat("\n========== Archaea db-RDA Results ==========\n")
  print(archaea_results$global_test)
  cat("\n---- Effects of individual variables ----\n")
  print(archaea_results$var_test)
}


#------
# Update the env_vars selection to include solid_disuelto
env_vars <- sample_data(ps)[, c("temp", "conduct", "pH", "oxy", "salinidad", "solid_disuelto")]
env_vars <- as(env_vars, "data.frame")

# Handle NaN/NA values - replace NaN with NA first
env_vars[] <- lapply(env_vars, function(x) ifelse(is.nan(x), NA, x))
env_vars <- na.omit(env_vars)  # Now properly removes rows with NA

# The rest of your analysis code remains the same, but now includes solid_disuelto:
analyze_domain <- function(ps_obj, domain_name, env_data) {
  # Subset to samples that have environmental data
  ps_obj <- prune_samples(rownames(env_data), ps_obj)
  
  if(ntaxa(ps_obj) < 2 || nsamples(ps_obj) < 3) {
    cat("Not enough data for", domain_name, "analysis\n")
    return(NULL)
  }
  
  # Hellinger transformation
  otu_mat <- as(otu_table(ps_obj), "matrix")
  hell_transformed <- decostand(otu_mat, method="hellinger")
  
  # Calculate Bray-Curtis distances
  bray_dist <- vegdist(t(hell_transformed), method="bray")
  
  # Perform db-RDA with all variables including solid_disuelto
  dbrda_result <- try(dbrda(bray_dist ~ ., data=env_data), silent=TRUE)
  
  if(inherits(dbrda_result, "try-error")) {
    cat("Error in db-RDA for", domain_name, "\n")
    return(NULL)
  }
  
  # Calculate variance explained
  var_explained <- round(100 * eigenvals(dbrda_result)/sum(eigenvals(dbrda_result)), 1)
  
  # Test significance 
  global_test <- anova(dbrda_result)
  var_test <- anova(dbrda_result, by="terms")
  
  # Perform NMDS for comparison
  nmds_result <- metaMDS(bray_dist, trymax=100)
  
  # Return results as a list
  return(list(
    dbrda = dbrda_result,
    global_test = global_test,
    var_test = var_test,
    var_explained = var_explained,
    nmds = nmds_result,
    bray_dist = bray_dist,
    hell_transformed = hell_transformed
  ))
}

# 1. Crear un vector de colores con nombres según los lagos
color_mapping <- sample_data(ps)$color_lago
names(color_mapping) <- sample_data(ps)$lago
color_mapping <- color_mapping[!duplicated(names(color_mapping))]

# Analyze each domain
bacteria_results <- analyze_domain(bacteria_ps, "Bacteria", env_vars)
archaea_results <- analyze_domain(archaea_ps, "Archaea", env_vars)

# Plot the results for bacteria
if(!is.null(bacteria_results)) {
  # Plot db-RDA
  dbrda_bact <- bacteria_results$dbrda
  
  # Extract scores and check column names
  site_scores <- scores(dbrda_bact, display="sites")
  print("Column names of site scores:")
  print(colnames(site_scores))
  
  # Create plot data - use actual column names
  site_data <- data.frame(site_scores, 
                          sample_data(bacteria_ps)[rownames(site_scores),])
  
  # Extract arrow data
  arrow_data <- scores(dbrda_bact, display="bp")
  print("Column names of biplot scores:")
  print(colnames(arrow_data))
  arrow_df <- data.frame(arrow_data)
  arrow_df$label <- rownames(arrow_df)
  
  # Create the plot using the correct column names
  # Replace CAP1 and CAP2 with the actual column names from your output
  axis1 <- colnames(site_scores)[1]  # First axis name
  axis2 <- colnames(site_scores)[2]  # Second axis name
  
  bact_plot <- ggplot() +
    geom_point(data=site_data, aes_string(x=axis1, y=axis2, color="lago"), size=3) +
    geom_text(data=site_data, aes_string(x=axis1, y=axis2, label="rownames(site_data)"), 
              hjust=-0.2, vjust=0, size=3) +
    geom_segment(data=arrow_df, aes_string(x=0, y=0, xend=paste(axis1, "*2", sep=""), 
                                           yend=paste(axis2, "*2", sep="")),
                 arrow=arrow(length=unit(0.2, "cm")), color="red") +
    geom_text(data=arrow_df, aes_string(x=paste(axis1, "*2.2", sep=""), 
                                        y=paste(axis2, "*2.2", sep=""), label="label"),
              color="red", size=4) +
    theme_bw() +
    labs(title="db-RDA for Bacterial Communities",
         subtitle=paste(axis1, ":", bacteria_results$var_explained[1], "%, ", 
                        axis2, ":", bacteria_results$var_explained[2], "%"),
         caption=paste("Global test p-value:", 
                       round(bacteria_results$global_test$`Pr(>F)`[1], 3))) +
    scale_color_manual(
      values = color_mapping,  # Usa el mapeo de colores exacto
      breaks = names(color_mapping)  # Asegura el orden correcto
    ) +
    coord_equal()
  
  bact_plot
  
  # Print significance test results
  cat("\n========== Bacteria db-RDA Results ==========\n")
  print(bacteria_results$global_test)
  cat("\n---- Effects of individual variables ----\n")
  print(bacteria_results$var_test)
}

# Plot the results for archaea
if(!is.null(archaea_results)) {
  # Plot db-RDA
  dbrda_arch <- archaea_results$dbrda
  
  # Extract scores and check column names
  site_scores <- scores(dbrda_arch, display="sites")
  print("Column names of archaea site scores:")
  print(colnames(site_scores))
  
  # Create plot data - use actual column names
  site_data <- data.frame(site_scores, 
                          sample_data(archaea_ps)[rownames(site_scores),])
  
  # Extract arrow data
  arrow_data <- scores(dbrda_arch, display="bp")
  print("Column names of archaea biplot scores:")
  print(colnames(arrow_data))
  arrow_df <- data.frame(arrow_data)
  arrow_df$label <- rownames(arrow_df)
  
  # Create the plot using the correct column names
  # Use the actual column names from your output
  axis1 <- colnames(site_scores)[1]  # First axis name
  axis2 <- colnames(site_scores)[2]  # Second axis name
  
  arch_plot <- ggplot() +
    geom_point(data=site_data, aes_string(x=axis1, y=axis2, color="lago"), size=3) +
    geom_text(data=site_data, aes_string(x=axis1, y=axis2, label="rownames(site_data)"), 
              hjust=-0.2, vjust=0, size=3) +
    geom_segment(data=arrow_df, aes_string(x=0, y=0, xend=paste(axis1, "*2", sep=""), 
                                           yend=paste(axis2, "*2", sep="")),
                 arrow=arrow(length=unit(0.2, "cm")), color="red") +
    geom_text(data=arrow_df, aes_string(x=paste(axis1, "*2.2", sep=""), 
                                        y=paste(axis2, "*2.2", sep=""), label="label"),
              color="red", size=4) +
    theme_bw() +
    labs(title="db-RDA for Archaeal Communities",
         subtitle=paste(axis1, ":", archaea_results$var_explained[1], "%, ", 
                        axis2, ":", archaea_results$var_explained[2], "%"),
         caption=paste("Global test p-value:", 
                       round(archaea_results$global_test$`Pr(>F)`[1], 3))) +
    scale_color_manual(
      values = color_mapping,  # Usa el mapeo de colores exacto
      breaks = names(color_mapping)  # Asegura el orden correcto
    ) + theme(legend.position = "none") +
    coord_equal()
  
  arch_plot
  
  # Print significance test results
  cat("\n========== Archaea db-RDA Results ==========\n")
  print(archaea_results$global_test)
  cat("\n---- Effects of individual variables ----\n")
  print(archaea_results$var_test)
  
  # 3. Combinar manteniendo 1 sola leyenda
  combined_plots <- arch_plot + bact_plot +
    plot_layout(guides = "collect") +  # Unifica leyendas
    plot_annotation(tag_levels = "A")  # Añade etiquetas (A), (B)
  
  # 4. Mostrar resultado
  print(combined_plots)
}
