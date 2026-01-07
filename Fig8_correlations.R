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

# Identificar ASVs con conteos cero
zero_count_asvs <- rowSums(count_tab) == 0

# Remove ASVs with zero counts and remove ASVs in tax_tab
filtered_count_tab <- count_tab[!zero_count_asvs, ]
filtered_tax_tab <- tax_tab[rownames(filtered_count_tab), ]

count_tab <- filtered_count_tab
tax_tab <- filtered_tax_tab

# Crear objeto phyloseq con tus datos
ps <- phyloseq(
  otu_table(as.matrix(count_tab), taxa_are_rows = TRUE),
  sample_data(sample_info_tab),
  tax_table(tax_tab)
)

# Agregar los conteos a nivel de género
ps_genus <- tax_glom(ps, taxrank = "genus", NArm = FALSE)

# Extraer la tabla de abundancia a nivel de género
genus_counts <- as.data.frame(t(otu_table(ps_genus)))

# Obtener los nombres de los géneros para mejorar la lectura
tax_genus <- as.data.frame(tax_table(ps_genus))
genus_names <- tax_genus$genus
genus_names[is.na(genus_names)] <- paste0("Unknown_", seq_along(genus_names[is.na(genus_names)]))
colnames(genus_counts) <- genus_names

# Transformación de los datos (CLR - recomendada para datos composicionales)
library(microbiome)
genus_counts_transformed <- microbiome::transform(genus_counts, "clr")

# Obtener los datos ambientales
env_data <- as.data.frame(sample_data(ps))
# Seleccionar solo las variables numéricas relevantes
env_vars <- env_data[, c("temp", "conduct", "pH", "oxy", "solid_disuelto", "salinidad")]
# Manejar NAs en variables ambientales
env_vars <- na.omit(env_vars)

# Asegurarse de que solo usamos las muestras que tienen datos completos
common_samples <- intersect(rownames(genus_counts_transformed), rownames(env_vars))
genus_filtered <- genus_counts_transformed[common_samples, ]
env_filtered <- env_vars[common_samples, ]

# Calcular correlaciones entre géneros y variables ambientales
library(Hmisc)
correlation_matrix <- rcorr(as.matrix(cbind(genus_filtered, env_filtered)), type = "pearson")

# Extraer las correlaciones específicas entre géneros y variables ambientales
n_genera <- ncol(genus_filtered)
n_env <- ncol(env_filtered)
cor_genus_env <- correlation_matrix$r[(n_genera+1):(n_genera+n_env), 1:n_genera]
p_genus_env <- correlation_matrix$P[(n_genera+1):(n_genera+n_env), 1:n_genera]

# Transponer para facilitar la visualización
cor_genus_env <- t(cor_genus_env)
p_genus_env <- t(p_genus_env)

# ============================================
# CORRECCIÓN POR COMPARACIONES MÚLTIPLES (FDR)
# ============================================

# Convertir matriz de p-valores a vector
p_values_vector <- as.vector(p_genus_env)

# Aplicar corrección Benjamini-Hochberg
p_adjusted <- p.adjust(p_values_vector, method = "BH")

# Reconstruir la matriz con p-valores ajustados
p_genus_env_adjusted <- matrix(p_adjusted, 
                                nrow = nrow(p_genus_env), 
                                ncol = ncol(p_genus_env),
                                dimnames = dimnames(p_genus_env))

# Comparar cuántas correlaciones se pierden
cat("====== COMPARACIÓN DE CORRECCIÓN FDR ======\n")
cat("Correlaciones significativas SIN corrección (p < 0.05):", sum(p_genus_env < 0.05, na.rm = TRUE), "\n")
cat("Correlaciones significativas CON corrección FDR (p.adj < 0.05):", sum(p_genus_env_adjusted < 0.05, na.rm = TRUE), "\n")
cat("Correlaciones perdidas:", sum(p_genus_env < 0.05, na.rm = TRUE) - sum(p_genus_env_adjusted < 0.05, na.rm = TRUE), "\n")
cat("============================================\n")

# Filtrar correlaciones significativas CON CORRECCIÓN FDR
sig_correlations <- cor_genus_env
sig_correlations[p_genus_env_adjusted >= 0.05] <- NA

# Preparar datos para visualización
cor_long <- as.data.frame(as.table(sig_correlations))
names(cor_long) <- c("Genus", "Environmental_Variable", "Correlation")
cor_long <- na.omit(cor_long)

# Verificar si quedaron correlaciones significativas
if(nrow(cor_long) == 0) {
  cat("\n¡ATENCIÓN! No quedaron correlaciones significativas después de la corrección FDR.\n")
  cat("Considera usar un umbral menos estricto (p.adj < 0.1) o reportar que no hay asociaciones significativas.\n")
  
  # Opción alternativa: usar p.adj < 0.1
  sig_correlations_relaxed <- cor_genus_env
  sig_correlations_relaxed[p_genus_env_adjusted >= 0.1] <- NA
  cor_long_relaxed <- as.data.frame(as.table(sig_correlations_relaxed))
  names(cor_long_relaxed) <- c("Genus", "Environmental_Variable", "Correlation")
  cor_long_relaxed <- na.omit(cor_long_relaxed)
  cat("Correlaciones significativas con p.adj < 0.1:", nrow(cor_long_relaxed), "\n")
}

# Ordenar por valor absoluto de correlación
cor_long <- cor_long[order(abs(cor_long$Correlation), decreasing = TRUE), ]

# Filtrar para quitar los géneros "Unknown_"
cor_long_no_unknown <- cor_long[!grepl("^Unknown_", cor_long$Genus), ]

# Ordenar por valor absoluto de correlación
cor_long_no_unknown <- cor_long_no_unknown[order(abs(cor_long_no_unknown$Correlation), decreasing = TRUE), ]

# Seleccionar los N géneros con las correlaciones más fuertes (sin Unknown)
top_n <- 30
if(nrow(cor_long_no_unknown) > top_n) {
  top_genera <- unique(cor_long_no_unknown$Genus)[1:min(top_n, length(unique(cor_long_no_unknown$Genus)))]
  cor_long_filtered <- cor_long_no_unknown[cor_long_no_unknown$Genus %in% top_genera, ]
} else {
  cor_long_filtered <- cor_long_no_unknown
}

# Verificar que hay datos para graficar
if(nrow(cor_long_filtered) > 0) {
  
  # Redondear los valores para mejor visualización
  cor_long_filtered$Correlation_rounded <- round(cor_long_filtered$Correlation, 2)
  
  # Crear heatmap con valores numéricos
  genus_corr <- ggplot(cor_long_filtered, aes(x = Environmental_Variable, y = Genus, fill = Correlation)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = Correlation_rounded), 
              color = ifelse(abs(cor_long_filtered$Correlation) > 0.5, "white", "black"),
              size = 3) +
    scale_fill_gradient2(
      low = "#2166AC", 
      mid = "white", 
      high = "#B2182B", 
      midpoint = 0,
      name = "Pearson\nCorrelation"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(face = "bold", size = 12),
      panel.grid = element_blank(),
      legend.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    ) +
    labs(x = "Environmental Variables", 
         y = "Genus",
         title = "Significant correlations between bacterial genera and environmental variables",
         subtitle = "FDR-corrected (Benjamini-Hochberg, p.adj < 0.05)")
  
  print(genus_corr)
  ggsave(file="Fig8_correlationsp_genus_FDR.svg", plot=genus_corr, width=8, height=6)
  
} else {
  cat("No hay suficientes datos para crear el heatmap después de la corrección FDR.\n")
}

# ============================================
# FUNCIÓN PARA HEATMAP POR DOMINIO (CON FDR)
# ============================================

generate_domain_correlation_heatmap <- function(ps, domain_name, n_top = 10) {
  # Filtrar por dominio
  tax_df <- as.data.frame(tax_table(ps))
  keep_taxa <- rownames(tax_df)[tax_df$domain == domain_name]
  ps_domain <- prune_taxa(keep_taxa, ps)
  
  # Verificar si hay datos para este dominio
  if(ntaxa(ps_domain) == 0) {
    message(paste0("No se encontraron taxa para el dominio: ", domain_name))
    return(NULL)
  }
  
  # Agregar los conteos a nivel de clase
  ps_class <- tax_glom(ps_domain, taxrank = "class", NArm = FALSE)
  
  # Extraer la tabla de abundancia a nivel de clase
  class_counts <- as.data.frame(t(otu_table(ps_class)))
  
  # Obtener los nombres de las clases
  tax_class <- as.data.frame(tax_table(ps_class))
  class_names <- tax_class$class
  class_names[is.na(class_names)] <- paste0("Unknown_", seq_along(class_names[is.na(class_names)]))
  colnames(class_counts) <- class_names
  
  # Calcular la abundancia total de cada clase
  class_abundance <- colSums(class_counts)
  
  # Verificar si hay suficientes clases
  n_avail <- min(n_top, length(class_abundance))
  if(n_avail == 0) {
    message(paste0("No se encontraron clases para el dominio: ", domain_name))
    return(NULL)
  }
  
  # Identificar las N clases más abundantes
  top_classes <- names(sort(class_abundance, decreasing = TRUE)[1:n_avail])
  class_counts_top <- class_counts[, top_classes]
  
  # Transformación de datos (CLR)
  if(requireNamespace("microbiome", quietly = TRUE)) {
    class_counts_transformed <- microbiome::transform(class_counts_top, "clr")
  } else {
    class_counts_transformed <- log1p(class_counts_top)
  }
  
  # Obtener los datos ambientales
  env_data <- as.data.frame(sample_data(ps))
  env_vars <- env_data[, c("temp", "conduct", "pH", "oxy", "solid_disuelto", "salinidad")]
  env_vars <- na.omit(env_vars)
  
  # Asegurarse de que solo usamos las muestras que tienen datos completos
  common_samples <- intersect(rownames(class_counts_transformed), rownames(env_vars))
  class_filtered <- class_counts_transformed[common_samples, ]
  env_filtered <- env_vars[common_samples, ]
  
  # Calcular correlaciones entre clases y variables ambientales
  correlation_matrix <- rcorr(as.matrix(cbind(class_filtered, env_filtered)), type = "pearson")
  
  # Extraer las correlaciones específicas entre clases y variables ambientales
  n_classes <- ncol(class_filtered)
  n_env <- ncol(env_filtered)
  cor_class_env <- correlation_matrix$r[(n_classes+1):(n_classes+n_env), 1:n_classes]
  p_class_env <- correlation_matrix$P[(n_classes+1):(n_classes+n_env), 1:n_classes]
  
  # Transponer para facilitar la visualización
  cor_class_env <- t(cor_class_env)
  p_class_env <- t(p_class_env)
  
  # ============================================
  # CORRECCIÓN FDR PARA CADA DOMINIO
  # ============================================
  p_values_vector <- as.vector(p_class_env)
  p_adjusted <- p.adjust(p_values_vector, method = "BH")
  p_class_env_adjusted <- matrix(p_adjusted, 
                                  nrow = nrow(p_class_env), 
                                  ncol = ncol(p_class_env),
                                  dimnames = dimnames(p_class_env))
  
  # Mostrar comparación
  cat(paste0("\n=== ", domain_name, " ===\n"))
  cat("Correlaciones significativas SIN corrección:", sum(p_class_env < 0.05, na.rm = TRUE), "\n")
  cat("Correlaciones significativas CON corrección FDR:", sum(p_class_env_adjusted < 0.05, na.rm = TRUE), "\n")
  
  # Preparar datos para visualización
  cor_long <- as.data.frame(as.table(cor_class_env))
  names(cor_long) <- c("Class", "Environmental_Variable", "Correlation")
  
  # Incluir información de p-value ajustado
  p_long <- as.data.frame(as.table(p_class_env_adjusted))
  names(p_long) <- c("Class", "Environmental_Variable", "P_Value_Adjusted")
  cor_long <- merge(cor_long, p_long, by = c("Class", "Environmental_Variable"))
  
  # Filtrar correlaciones significativas con FDR
  cor_long_filtered <- cor_long[cor_long$P_Value_Adjusted < 0.05, ]
  
  # Si no hay correlaciones significativas, mostrar todas con nota
  if(nrow(cor_long_filtered) == 0) {
    cor_long_filtered <- cor_long
    subtitle <- paste0("Top ", n_avail, " most abundant classes (no significant correlations after FDR correction)")
  } else {
    subtitle <- paste0("Top ", n_avail, " most abundant classes (FDR-corrected, p.adj < 0.05)")
  }
  
  # Redondear los valores para mejor visualización
  cor_long_filtered$Correlation_rounded <- round(cor_long_filtered$Correlation, 2)
  
  # Crear heatmap con valores numéricos
  p <- ggplot(cor_long_filtered, aes(x = Environmental_Variable, y = Class, fill = Correlation)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = Correlation_rounded), 
              color = ifelse(abs(cor_long_filtered$Correlation) > 0.5, "white", "black"),
              size = 3) +
    scale_fill_gradient2(
      low = "#2166AC", 
      mid = "white", 
      high = "#B2182B", 
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson\nCorrelation"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(face = "bold", size = 12),
      panel.grid = element_blank(),
      legend.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    ) +
    labs(title = paste0("Correlation between classes from ", domain_name, " and environmental variables"),
         subtitle = subtitle,
         x = "Environmental Variables", y = "Class")
  
  return(p)
}

# Generar heatmaps para bacteria y archaea
heatmap_bacteria <- generate_domain_correlation_heatmap(ps, "Bacteria", 10)
heatmap_archaea <- generate_domain_correlation_heatmap(ps, "Archaea", 10)

# Mostrar los heatmaps
if(!is.null(heatmap_bacteria)) {
  print(heatmap_bacteria)
  ggsave("heatmap_correlacion_bacterias_FDR.pdf", heatmap_bacteria, width = 10, height = 8)
}

if(!is.null(heatmap_archaea)) {
  print(heatmap_archaea)
  ggsave("heatmap_correlacion_arqueas_FDR.pdf", heatmap_archaea, width = 10, height = 8)
}

# ============================================
# GUARDAR TABLAS CON P-VALORES AJUSTADOS
# ============================================
write.csv(cor_genus_env, "correlaciones_generos_ambiente.csv")
write.csv(p_genus_env, "pvalores_SIN_correccion.csv")
write.csv(p_genus_env_adjusted, "pvalores_CON_correccion_FDR.csv")

cat("\n¡Análisis completado! Se guardaron los archivos con y sin corrección para comparar.\n")
