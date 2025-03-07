# Load the packages to verify installation
library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)
library(viridis)
library(dendextend)
library(tidyr)
library(dplyr)

tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T, row.names=1, check.names=F, sep="\t"))
print("Successfully read taxonomy file")
sample_info_tab <- read.csv("pantano_metadata.csv", header=T, row.names=1, check.names=F)
print("Successfully read metadata file")
count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1, check.names=F, sep="\t")[ , -c(7:9)]
print("Successfully read counts file")

colnames(count_tab) <- gsub("reads/", "", colnames(count_tab))

count_tab <- count_tab + 1 # Because of the zeros

deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~mes) # Cambiarlo para que se introduzca desde param_file
deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

vst_trans_count_tab <- assay(deseq_counts_vst)
euc_dist <- dist(t(vst_trans_count_tab)) #Decidir también por el tipo de distancia para la corrida.
euc_clust <- hclust(euc_dist, method="ward.D2")

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_tab$color_lago[order.dendrogram(euc_dend)]) # Cambiar esto por variables de columna de color
labels_colors(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.")

vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# Generate Ordinate Plot
method_pord <- "MDS"
vst_pcoa <- ordinate(vst_physeq, method=method_pord, distance="bray")
eigen_vals <- vst_pcoa$values$Eigenvalues

otu_table(physeq) <- otu_table(decostand(otu_table(physeq), method="hellinger"), 
                               taxa_are_rows=TRUE)
bray_pcoa <- ordinate(physeq, method="PCoA", distance="bray")

# Set color for ordinate plot with parameter list
set_color <- "lago"
col_color <- "color_lago"

pord1 <- plot_ordination(vst_physeq, vst_pcoa, color = "lago") + 
  geom_point(size = 1) + 
  labs(col = "type") + 
  geom_text(aes(label = rownames(sample_info_tab)), hjust = 0.3, vjust = -0.4) + 
  coord_fixed(ratio = sqrt(eigen_vals[2] / eigen_vals[1])) + 
  ggtitle(paste("Ordinates Graphic -", method_pord)) +  # Aquí incluimos el método
  scale_color_manual(values = unique(sample_info_tab$color_lago[order(sample_info_tab$lago)])) + 
  theme_bw() + 
  theme(legend.position = "none")

pord1

####
# Create phyloseq object from your imported data
physeq <- phyloseq(otu_table(count_tab, taxa_are_rows=TRUE),
                   tax_table(tax_tab),
                   sample_data(sample_info_tab))

# Apply Hellinger transformation
# Extract the OTU table from your phyloseq object
otu <- otu_table(physeq)
# Apply Hellinger transformation using vegan's decostand function
hell_transformed <- decostand(otu, method="hellinger")
# Replace the OTU table in your phyloseq object with the transformed data
otu_table(physeq) <- otu_table(hell_transformed, taxa_are_rows=TRUE)

# Now perform ordination with Bray-Curtis distance
hellinger_bray_pcoa <- ordinate(physeq, method="PCoA", distance="bray")

# To visualize the results
pord1 <- plot_ordination(physeq, hellinger_bray_pcoa, color = "lago") + 
  geom_point(size = 1) + 
  labs(col = "type") + 
  geom_text(aes(label = rownames(sample_info_tab)), hjust = 0.3, vjust = -0.4) + 
  coord_fixed(ratio = sqrt(eigen_vals[2] / eigen_vals[1])) + 
  ggtitle(paste("Ordinates Graphic - PCoA")) +  # Aquí incluimos el método
  scale_color_manual(values = unique(sample_info_tab$color_lago[order(sample_info_tab$lago)])) + 
  theme_bw() + 
  theme(legend.position = "none")

pord1
####