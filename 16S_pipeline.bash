#!/bin/bash

# Activating conda
source ~/anaconda3/etc/profile.d/conda.sh  # Load Conda
conda activate dada2_env

# Exit on error
set -e

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_dir> -d <database_file> -o <output_dir> -t <threads> -v <overlap>"
    echo "  -i  Directory containing input FASTQ files (required)"
    echo "  -d  Path to the compressed database file (required)"
    echo "  -o  Directory to save output files (default: './output')"
    echo "  -t  Number of threads to use (default: 4)"
    echo "  -v  Minimum overlap required (default: 12)"
    exit 1
}

# Default values
output_dir="./output"
threads=4
overlap=12

# Parse command-line arguments
while getopts ":i:d:o:t:v:" opt; do
    case $opt in
        i) input_dir="$OPTARG" ;;
        d) database="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        v) overlap="$OPTARG" ;;
        *) usage ;;
    esac
done

# Ensure required arguments are provided
if [[ -z "$input_dir" || -z "$database" ]]; then
    usage
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Generate R script for DADA2 and Phyloseq analysis
r_script="$output_dir/dada2_pipeline_phyloseq.R"
cat <<EOF > "$r_script"
library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan)

# Input and output directories
input_dir <- "$input_dir"
output_dir <- "$output_dir"
db_file <- "$database"
threads <- $threads
overlap <- $overlap

# Ensure the database file exists
if (!file.exists(db_file)) {
    stop("Database file does not exist: ", db_file)
}

# List all FASTQ files in the input directory
fnFs <- sort(list.files(input_dir, pattern="_R1_.*fastq", full.names=TRUE))
fnRs <- sort(list.files(input_dir, pattern="_R2_.*fastq", full.names=TRUE))

# Divide samples into four equal groups for quality plots
sample_count <- length(fnFs)
group_size <- ceiling(sample_count / 4)

cat("\nGenerating quality plots for 4 sample groups...\n")
for (i in 1:4) {
  start_index <- ((i - 1) * group_size) + 1
  end_index <- min(i * group_size, sample_count)
  
  if (start_index <= sample_count) {
    fnFs_subset <- fnFs[start_index:end_index]
    fnRs_subset <- fnRs[start_index:end_index]
    
    pdf(file.path(output_dir, paste0("quality_group_", i, ".pdf")))
    plotQualityProfile(fnFs_subset)
    plotQualityProfile(fnRs_subset)
    dev.off()
  }
}
cat("Quality plots saved as PDFs for each group.\n")

# Prompt user for truncLen values
cat("\nInspect the quality profiles and choose truncLen values (e.g., 240,160):\n")
truncLen_input <- readline(prompt="Enter truncLen values (comma-separated): ")
truncLen <- as.integer(strsplit(truncLen_input, ",")[[1]])

# Filter and trim files
filtFs <- file.path(output_dir, basename(fnFs))
filtRs <- file.path(output_dir, basename(fnRs))
filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncLen, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=threads)

# Learn error rates
errF <- learnErrors(filtFs, multithread=threads)
errR <- learnErrors(filtRs, multithread=threads)

# Dereplicate sequences
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

# Sample names
sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

# Infer sequence variants
dadaFs <- dada(derepFs, err=errF, multithread=threads)
dadaRs <- dada(derepRs, err=errR, multithread=threads)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=overlap, verbose=TRUE)

# Make sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads, verbose=TRUE)

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, db_file, multithread=threads)

# Save results
saveRDS(seqtab.nochim, file.path(output_dir, "seqtab_nochim.rds"))
saveRDS(taxa, file.path(output_dir, "taxa.rds"))
write.csv(taxa, file.path(output_dir, "taxonomy.csv"))

# Interactive Phyloseq Analysis
cat("\nLoading data into Phyloseq...\n")
physeq <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows=FALSE),
  tax_table(taxa),
  sample_data(data.frame(SampleID=sample_names, row.names=sample_names))
)

# User interaction for distance index and color options
cat("\nAvailable distance indices: bray, jaccard, unifrac, euclidean\n")
distance_index <- readline(prompt="Choose a distance index for beta diversity (default: bray): ")
if (distance_index == "") distance_index <- "bray"

cat("\nEnter a string of colors (comma-separated) for plots (e.g., '#FF5733,#33FF57,#3357FF'):\n")
color_string <- readline(prompt="Enter custom colors: ")
custom_colors <- unlist(strsplit(color_string, ","))

# Alpha diversity
cat("\nGenerating alpha diversity plot...\n")
p_alpha <- plot_richness(physeq, x="SampleID", measures=c("Shannon", "Simpson")) +
  theme_minimal() +
  scale_color_manual(values=custom_colors)
ggsave(file.path(output_dir, "alpha_diversity.png"), p_alpha)

# Beta diversity (distance matrix and PCoA)
cat("\nCalculating beta diversity and generating PCoA plot...\n")
dist_matrix <- distance(physeq, method=distance_index)
pcoa <- ordinate(physeq, method="PCoA", distance=dist_matrix)
p_beta <- plot_ordination(physeq, pcoa, color="SampleID") +
  theme_minimal() +
  scale_color_manual(values=custom_colors)
ggsave(file.path(output_dir, "beta_diversity_pcoa.png"), p_beta)

# Taxonomy barplot
cat("\nGenerating taxonomy barplot...\n")
p_tax <- plot_bar(physeq, fill="Phylum") +
  theme_minimal() +
  scale_fill_manual(values=custom_colors)
ggsave(file.path(output_dir, "taxonomy_barplot.png"), p_tax)

cat("\nAnalysis complete. Results saved in output directory.\n")
EOF

# Print message
echo "Generated R script for DADA2 and Phyloseq pipeline at: $r_script"
echo "Run it in R to process your 16S data and generate interactive plots."