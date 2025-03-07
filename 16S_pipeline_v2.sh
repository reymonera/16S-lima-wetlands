#!/bin/bash

# Activating conda - use a more flexible approach
eval "$(conda shell.bash hook)"
conda activate dada2_env || { echo "Failed to activate conda environment"; exit 1; }

# Exit on error
set -e

# Function to display usage information
usage() {
    echo "Usage: $0 -i <input_dir> -d <database_file> -o <output_dir> -t <threads> -p <params_file> -m <metadata_file>"
    echo "  -i  Directory containing input FASTQ files (required)"
    echo "  -d  Path to the compressed database file (required)"
    echo "  -o  Directory to save output files (default: './output')"
    echo "  -t  Number of threads to use (default: 4)"
    echo "  -p  Parameter file containing variable sets (required)"
    echo "  -m  Metadata file containing sample data (required)"
    exit 1
}

# Default values
output_dir="./output"
threads=4

# Parse command-line argument
while getopts ":i:d:o:t:p:m:" opt; do
    case $opt in
        i) input_dir="$OPTARG" ;;
        d) database="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        t) threads="$OPTARG" ;;
        p) params_file="$OPTARG" ;;
        m) metadata_file="$OPTARG" ;;
        *) usage ;;
    esac
done

# Ensure required arguments are provided
if [[ -z "$input_dir" || -z "$database" || -z "$params_file"  || -z "$metadata_file" ]]; then
    usage
fi

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

ls $input_dir*_R1_001.fastq.gz | cut -f1 -d "_" > $output_dir/samples
cp $database $output_dir
cp $metadata_file $output_dir

# Pre-create run directories
num_runs=$(grep -c "overlap=" "$params_file")
for ((i=1; i<=$num_runs; i++)); do
    mkdir -p "$output_dir/run_$i"
done

# Generate R script for DADA2 and Phyloseq analysis
r_script="$output_dir/dada2_analysis.R"
cat <<'EOF' > "$r_script"
library(dada2)
library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)
library(viridis)
library(dendextend)
library(decontam)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
db_file <- args[3]
threads <- as.numeric(args[4])
params_file <- args[5]
metadata_file <- args[6]

# Store original directory
original_dir <- getwd()
#print(original_dir)
# Read parameter sets from file
params_text <- readLines(params_file)
param_sets <- list()
current_set <- character()

# Parse parameter sets
for (line in params_text) {
    if (line == "") {
        if (length(current_set) > 0) {
            param_sets[[length(param_sets) + 1]] <- current_set
            current_set <- character()
        }
    } else {
        current_set <- c(current_set, line)
    }
}
# Add the last set if it exists
if (length(current_set) > 0) {
    param_sets[[length(param_sets) + 1]] <- current_set
}

# Function to parse parameter set
parse_params <- function(param_set) {
    params <- list()
    for (line in param_set) {
        parts <- strsplit(line, "=")[[1]]
        key <- trimws(parts[1])
        value <- trimws(parts[2])
        params[[key]] <- value
    }
    return(params)
}

# Setup initial data
samples <- scan(file.path(output_dir, "samples"), what="character")
print(file.path(output_dir, "samples"))
print(samples)
#forward_reads <- file.path(input_dir, paste0(samples, "_S1_L002_R1_001.fastq.gz"))
#reverse_reads <- file.path(input_dir, paste0(samples, "_S1_L002_R2_001.fastq.gz"))

#samples <- scan("samples", what="character")
forward_reads <- paste0(original_dir, "/", samples, "_S1_L002_R1_001.fastq.gz")
reverse_reads <- paste0(original_dir, "/", samples, "_S1_L002_R2_001.fastq.gz")

# Generate quality plots (only once)
sample_count <- length(forward_reads)
group_size <- ceiling(sample_count / 4)

cat("\nGenerating quality plots for each sample...\n")

# Initial quality plots in main output directory
setwd(output_dir)

forward_quality_files <- c()
reverse_quality_files <- c()

for (i in seq(1, sample_count, 4)) {
  if (sample_count <= group_size) {
    start = i
    end = sample_count
  } else if (sample_count >= group_size) {
    start = i
    end = sample_count
  } else {
    start = i
    end = i + 3
  }
  
  forward_quality_files <- c(forward_quality_files, paste0("forward_quality_", start, "-", end, ".pdf"))
  reverse_quality_files <- c(reverse_quality_files, paste0("reverse_quality_", start, "-", end, ".pdf"))
}

if (!all(file.exists(forward_quality_files)) || !all(file.exists(reverse_quality_files))) {
    for (i in seq(1, sample_count, 4)) {
        if (sample_count <= group_size) {
            start = i
            end = sample_count
        } else if (sample_count >= group_size) {
            start = i
            end = sample_count
        } else {
            start = i
            end = i + 3
        }
        
        pdf(paste0("forward_quality_", start,"-", end,".pdf"))
        print(plotQualityProfile(forward_reads[start:end]))
        dev.off()
        
        pdf(paste0("reverse_quality_", start,"-", end, ".pdf"))
        print(plotQualityProfile(reverse_reads[start:end]))
        dev.off()
    }
}

# Process each parameter set
for (run_index in seq_along(param_sets)) {
    cat(sprintf("\nProcessing run %d of %d\n", run_index, length(param_sets)))
    
    # Use pre-created run directory
    run_dir <- file.path(sprintf("run_%d", run_index))
    cat(sprintf("\nUsing directory: %s\n", run_dir))
    
    if (!dir.exists(run_dir)) {
        stop(sprintf("Run directory does not exist: %s", run_dir))
    }
    
    tryCatch({
        setwd(run_dir)
    }, error = function(e) {
        stop(sprintf("Failed to change to directory %s: %s", run_dir, e$message))
    })
    
    # Parse parameters for this run
    params <- parse_params(param_sets[[run_index]])

    # Set up filtered read paths for this run
    filtered_forward_reads <- file.path(run_dir, paste0(basename(samples), "_sub_R1_filtered.fq.gz"))
    filtered_reverse_reads <- file.path(run_dir, paste0(basename(samples), "_sub_R2_filtered.fq.gz"))

    # Check if filtered reads DON'T already exist
    if (!all(file.exists(filtered_forward_reads)) || !all(file.exists(filtered_reverse_reads))) {
        # Get parameters
        overlap <- as.numeric(params$overlap)
        truncLen <- c(as.numeric(params$trunclen_fr), as.numeric(params$trunclen_rev))
        cat(sprintf("\nParameters for run %d:\n", run_index))
        cat(sprintf("Overlap: %d\n", overlap))
        cat(sprintf("TruncLen Forward: %d, Reverse: %d\n", truncLen[1], truncLen[2]))
        
        # Filter and trim
        filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                                    reverse_reads, filtered_reverse_reads,
                                    maxN=0, maxEE=c(2,2), rm.phix=TRUE,
                                    minLen=175, truncLen=truncLen)

        # Quality plots only if we just created new filtered reads
        for (i in seq(1, sample_count, 4)) {
            if (sample_count <= group_size) {
                start = i
                end = sample_count
            } else if (sample_count >= group_size) {
                start = i
                end = sample_count
            } else {
                start = i
                end = i + 3
            }
            pdf(paste0("filtered_forward_quality_", start,"-", end,".pdf"))
            print(plotQualityProfile(filtered_forward_reads[start:end]))
            dev.off()

            pdf(paste0("filtered_reverse_quality_", start,"-", end, ".pdf"))
            print(plotQualityProfile(filtered_reverse_reads[start:end]))
            dev.off()
        }
    } else {
        cat("Filtered reads already exist, skipping filtering and quality plotting steps...\n")
    }

    # Check if final output files exist
    if (!file.exists("ASVs_counts-no-contam.tsv") || !file.exists("ASVs_taxonomy-no-contam.tsv")) {
        # Learn errors
        err_forward_reads <- learnErrors(filtered_forward_reads, multithread=threads)
        err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=threads)
        
        # Generate error plots
        pdf("error_forward_plot.pdf")
        print(plotErrors(err_forward_reads, nominalQ=TRUE))
        dev.off()
        pdf("error_reverse_plot.pdf")
        print(plotErrors(err_reverse_reads, nominalQ=TRUE))
        dev.off()
        
        # Save run parameters
        write.table(data.frame(
            parameter = names(params),
            value = unlist(params)
        ), "run_parameters.txt", row.names = FALSE, quote = FALSE, sep = "\t")
        
        # Derreplication
        derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
        names(derep_forward) <- samples
        derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
        names(derep_reverse) <- samples
        
        # Inferring ASVs with DADAs algorithm
        pool <- as.character(params$pool)
        dada_forward <- dada(derep_forward, err=err_forward_reads, pool=pool, multithread=threads)
        dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool=pool, multithread=threads)
        
        # Merge amplicons
        merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, trimOverhang=TRUE, minOverlap=overlap)
        seqtab <- makeSequenceTable(merged_amplicons)
        
        # Check how much you have
        cat(sprintf("\nCounts: %d : %d\n", dim(seqtab)[1], dim(seqtab)[2]))
        seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)
        cat(sprintf("\nPercentage left after chimera removal: %.2f%%\n", (sum(seqtab.nochim) / sum(seqtab)) * 100))
        
        # Summary table
        getN <- function(x) sum(getUniques(x))
        summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1], filtered=filtered_out[,2],
            dada_f=sapply(dada_forward, getN), dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
            nonchim=rowSums(seqtab.nochim), final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
        cat("\n=== Summary Table ===\n")
        print(format(summary_tab, justify = "left"))
        
        # Assigning taxa
        minBoot <- as.numeric(params$minBoot)
        db_tax_path <- file.path("..", db_file)
        at_taxa <- assignTaxonomy(seqtab.nochim, db_tax_path, minBoot=minBoot, multithread=threads)
        
        # Assigning ASV names
        asv_seqs <- colnames(seqtab.nochim)
        asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
        for (i in 1:dim(seqtab.nochim)[2]) {
            asv_headers[i] <- paste(">ASV", i, sep="_")
        }
        
        # Write output files
        asv_fasta <- c(rbind(asv_headers, asv_seqs))
        write(asv_fasta, "ASVs.fa")
        
        asv_tab <- t(seqtab.nochim)
        row.names(asv_tab) <- sub(">", "", asv_headers)
        write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)
        
        ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
        colnames(at_taxa) <- ranks
        at_taxa[startsWith(at_taxa, "unclassified_")] <- NA
        rownames(at_taxa) <- gsub(pattern=">", replacement="", x=asv_headers)
        write.table(at_taxa, "ASVs_taxonomy.tsv", sep="\t", quote=FALSE, col.names=NA)

        # Decontamination
        colnames(asv_tab) # en estas muestras, los blancos son los primeros 4

        # Hay que hacer una forma de que se introduzcan los nombres de los blancos, de modo tal que se reordene y se elimine.

        vector_for_decontam <- c(rep(TRUE, 6), rep(FALSE, 3), rep(TRUE, 19))
        contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)
        table(contam_df$contaminant)

        contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

        contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
        dont_want <- sort(c(contam_indices, contam_indices + 1))

        asv_fasta_no_contam <- asv_fasta[- dont_want]
        asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]
        asv_tax_no_contam <- at_taxa[!row.names(at_taxa) %in% contam_asvs, ]

        write(asv_fasta_no_contam, "ASVs-no-contam.fa")
        write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv", sep="\t", quote=F, col.names=NA)
        write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv", sep="\t", quote=F, col.names=NA)

    } else {
        cat("Final ASV tables already exist, skipping DADA2 processing steps...\n")
    }

    # Assigning it to variables for Phyloseq
    print("Asignando a tabs...")
    tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T, row.names=1, check.names=F, sep="\t"))
    print("Successfully read taxonomy file")
    sample_info_tab <- read.csv(file.path("..", metadata_file), header=T, row.names=1, check.names=F)
    print("Successfully read metadata file")
    count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1, check.names=F, sep="\t")[ , -c(7:9)]
    print("Successfully read counts file")

    print(sample_info_tab)

    print("Ya se guardó en tab variables")
    print(input_dir)
    colnames(count_tab) <- gsub(input_dir, "", colnames(count_tab))

    print(head(count_tab))

    count_tab <- count_tab + 1 # Because of the zeros

    # DESeq normalization
    deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~mes) # Cambiarlo para que se introduzca desde param_file
    deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

    print("Ya pasó DeSeq")

    vst_trans_count_tab <- assay(deseq_counts_vst)
    euc_dist <- dist(t(vst_trans_count_tab)) #Decidir también por el tipo de distancia para la corrida.
    euc_clust <- hclust(euc_dist, method="ward.D2")

    print("Ya pasó euc_clust")

    # Transforming it into a dendogram
    euc_dend <- as.dendrogram(euc_clust, hang=0.1)
    dend_cols <- as.character(sample_info_tab$color_lago[order.dendrogram(euc_dend)]) # Cambiar esto por variables de columna de color
    labels_colors(euc_dend) <- dend_cols

    pdf("euclidean_distance.pdf")
    print(plot(euc_dend, ylab="VST Euc. dist."))
    dev.off()

    print("Ya pasó la euc gráfica")

    vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
    sample_info_tab_phy <- sample_data(sample_info_tab)
    vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

    # Generate Ordinate Plot
    method_pord <- as.character(params$method_pord)
    vst_pcoa <- ordinate(vst_physeq, method=method_pord, distance="euclidean")
    eigen_vals <- vst_pcoa$values$Eigenvalues

    # Set color for ordinate plot with parameter list
    set_color <- as.character(params$set_color)
    col_color <- as.character(params$set_col_color)

    pord1 <- plot_ordination(vst_physeq, vst_pcoa, color = "lago") + 
        geom_point(size = 1) + 
        labs(col = "type") + 
        geom_text(aes(label = rownames(sample_info_tab)), hjust = 0.3, vjust = -0.4) + 
        coord_fixed(ratio = sqrt(eigen_vals[2] / eigen_vals[1])) + 
        ggtitle(paste("Ordinates Graphic -", method_pord)) +  # Aquí incluimos el método
        scale_color_manual(values = unique(sample_info_tab$color_lago[order(sample_info_tab$lago)])) + 
        theme_bw() + 
        theme(legend.position = "none")

    pdf("beta_ordination.pdf")
    print(pord1)
    dev.off()

    # Set variables for alpha-diversity
    count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
    tax_tab_phy <- tax_table(tax_tab)
    ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

    alpha_x <- as.character(params$alpha_x)
    alpha_color <- as.character(params$alpha_color)
    alpha_col_color <- as.character(params$alpha_col_color)
    alpha_measure_1 <- as.character(params$alpha_measure_1)
    alpha_measure_2 <- as.character(params$alpha_measure_2)

    palpha1 <- plot_richness(ASV_physeq, x = alpha_x, color = alpha_color, 
                         measures = c(alpha_measure_1, alpha_measure_2)) + 
    scale_color_manual(values = unique(sample_info_tab[[alpha_col_color]][order(sample_info_tab[[alpha_color]])])) + 
    theme_bw() + 
    theme(legend.title = element_blank(), 
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    pdf("alpha_diversity.pdf")
    print(palpha1)
    dev.off()

    tax_category <- as.character(params$tax_category)
    tax_category_color <- as.character(params$tax_category_color)

    if(exists("params") && !is.null(params$tax_min_abundance)) {
        tax_min_abundance <- as.numeric(params$tax_min_abundance)
    } else {
        tax_min_abundance <- 5
    }

    if(exists("params") && !is.null(params$tax_rank_specific)) {
        tax_rank_specific <- as.character(params$tax_category_color)
    } else {
        tax_rank_specific <- NULL
    }

    ### ============================ MOST ABUNDANT FUNCTION =================================================
    plot_taxa_by_rank <- function(tax_rank = "phylum", min_abundance = tax_min_abundance, 
                             special_group = NULL, count_tab = NULL, physeq_obj = NULL,
                             break_special_group = TRUE,
                             tax_category_param = tax_category,
                             tax_category_color_param = tax_category_color) {
    
        # Check for global objects first
        if(is.null(physeq_obj) && exists("ASV_physeq", envir = .GlobalEnv)) {
            physeq_obj <- get("ASV_physeq", envir = .GlobalEnv)
        }
        
        if(is.null(count_tab) && exists("count_tab", envir = .GlobalEnv)) {
            count_tab <- get("count_tab", envir = .GlobalEnv)
        }
        
        # Validate we have required objects after checking globals
        if(is.null(physeq_obj)) {
            stop("No phyloseq object found. Please provide 'physeq_obj' or make sure 'ASV_physeq' exists.")
        }
        
        # Get count_tab from physeq_obj if not provided
        if(is.null(count_tab)) {
            message("Using count table from the phyloseq object.")
            count_tab <- as(otu_table(physeq_obj), "matrix")
        }
        
        # Validate taxonomic rank
        valid_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
        if(!tax_rank %in% valid_ranks) {
            stop(paste("tax_rank must be one of:", paste(valid_ranks, collapse=", ")))
        }
        
        # Create taxa counts at the specified rank
        taxa_counts_tab <- otu_table(tax_glom(physeq_obj, taxrank = tax_rank))
        
        # Get taxonomy vector for the specified rank
        taxa_tax_vec <- as.vector(tax_table(tax_glom(physeq_obj, taxrank = tax_rank))[, tax_rank])
        rownames(taxa_counts_tab) <- taxa_tax_vec
        
        # Calculate unclassified sequences
        unclassified_tax_counts <- colSums(count_tab) - colSums(taxa_counts_tab)
        
        # Add unclassified to the table
        taxa_and_unidentified_counts_tab <- rbind(taxa_counts_tab, "Unclassified" = unclassified_tax_counts)
        
        # Handle special taxonomic group if specified (e.g., Proteobacteria) and break_special_group is TRUE
        if(!is.null(special_group) && special_group %in% rownames(taxa_and_unidentified_counts_tab) && break_special_group) {
            # First remove the special group, we'll add it back broken down
            temp_major_taxa_counts_tab <- taxa_and_unidentified_counts_tab[
            !rownames(taxa_and_unidentified_counts_tab) %in% special_group, ]
            
            # Find the next taxonomic rank
            next_rank_idx <- which(valid_ranks == tax_rank) + 1
            if(next_rank_idx <= length(valid_ranks)) {
            next_rank <- valid_ranks[next_rank_idx]
            
            # Get counts for the next rank
            next_rank_counts_tab <- otu_table(tax_glom(physeq_obj, taxrank = next_rank))
            
            # Get taxonomy for both current and next rank
            next_rank_tax_tab <- tax_table(tax_glom(physeq_obj, taxrank = next_rank))
            
            # Create vectors for both taxonomic levels
            current_rank_vec <- next_rank_tax_tab[, tax_rank]
            next_rank_vec <- next_rank_tax_tab[, next_rank]
            rows_tmp <- rownames(next_rank_tax_tab)
            
            # Create a data frame with both taxonomic levels
            next_rank_tax_df <- data.frame(
                "current_rank" = current_rank_vec, 
                "next_rank" = next_rank_vec, 
                row.names = rows_tmp
            )
            
            # Get children of the special group
            special_children_vec <- as.vector(
                next_rank_tax_df[next_rank_tax_df$current_rank == special_group, "next_rank"]
            )
            
            # Set next rank as row names
            rownames(next_rank_counts_tab) <- as.vector(next_rank_tax_df$next_rank)
            
            # Extract counts for the special group's children
            special_children_counts_tab <- next_rank_counts_tab[
                rownames(next_rank_counts_tab) %in% special_children_vec, ]
            
            # Calculate counts for members of special group not resolved to the next rank
            special_unresolved_counts <- taxa_and_unidentified_counts_tab[
                rownames(taxa_and_unidentified_counts_tab) %in% special_group, ] - 
                colSums(special_children_counts_tab)
            
            # Combine everything
            unresolved_rowname <- paste0("Unresolved_", special_group)
            combined_data <- rbind(temp_major_taxa_counts_tab, special_children_counts_tab)
            major_taxa_counts_tab <- rbind(combined_data, special_unresolved_counts)
            rownames(major_taxa_counts_tab)[nrow(major_taxa_counts_tab)] <- unresolved_rowname
            } else {
            # If we can't break down further, just use the original table
            major_taxa_counts_tab <- taxa_and_unidentified_counts_tab
            }
        } else {
            # If no special group handling, use the original table
            major_taxa_counts_tab <- taxa_and_unidentified_counts_tab
        }
        
        # Check that we didn't lose any counts
        if(!identical(colSums(major_taxa_counts_tab), colSums(count_tab))) {
            warning("Some counts may have been lost in the taxonomic aggregation process")
        }
        
        # Generate proportions table
        major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x) * 100)
        
        # Filter taxa below threshold
        temp_filt_major_taxa_proportions_tab <- data.frame(
            major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > min_abundance, ]
        )
        
        # Add "Other" category for filtered taxa
        filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
        filt_major_taxa_proportions_tab <- rbind(
            temp_filt_major_taxa_proportions_tab, 
            "Other" = filtered_proportions
        )
        
        # Prepare for plotting
        filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab
        filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- rownames(filt_major_taxa_proportions_tab_for_plot)
        
        # Transform to long format
        library(tidyr)
        library(dplyr)
        filt_major_taxa_proportions_tab_for_plot.g <- pivot_longer(
            filt_major_taxa_proportions_tab_for_plot, 
            !Major_Taxa, 
            names_to = "Sample", 
            values_to = "Proportion"
        ) %>% data.frame()
        
        # Get sample metadata
        sample_info_tab <- sample_data(physeq_obj)
        
        # Use the parameters instead of global variables
        sample_info_for_merge <- data.frame(
            "Sample" = rownames(sample_info_tab),
            "char" = sample_info_tab[[tax_category_param]],
            "color" = sample_info_tab[[tax_category_color_param]], 
            stringsAsFactors = FALSE
        )
        
        # Merge with metadata
        filt_major_taxa_proportions_tab_for_plot.g2 <- merge(
            filt_major_taxa_proportions_tab_for_plot.g, 
            sample_info_for_merge
        )
        
        # Create plot
        library(ggplot2)
        p <- ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Major_Taxa, Proportion)) +
            geom_jitter(aes(color = factor(char), shape = factor(char)), 
                        size = 2, width = 0.15, height = 0) +
            scale_color_manual(values = unique(
            filt_major_taxa_proportions_tab_for_plot.g2$color[
                order(filt_major_taxa_proportions_tab_for_plot.g2$char)
            ]
            )) +
            geom_boxplot(fill = NA, outlier.color = NA) + 
            theme_bw() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                legend.title = element_blank()) +
            labs(
            x = paste0(toupper(substr(tax_rank, 1, 1)), substr(tax_rank, 2, nchar(tax_rank))), 
            y = "% of 16S rRNA gene copies recovered", 
            title = paste("All samples -", tax_rank, "level")
            )
        
        # Return both the data and the plot
        return(list(
            data = filt_major_taxa_proportions_tab_for_plot.g2,
            plot = p,
            proportions_table = filt_major_taxa_proportions_tab
        ))
    }
    ### ============================ MOST ABUNDANT FUNCTION =================================================

    ptaxa1 <- plot_taxa_by_rank(
        tax_rank = "phylum", 
        min_abundance = tax_min_abundance,
        tax_category_param = tax_category,
        tax_category_color_param = tax_category_color
    )

    if(exists("params") && !is.null(params$tax_rank_specific)) {
        tax_rank <- as.character(params$tax_rank_specific)
        print(tax_rank)
        ptaxa2 <- plot_taxa_by_rank(
            tax_rank = tax_rank, 
            min_abundance = tax_min_abundance,
            tax_category_param = tax_category,
            tax_category_color_param = tax_category_color
        )

        pdf("taxa_abundance_specific.pdf")
        print(ptaxa2)
        dev.off()
        } 

    pdf("taxa_abundance.pdf")
    print(ptaxa1)
    dev.off()


    # Return to original directory
    setwd("..")
}

EOF

# Execute the R script with parameters
Rscript "$r_script" "$input_dir" "$output_dir" "$database" "$threads" "$params_file" "$metadata_file"


