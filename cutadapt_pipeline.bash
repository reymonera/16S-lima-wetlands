#!/bin/bash

# Activate the Conda environment
source /home/marlen/miniforge3/etc/profile.d/conda.sh
conda init
conda activate dada2_env

# Function to display usage
usage() {
   echo "Usage: $0 -b BARCODE_TSV -d INPUT_DIRECTORY"
   exit 1
}

# Parse command-line arguments
while getopts ":b:d:" opt; do
   case $opt in
       b) BARCODE_TSV="$OPTARG"
          ;;
       d) INPUT_DIRECTORY="$OPTARG"
          ;;
       \?) echo "Invalid option: -$OPTARG" >&2
           usage
           ;;
       :) echo "Option -$OPTARG requires an argument." >&2
          usage
          ;;
   esac
done

# Check if required arguments are provided
[[ -z "$BARCODE_TSV" || -z "$INPUT_DIRECTORY" ]] && usage

mkdir -p nonprim_reads

# Process each sample
while IFS=$'\t' read -r sampleid barcode primer barcodename revprimer projectname description; do
   # Skip header line
   [[ "$sampleid" == "#SampleID" ]] && continue
   
   # Find matching R1 and R2 files
   R1=$(find "$INPUT_DIRECTORY" -name "*${sampleid}*R1*.fastq*")
   R2=$(find "$INPUT_DIRECTORY" -name "*${sampleid}*R2*.fastq*")
   
   # Skip if files not found
   # -g y -G se refieren a la orientación, y ambos son 5', así sea
   # el forward o el reverse. No se encontraron reversos complementarios
   # en el final de los reads, así que solo se quitó la primera
   # secuencia que hacía match en los reads que se procedieron a quitar.
   
   [[ -z "$R1" || -z "$R2" ]] && continue
   # -g "^$barcode" -G "^$barcode" \
   cutadapt -g "^$primer" -G "^$revprimer" \
            -o "nonprim_reads/$(basename "$R1")" \
            -p "nonprim_reads/$(basename "$R2")" \
            "$R1" "$R2"
done < "$BARCODE_TSV"