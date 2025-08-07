#!/usr/bin/env Rscript

#=====================================================
# Description:  Generate a phylogenetic distance matrix from a list of NCBI TaxIDs.
# Author:       SHP
# Date:         2025
# Revised:      2025-08-07
#=====================================================

#=====================================================
# Preamble: Package Management
#=====================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(argparse, phylostratr, data.table, ape, anicon)

#=====================================================
# Argument Parser
#=====================================================
parser <- ArgumentParser(description = "Generate a phylogenetic distance matrix from a list of NCBI TaxIDs.")

parser$add_argument("-i", "--input", 
                    type = "character", 
                    required = TRUE, 
                    help = "Input CSV file with a 'taxId' column.")

parser$add_argument("-o", "--output", 
                    type = "character", 
                    required = TRUE, 
                    help = "Path for the output CSV file containing the distance matrix.")

args <- parser$parse_args()

#=====================================================
# Main Logic
#=====================================================

# ---- 1. Input Validation ----
if (!file.exists(args$input)) {
  stop(paste("Input file not found:", args$input))
}

# ---- 2. Read Data ----
message("--> Reading taxon IDs from: ", args$input)
tryCatch({
  dat <- fread(args$input)
}, error = function(e) {
  stop(paste("Failed to read input file. Error:", e$message))
})

if (!"taxId" %in% names(dat)) {
  stop("Input file must contain a column named 'taxId'.")
}

# ---- 3. Build Tree ----
tax_ids <- unique(na.omit(dat$taxId))
message(paste("--> Building NCBI phylogenetic tree for", length(tax_ids), "unique taxa..."))

tryCatch({
  tree <- ncbi_tree(tax_ids)
}, error = function(e) {
  stop(paste("Failed to build the NCBI tree. Error:", e$message,
             "\nThis may be due to invalid taxon IDs or NCBI connectivity issues."))
})

# ---- 4. Calculate Distance ----
message("--> Calculating cophenetic distance matrix...")
tree_dist_matrix <- cophenetic.phylo(tree)

# ---- 5. Format and Save Output ----
message("--> Formatting and saving results to: ", args$output)

# Convert the matrix to a data.table
tree_dist_dt <- as.data.table(tree_dist_matrix, keep.rownames = "taxid_1")

# Melt the wide-format distance matrix to a long (tidy) format
tree_dist_long <- melt(tree_dist_dt, 
                       id.vars = "taxid_1", 
                       variable.name = "taxid_2", 
                       value.name = "cophenetic_distance")

# Remove self-comparisons
tree_dist_long <- tree_dist_long[taxid_1 != taxid_2]

# Save the formatted data
tryCatch({
  fwrite(tree_dist_long, args$output)
}, error = function(e) {
  stop(paste("Failed to write to output file. Error:", e$message))
})

message("✓ All done.")
