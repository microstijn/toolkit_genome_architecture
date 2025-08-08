#!/usr/bin/env julia

#-----------------------------------------------------------------
#   Description:      Calculate gene overlap from GFF files.
#   Author:           SHP
#   Date:             2025
#   Revised:          2025-08-08
#-----------------------------------------------------------------

#-----------------------------------------------------------------
# preamble
#-----------------------------------------------------------------
using Pkg
Pkg.activate(".") # Activate environment in current directory

using ArgParse
using CSV
using DataFrames
using GFF3 # Corrected module name
using GenomicFeatures
using IntervalTrees
using Glob

#-----------------------------------------------------------------
# Helper Function for Genome Structure Analysis
#-----------------------------------------------------------------

"""
    find_gff_files_recursive(dir_path::String)

Recursively finds all GFF files in a directory and its subdirectories.
"""
function find_gff_files_recursive(dir_path::String)
    gff_files = String[]
    for (root, dirs, files) in walkdir(dir_path)
        for file in files
            if endswith(lowercase(file), ".gff3") || endswith(lowercase(file), ".gff")
                push!(gff_files, normpath(joinpath(root, file)))
            end
        end
    end
    return gff_files
end


"""
    get_genome_structure(file_paths::Vector{String})

Analyzes a list of GFF files to calculate gene overlaps and other metrics.

Args:
- file_paths: A vector of strings, where each string is a path to a GFF file.

Returns:
A DataFrame containing the calculated genome architecture metrics for each file.
"""
function get_genome_structure(file_paths::Vector{String})
    # Pre-allocate vectors for efficiency
    gff_files = String[]
    contig_counts = Union{Missing, Int}[]
    gene_counts = Union{Missing, Int}[]
    mean_overlap_lengths = Union{Missing, Float64}[]
    total_overlap_lengths = Union{Missing, Int}[]
    num_overlaps = Union{Missing, Int}[]

    for file_path in file_paths
        println("  Processing $file_path...")
        
        try
            # Check if the file is empty before trying to read it
            if stat(file_path).size == 0
                @warn "Skipping empty file: $file_path"
                push!(gff_files, basename(file_path))
                push!(contig_counts, missing)
                push!(gene_counts, missing)
                push!(mean_overlap_lengths, missing)
                push!(total_overlap_lengths, missing)
                push!(num_overlaps, missing)
                continue
            end

            # Read GFF file
            gff_records = GFF3.read(file_path)
            
            # Group records by sequence ID (contig/chromosome)
            grouped_records = Dict{String, Vector{GFF3.Record}}()
            for record in gff_records
                contig_name = GFF3.seqname(record)
                if !haskey(grouped_records, contig_name)
                    grouped_records[contig_name] = GFF3.Record[]
                end
                push!(grouped_records[contig_name], record)
            end

            contig_count = length(keys(grouped_records))
            
            # Initialize metrics for this file
            total_genes = 0
            file_total_overlap_length = 0
            file_num_overlaps = 0
            
            for (contig_name, records) in grouped_records
                # Filter for gene features
                genes = filter(record -> GFF3.featuretype(record) == "gene", records)
                total_genes += length(genes)

                # Build an IntervalTree for efficient overlap searching
                gene_intervals = GenomicFeatures.Interval{Nothing}[]
                for gene in genes
                    push!(gene_intervals, GenomicFeatures.Interval(GFF3.seqname(gene), GFF3.start(gene), GFF3.stop(gene)))
                end
                tree = IntervalTree(gene_intervals)

                # Find overlaps for each gene
                for gene_a in genes
                    # Exclude self-overlap
                    overlaps = eachoverlap(tree, GenomicFeatures.Interval(GFF3.seqname(gene_a), GFF3.start(gene_a), GFF3.stop(gene_a)))
                    for gene_b_interval in overlaps
                        if gene_a != gene_b_interval
                            overlap_start = max(GFF3.start(gene_a), leftposition(gene_b_interval))
                            overlap_stop = min(GFF3.stop(gene_a), rightposition(gene_b_interval))
                            
                            overlap_length = overlap_stop - overlap_start + 1
                            if overlap_length > 0
                                file_total_overlap_length += overlap_length
                                file_num_overlaps += 1
                            end
                        end
                    end
                end
            end

            mean_overlap_length = file_num_overlaps > 0 ? file_total_overlap_length / file_num_overlaps : 0.0

            # Push results to vectors
            push!(gff_files, basename(file_path))
            push!(contig_counts, contig_count)
            push!(gene_counts, total_genes)
            push!(mean_overlap_lengths, mean_overlap_length)
            push!(total_overlap_lengths, file_total_overlap_length)
            push!(num_overlaps, file_num_overlaps)

        catch e
            @error "Failed to process $file_path" exception=(e, catch_backtrace())
            # Add placeholders for failed files
            push!(gff_files, basename(file_path))
            push!(contig_counts, missing)
            push!(gene_counts, missing)
            push!(mean_overlap_lengths, missing)
            push!(total_overlap_lengths, missing)
            push!(num_overlaps, missing)
        end
    end

    # Create a DataFrame from the collected data
    return DataFrame(
        file_name = gff_files,
        contigs = contig_counts,
        genes = gene_counts,
        total_overlap_length = total_overlap_lengths,
        num_overlaps = num_overlaps,
        mean_overlap_length = mean_overlap_lengths
    )
end

#-----------------------------------------------------------------
# Argument Parsing
#-----------------------------------------------------------------
function parse_commandline()
    s = ArgParseSettings(description="Calculate genome architecture metrics from GFF files.")
    @add_arg_table! s begin
        "--inputdir", "-i"
            help = "Path to the input directory containing GFF files"
            required = true
        "--output", "-o"
            help = "Path for the output CSV file"
            required = true
    end
    return parse_args(s)
end


#-----------------------------------------------------------------
# Main Logic
#-----------------------------------------------------------------
function main()
    args = parse_commandline()
    
    input_dir = normpath(args["inputdir"])
    output_path = normpath(args["output"])

    if !isdir(input_dir)
        @error "Input directory not found: $input_dir"
        return
    end

    # Use the new recursive function to find all GFF files
    file_paths = find_gff_files_recursive(input_dir)
    if isempty(file_paths)
        @warn "No GFF files found in $input_dir. Exiting."
        return
    end

    println("Searching for GFF files in $input_dir...")
    println("Found $(length(file_paths)) GFF files. Starting analysis...")

    # Calculate genome structure metrics
    df_metrics = get_genome_structure(file_paths)

    # Write results to output file
    println("Writing results to: $output_path")
    CSV.write(output_path, df_metrics)

    println("Done.")
end

main()
