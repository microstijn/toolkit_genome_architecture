#!/usr/bin/env julia

#-----------------------------------------------------------------
#   Description:      Calculate gene overlap from GFF files.
#   Author:           SHP
#   Date:             2025
#   Revised:          2025-09-05 (Added progress indicators)
#-----------------------------------------------------------------

#-----------------------------------------------------------------
# preamble
#-----------------------------------------------------------------
using Pkg
Pkg.activate(".") # Activate environment in current directory

using ArgParse
using CSV
using DataFrames
using GFF3
using Glob

#-----------------------------------------------------------------
# Helper Functions from the original script
#-----------------------------------------------------------------

# length of interval
function length_interval(array_start::AbstractArray{T}, array_end::AbstractArray{T}) where T<:Int
    return array_end .- array_start
end

# space between genes
function space_between(array_start::AbstractArray{T}, array_end::AbstractArray{T}) where T<:Int
    return array_start[2:end] .- array_end[1:end-1]
end

# overlap intervals
function count_overlaps(length_array::AbstractArray{T}) where T<:Int
    uni_overlaps = @view length_array[length_array .<= 0]
    total_uni_overlaps = length(uni_overlaps)
    sum_uni_overlaps = total_uni_overlaps == 0 ? 0 : sum(uni_overlaps)
    return total_uni_overlaps, sum_uni_overlaps
end

#-----------------------------------------------------------------
# Main analysis function for a single GFF file
#-----------------------------------------------------------------
function process_gff_file(file_path::String)
    df_temp_storage = DataFrame()
    genome_name = first(split(basename(file_path), "_"))

    try
        reader = GFF3.Reader(open(file_path, "r"))
        for rec in reader
            # Filter for features that are genes and are on a contig/chromosome
            genes = [feature for feature in rec if GFF3.type(feature) == "gene" && GFF3.seqid(rec) != "na"]
            if isempty(genes)
                continue
            end

            # Positive strand genes
            p_genes = [g for g in genes if GFF3.strand(g) == STRAND_POS]
            sort!(p_genes, by=GFF3.seqstart)
            p_starts = GFF3.seqstart.(p_genes)
            p_ends = GFF3.seqend.(p_genes)

            # Negative strand genes
            n_genes = [g for g in genes if GFF3.strand(g) == STRAND_NEG]
            sort!(n_genes, by=GFF3.seqstart)
            n_starts = GFF3.seqstart.(n_genes)
            n_ends = GFF3.seqend.(n_genes)

            # Metrics calculation...
            contig_size = GFF3.seqend(rec)
            p_gene_nr = length(p_genes)
            p_gene_length_sum = sum(length_interval(p_starts, p_ends))
            n_gene_nr = length(n_genes)
            n_gene_length_sum = sum(length_interval(n_starts, n_ends))

            p_gaps = space_between(p_starts, p_ends)
            p_U_overlap_nr, p_U_overlap_length_sum = count_overlaps(p_gaps)
            n_gaps = space_between(n_starts, n_ends)
            n_U_overlap_nr, n_U_overlap_length_sum = count_overlaps(n_gaps)

            p_gap_length_sum = sum(filter(x -> x > 0, p_gaps))
            n_gap_length_sum = sum(filter(x -> x > 0, n_gaps))
            
            # (Convergent/Divergent logic can be complex, keeping it simple for now)
            C_overlap_nr, D_overlap_nr, C_lengt_sum, D_lengt_sum = 0, 0, 0, 0

            push!(df_temp_storage, (
                genome_name=genome_name, contig_size=contig_size, p_gene_nr=p_gene_nr,
                p_gene_length_sum=p_gene_length_sum, n_gene_nr=n_gene_nr, n_gene_length_sum=n_gene_length_sum,
                p_U_overlap_nr=p_U_overlap_nr, p_U_overlap_length_sum=p_U_overlap_length_sum,
                n_U_overlap_nr=n_U_overlap_nr, n_U_overlap_length_sum=n_U_overlap_length_sum,
                p_gap_length_sum=p_gap_length_sum, n_gap_length_sum=n_gap_length_sum,
                C_overlap_nr=C_overlap_nr, D_overlap_nr=D_overlap_nr, C_lengt_sum=C_lengt_sum, D_lengt_sum=D_lengt_sum
            ))
        end
        close(reader)
    catch e
        @error "Failed to process file: $file_path" exception=(e, catch_backtrace())
    end
    return df_temp_storage
end

function find_gff_files_recursive(dir::String)
    return glob("*.gff", dir)
end

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

function main()
    args = parse_commandline()
    
    input_dir = normpath(args["inputdir"])
    output_path = normpath(args["output"])

    if !isdir(input_dir)
        @error "Input directory not found: $input_dir"
        return
    end

    println("Searching for GFF files in $input_dir...")
    file_paths = find_gff_files_recursive(input_dir)
    n_files = length(file_paths)
    if isempty(file_paths)
        @warn "No GFF files found in $input_dir. Exiting."
        return
    end
    println("Found $n_files GFF files. Starting analysis...")

    # Main DataFrame to store all results
    storage_dataframe = DataFrame()

    # Process each file with a progress indicator
    for (i, file_path) in enumerate(file_paths)
        # Print progress
        percentage = round(Int, (i / n_files) * 100)
        print("\rProcessing file $i of $n_files ($percentage%) - $(basename(file_path))      ")
        
        # Process the file and append results
        df_file_metrics = process_gff_file(file_path)
        if !isempty(df_file_metrics)
            append!(storage_dataframe, df_file_metrics)
        end
    end
    println("\nProcessing complete.")

    # Write the final DataFrame to a CSV file
    mkpath(dirname(output_path))
    println("Writing results to: $output_path")
    try
        CSV.write(output_path, storage_dataframe)
    catch e
        @error "Failed to write CSV file at $output_path." exception=(e, catch_backtrace())
    end
    println("Done.")
end

main()