#!/usr/bin/env julia

#-----------------------------------------------------------------
#   Description:      Calculate gene overlap from GFF files.
#   Author:           SHP
#   Date:             2025
#   Revised:          2025-08-09
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
using GenomicFeatures
using IntervalTrees
using Glob

#-----------------------------------------------------------------
# Helper Functions from the original script
#-----------------------------------------------------------------

# length of interval
function length_interval(array_start::AbstractArray{T}, array_end::AbstractArray{T}) where T<:Int
    array_end .- array_start
end

# space between genes 
function space_between(array_start::AbstractArray{T}, array_end::AbstractArray{T}) where T<:Int
    array_start[2:end] .- array_end[1:end-1]
end

# overlap intervals
function count_overlaps(length_array::AbstractArray{T}) where T<:Int
    uni_overlaps = @view length_array[length_array .<= 0]
    total_uni_overlaps = length(uni_overlaps)
    if total_uni_overlaps == 0
        sum_uni_overlaps = 0
    else
        sum_uni_overlaps = sum(uni_overlaps)
    end
    return total_uni_overlaps, sum_uni_overlaps
end

function total_interval_length(length_array::AbstractArray{T}) where T<:Int
    n = length(length_array)
    if n == 0
        tot_len = 0
    else
        tot_len = sum(length_array)
    end
    return tot_len
end

function intervalify(s::AbstractArray{T}, e::AbstractArray{T}) where T<:Int
    [GenomicFeatures.Interval(string(i), a, b) for (i, (a, b)) in enumerate(zip(s, e))]
end

function bidirectional_overlaps(set1, set2)
    overlap_vect_divergent = Int64[]
    overlap_vect_convergent = Int64[]
    for overlap in eachoverlap(set1, set2)
        start_a = leftposition(overlap[1])
        start_b = leftposition(overlap[2])
        if start_a < start_b
            length_overlap_con = start_b - rightposition(overlap[1])
            push!(overlap_vect_convergent, length_overlap_con)
        else
            length_overlap_di = rightposition(overlap[2]) - start_a
            push!(overlap_vect_divergent, length_overlap_di)
        end
    end
    len_con = length(overlap_vect_convergent)
    len_di = length(overlap_vect_divergent)
    sum_vect_con = sum(overlap_vect_convergent)
    sum_vect_di = sum(overlap_vect_divergent)
    return len_con, len_di, sum_vect_con, sum_vect_di
end

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

#-----------------------------------------------------------------
# Main Logic
#-----------------------------------------------------------------

function get_genome_structure(file_paths::Vector{String})

    storage_dataframe = DataFrame(
        genome_name = String[],
        contig_size = Union{Missing, Int64}[],
        p_gene_nr = Union{Missing, Int64}[],
        p_gene_length_sum = Union{Missing, Int64}[],
        n_gene_nr = Union{Missing, Int64}[],
        n_gene_length_sum = Union{Missing, Int64}[],
        p_U_overlap_nr = Union{Missing, Int64}[],
        p_U_overlap_length_sum = Union{Missing, Int64}[],
        n_U_overlap_nr = Union{Missing, Int64}[],
        n_U_overlap_length_sum = Union{Missing, Int64}[],
        p_gap_length_sum = Union{Missing, Int64}[],
        n_gap_length_sum = Union{Missing, Int64}[],
        C_overlap_nr = Union{Missing, Int64}[],
        D_overlap_nr = Union{Missing, Int64}[],
        C_lengt_sum = Union{Missing, Int64}[],
        D_lengt_sum = Union{Missing, Int64}[]
    )
    
    n = length(file_paths)
    
    for (enum, file) in enumerate(file_paths)
        # Update and print a simple progress bar
        percentage = round(Int, (enum / n) * 100)
        progress_bar = "[" * repeat("=", percentage) * repeat(" ", 100 - percentage) * "]"
        print("\rProgress: $progress_bar $percentage% (File $enum of $n)")

        try
            # Get the unique genome name from the directory path
            path_parts = splitpath(file)
            genome_name = path_parts[end-1]
            
            reader = open(GFF3.Reader, file)
            output_df = DataFrame(
                contig_nr = Int64[],
                origin = String[],
                strand = GenomicFeatures.Strand[],
                seqend = Int64[],
                seqstart = Int64[]
            );
            
            contig_nr = 0
            for record in reader
                featuretype_o = GFF3.featuretype(record)
                if featuretype_o == "region"
                    contig_nr += 1
                end
                strand_o = GFF3.strand(record)
                seqend_o = GFF3.seqend(record)
                seqstart_o = GFF3.seqstart(record)
                push!(output_df, [contig_nr, featuretype_o, strand_o, seqend_o, seqstart_o])
            end
            close(reader)

            for g in groupby(output_df, :contig_nr)
                p_strand = @view g[[i == GenomicFeatures.STRAND_POS for i in g.strand] .&& g.origin .== "gene", :]
                n_strand = @view g[[i == GenomicFeatures.STRAND_NEG for i in g.strand] .&& g.origin .== "gene", :]
                
                p_nr_genes = size(p_strand, 1)
                n_nr_genes = size(n_strand, 1)
    
                genome_length = g.seqend[g.origin .== "region"][1]
            
                p_gene_lengths = length_interval(p_strand.seqstart, p_strand.seqend)
                n_gene_lengths = length_interval(n_strand.seqstart, n_strand.seqend)
            
                p_gene_length_sum = total_interval_length(p_gene_lengths)
                n_gene_length_sum = total_interval_length(n_gene_lengths)
            
                p_gap_lengths = space_between(p_strand.seqstart, p_strand.seqend)
                n_gap_lengths = space_between(n_strand.seqstart, n_strand.seqend)
            
                p_U_overlap_nr, p_U_overlap_length_sum = count_overlaps(p_gap_lengths)
                n_U_overlap_nr, n_U_overlap_length_sum = count_overlaps(n_gap_lengths)
            
                p_gap_length_sum = total_interval_length(p_gap_lengths)
                n_gap_length_sum = total_interval_length(n_gap_lengths)
            
                p_intervals = intervalify(p_strand.seqstart, p_strand.seqend)
                n_intervals = intervalify(n_strand.seqstart, n_strand.seqend)
            
                sort!(p_intervals)
                sort!(n_intervals)
            
                C_overlap_nr, D_overlap_nr, C_lengt_sum, D_lengt_sum = bidirectional_overlaps(n_intervals, p_intervals)
            
                push!(
                    storage_dataframe,
                    [genome_name, genome_length, p_nr_genes, p_gene_length_sum, n_nr_genes, n_gene_length_sum, p_U_overlap_nr, p_U_overlap_length_sum, n_U_overlap_nr, n_U_overlap_length_sum, p_gap_length_sum, n_gap_length_sum, C_overlap_nr, D_overlap_nr, C_lengt_sum, D_lengt_sum]
                )
            end
        catch e
            @error "Failed to process $file" exception=(e, catch_backtrace())
        end
    end

    # Final print to ensure a clean line
    println("\nProcessing complete.")

    return storage_dataframe
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

    file_paths = find_gff_files_recursive(input_dir)
    if isempty(file_paths)
        @warn "No GFF files found in $input_dir. Exiting."
        return
    end

    println("Searching for GFF files in $input_dir...")
    println("Found $(length(file_paths)) GFF files. Starting analysis...")

    df_metrics = get_genome_structure(file_paths)

    mkpath(dirname(output_path))
    
    println("Writing results to: $output_path")
    try
        CSV.write(output_path, df_metrics)
    catch e
        @error "Failed to write CSV file: $output_path. Please check file permissions or if the file is open in another application." exception=(e, catch_backtrace())
        exit(1)
    end

    println("Done.")
end

main()
