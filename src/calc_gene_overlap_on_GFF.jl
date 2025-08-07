#-----------------------------------------------------------------
#   Description:    Determine overlaps in any GFF file.
#   Author:         SHP
#   Date:           2022-2023
#   Revised:        2025-08-07
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

#-----------------------------------------------------------------
# Functions (Original logic retained as it is sound)
#-----------------------------------------------------------------

# length of interval
length_interval(array_start, array_end) = array_end .- array_start

# space between genes
space_between(array_start, array_end) = array_start[2:end] .- array_end[1:end-1]

# count unidirectional overlaps
function count_overlaps(gap_array)
    overlaps = @view gap_array[gap_array .<= 0]
    return length(overlaps), isempty(overlaps) ? 0 : sum(overlaps)
end

# sum of interval lengths
total_interval_length(length_array) = isempty(length_array) ? 0 : sum(length_array)

# create Interval collection
intervalify(s, e) = IntervalCollection([Interval(a, b) for (a, b) in zip(s, e)])

# calculate convergent/divergent overlaps
function bidirectional_overlaps(set1, set2)
    convergent_overlaps = Int[]
    divergent_overlaps = Int[]
    
    for (interval1, interval2) in eachoverlap(set1, set2)
        # Assuming set1 is negative strand and set2 is positive strand
        # Convergent: --> <-- (end of positive overlaps start of negative)
        if strand(interval1) == STRAND_NEG && strand(interval2) == STRAND_POS
             push!(convergent_overlaps, min(last(interval1), last(interval2)) - max(first(interval1), first(interval2)))
        # Divergent: <-- --> (end of negative overlaps start of positive)
        elseif strand(interval1) == STRAND_POS && strand(interval2) == STRAND_NEG
             push!(divergent_overlaps, min(last(interval1), last(interval2)) - max(first(interval1), first(interval2)))
        end
    end
    
    len_con = length(convergent_overlaps)
    len_di = length(divergent_overlaps)
    sum_con = total_interval_length(convergent_overlaps)
    sum_di = total_interval_length(divergent_overlaps)
    
    return len_con, len_di, sum_con, sum_di
end


"""
    list_gff_files(dir::String)
Recursively finds all GFF files in a directory.
"""
function list_gff_files(dir::String)
    gff_files = String[]
    for (root, _, files) in walkdir(dir)
        for file in files
            if endswith(lowercase(file), ".gff") || endswith(lowercase(file), ".gff3")
                push!(gff_files, joinpath(root, file))
            end
        end
    end
    return gff_files
end

#-----------------------------------------------------------------
# Main Analysis Function
#-----------------------------------------------------------------
function get_genome_structure(file_paths::AbstractArray{String})
    storage_dataframe = DataFrame(
        genome_name = String[],
        contig_id = String[],
        contig_size = Int[],
        p_gene_nr = Int[], p_gene_length_sum = Int[],
        n_gene_nr = Int[], n_gene_length_sum = Int[],
        p_U_overlap_nr = Int[], p_U_overlap_length_sum = Int[],
        n_U_overlap_nr = Int[], n_U_overlap_length_sum = Int[],
        p_gap_length_sum = Int[], n_gap_length_sum = Int[],
        C_overlap_nr = Int[], D_overlap_nr = Int[],
        C_length_sum = Int[], D_length_sum = Int[]
    )
    
    n_files = length(file_paths)
    for (enum, file_path) in enumerate(file_paths)
        if mod(enum, 100) == 0
            @info("Processing file $enum of $n_files: $file_path")
        end
        
        # Robustly get genome name from filename
        genome_name = first(split(basename(file_path), r"(\.gff|\.gff3)"))
        
        reader = open(GFF3.Reader, file_path)
        # Group records by contig/sequence ID
        for contig_records in GFF3.grouprecords(reader, "FASTA_ID")
            
            # Extract region/contig info
            region = filter(r -> GFF3.featuretype(r) == "region", contig_records)
            if isempty(region)
                @warn "No 'region' feature found for a contig in $file_path. Skipping contig."
                continue
            end
            contig_id = GFF3.seqid(first(region))
            contig_length = GFF3.seqend(first(region))

            # Filter for genes on positive and negative strands
            genes = filter(r -> GFF3.featuretype(r) == "gene", contig_records)
            p_strand = sort(filter(r -> GFF3.strand(r) == STRAND_POS, genes), by=GFF3.seqstart)
            n_strand = sort(filter(r -> GFF3.strand(r) == STRAND_NEG, genes), by=GFF3.seqstart)

            # Gene counts
            p_nr_genes = length(p_strand)
            n_nr_genes = length(n_strand)

            # Gene lengths
            p_gene_lengths = length_interval.([GFF3.seqstart(r) for r in p_strand], [GFF3.seqend(r) for r in p_strand])
            n_gene_lengths = length_interval.([GFF3.seqstart(r) for r in n_strand], [GFF3.seqend(r) for r in n_strand])
            p_gene_length_sum = total_interval_length(p_gene_lengths)
            n_gene_length_sum = total_interval_length(n_gene_lengths)

            # Unidirectional Overlaps & Gaps
            p_gaps = space_between([GFF3.seqstart(r) for r in p_strand], [GFF3.seqend(r) for r in p_strand])
            n_gaps = space_between([GFF3.seqstart(r) for r in n_strand], [GFF3.seqend(r) for r in n_strand])
            p_U_overlap_nr, p_U_overlap_length_sum = count_overlaps(p_gaps)
            n_U_overlap_nr, n_U_overlap_length_sum = count_overlaps(n_gaps)
            p_gap_length_sum = sum(filter(x -> x > 0, p_gaps))
            n_gap_length_sum = sum(filter(x -> x > 0, n_gaps))

            # Bidirectional Overlaps
            p_intervals = intervalify([GFF3.seqstart(r) for r in p_strand], [GFF3.seqend(r) for r in p_strand])
            n_intervals = intervalify([GFF3.seqstart(r) for r in n_strand], [GFF3.seqend(r) for r in n_strand])
            C_overlap_nr, D_overlap_nr, C_length_sum, D_length_sum = bidirectional_overlaps(p_intervals, n_intervals)

            push!(storage_dataframe, 
                (genome_name, contig_id, contig_length, p_nr_genes, p_gene_length_sum, n_nr_genes, n_gene_length_sum, 
                p_U_overlap_nr, p_U_overlap_length_sum, n_U_overlap_nr, n_U_overlap_length_sum, 
                p_gap_length_sum, n_gap_length_sum, C_overlap_nr, D_overlap_nr, C_length_sum, D_length_sum)
            )
        end
        close(reader)
    end
    return storage_dataframe
end

#-----------------------------------------------------------------
# Argument Parsing and Execution
#-----------------------------------------------------------------
function parse_commandline()
    s = ArgParseSettings(description="Calculate gene overlap and genome structure metrics from GFF files.")
    @add_arg_table! s begin
        "--inputdir", "-d"
            help = "Input directory containing GFF files"
            required = true
        "--output", "-o"
            help = "Path for the output CSV file"
            required = true
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()

    if !isdir(args["inputdir"])
        @error "Input directory not found: $(args["inputdir"])"
        return
    end

    println("Searching for GFF files in $(args["inputdir"])...")
    gff_files = list_gff_files(args["inputdir"])

    if isempty(gff_files)
        @warn "No GFF files found. Exiting."
        return
    end
    
    println("Found $(length(gff_files)) GFF files. Starting analysis...")
    struct_genome = get_genome_structure(gff_files)
    
    println("Writing results to: $(args["output"])")
    CSV.write(args["output"], struct_genome)
    println("Done.")
end

main()
