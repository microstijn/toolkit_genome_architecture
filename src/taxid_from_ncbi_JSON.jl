#-----------------------------------------------------------------
#   Description:    From JSONL from NCBI retrieve INFO.
#   Author:         SHP
#   Date:           2025
#   Revised:        2025-08-07
#-----------------------------------------------------------------
#= 
Notes:
This data comes from downloading data from the NCBI using the command line toolkit.
I used the following command to download all refence gffs:
datasets download genome taxon 2 --assembly-source refseq --reference --include gff3,gtf,seq-report --dehydrated --filename bacteria_reference.zip
datasets download genome taxon 2157 --assembly-source refseq --reference --include gff3,gtf,seq-report --dehydrated --filename archaea_reference.zip
datasets rehydrate --directory archaea_reference/ 
datasets rehydrate --directory bacteria_reference/ 

The command download all gff3, gts and seq reports.
I am to query the accompagnying JSON file for info on the assmblies etc. 

The taxdump is obtained from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
=#

#-----------------------------------------------------------------
# preamble: packages
#-----------------------------------------------------------------
using Pkg
Pkg.activate(".") # Activate environment in current directory

using ArgParse
using CSV
using DataFrames
using JSON3
using Taxonomy

#-----------------------------------------------------------------
# Helper Functions
#-----------------------------------------------------------------

"""
    get_nested(obj, keys...; default=missing)

Safely retrieve a nested value from a JSON object.
"""
function get_nested(obj, keys...; default=missing)
    val = obj
    for k in keys
        if !isnothing(val) && haskey(val, k)
            val = val[k]
        else
            return default
        end
    end
    return val
end

"""
    getFullTaxonomy(tax_ids, db)

For a list of taxon IDs, retrieve the full lineage for each.
"""
function getFullTaxonomy(tax_ids, db::Taxonomy.DB)
    ranks = [:superkingdom, :phylum, :class, :order, :family, :genus, :species]
    
    # We create a vector of Lineage objects directly, which guarantees
    # the same length as the input tax_ids.
    full_lineages = map(tax_id -> begin
        if ismissing(tax_id)
            return nothing
        end
        try
            return Lineage(Taxon(tax_id, db))
        catch
            return nothing
        end
    end, tax_ids)

    # Convert the lineages into a DataFrame with consistent columns
    tax_data = Dict{Symbol, Vector{Union{Missing, String}}}()
    for rank in ranks
        # The correct way to check for a rank in a Lineage object is using `in`.
        tax_data[rank] = [isnothing(lineage) ? missing : (rank in lineage ? string(lineage[rank]) : missing) for lineage in full_lineages]
    end

    return DataFrame(tax_data)
end


#-----------------------------------------------------------------
# Argument Parsing
#-----------------------------------------------------------------
function parse_commandline()
    s = ArgParseSettings(description="Extract assembly metadata from NCBI JSONL files and add full taxonomy.")
    @add_arg_table! s begin
        "--jsonl", "-j"
            help = "Input assembly_data_report.jsonl file(s)"
            nargs = '+'
            required = true
        "--nodes", "-n"
            help = "Path to NCBI taxdump nodes.dmp file"
            required = true
        "--names", "-m"
            help = "Path to NCBI taxdump names.dmp file"
            required = true
        "--output", "-o"
            help = "Output directory for the resulting CSV files"
            required = true
    end
    return parse_args(s)
end

#-----------------------------------------------------------------
# Main Logic
#-----------------------------------------------------------------
function main()
    args = parse_commandline()
    
    # Normalize paths before using them.
    nodes_dmp = normpath(args["nodes"])
    names_dmp = normpath(args["names"])

    println("Loading taxonomy database...")
    tax_db = Taxonomy.DB(nodes_dmp, names_dmp)
    
    for jsonl_path in args["jsonl"]
        jsonl_path = normpath(jsonl_path)
        if !isfile(jsonl_path)
            @error "File not found: $jsonl_path. Skipping."
            continue
        end

        println("Processing file: $jsonl_path")

        # Pre-allocate column vectors for performance
        organismName = Union{Missing, String}[]
        taxId = Union{Missing, Int}[]
        checkmSpeciesTaxId = Union{Missing, Int}[]
        accession = Union{Missing, String}[]
        completeness = Union{Missing, Float64}[]
        gcPercent = Union{Missing, Float64}[]
        
        # Correctly handle line-delimited JSON by reading the file line-by-line
        total_lines = length(readlines(jsonl_path)) # Count lines for progress bar
        open(jsonl_path, "r") do io
            for (i, line) in enumerate(eachline(io))
                record = JSON3.read(line)
                push!(organismName, get_nested(record, :organism, :organismName))
                push!(taxId, get_nested(record, :organism, :taxId))
                push!(accession, get_nested(record, :accession))
                push!(checkmSpeciesTaxId, get_nested(record, :checkmInfo, :checkmSpeciesTaxId))
                push!(completeness, get_nested(record, :checkmInfo, :completeness))
                push!(gcPercent, get_nested(record, :assemblyStats, :gcPercent))

                # Update and print a simple progress bar
                percentage = round(Int, (i / total_lines) * 100)
                progress_bar = "[" * repeat("=", percentage) * repeat(" ", 100 - percentage) * "]"
                print("\rProgress: $progress_bar $percentage% (Record $i of $total_lines)")
            end
        end
        println("\nProcessing complete.")

        # Create DataFrame from collected data
        df = DataFrame(
            organismName = organismName,
            taxId = taxId,
            checkmSpeciesTaxId = checkmSpeciesTaxId,
            accession = accession,
            completeness = completeness,
            gcPercent = gcPercent
        )

        println("Retrieving full taxonomy for $(nrow(df)) entries...")
        full_tax = getFullTaxonomy(df.taxId, tax_db)
        
        # Combine the dataframes
        df_final = hcat(df, full_tax)

        # Write to output file
        output_filename = first(split(basename(jsonl_path), '.')) * "_TaxId.csv"
        output_path = normpath(joinpath(args["output"], output_filename))
        
        println("Writing results to: $output_path")
        try
            CSV.write(output_path, df_final)
        catch e
            @error "Failed to write CSV file: $output_path" exception=(e, catch_backtrace())
            exit(1)
        end
    end
    println("Done.")
end

main()