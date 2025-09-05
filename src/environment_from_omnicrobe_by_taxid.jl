#!/usr/bin/env julia

#-----------------------------------------------------------------
#   Description:      Retrieve environment info from omnicrobe and output in a tidy format.
#   Author:           SHP
#   Date:             2025
#   Revised:          2025-09-05
#-----------------------------------------------------------------

#-----------------------------------------------------------------
# preamble
#-----------------------------------------------------------------
using Pkg
Pkg.activate(".") # Activate environment in current directory

using ArgParse
using CSV
using DataFrames
using HTTP
using JSON3

#-----------------------------------------------------------------
# Helper Function for API call
#-----------------------------------------------------------------
"""
    fetch_omnicrobe_env(taxid, base_url)

Asynchronously fetches environment data for a given taxid from the OmniMicrobe API.
Returns a tuple containing the original taxid and the fetched data.
"""
function fetch_omnicrobe_env(taxid::Int, base_url::String)
    full_url = string(base_url, taxid)
    try
        # Perform non-blocking GET request
        r = HTTP.get(full_url, status_exception=false) # Prevent exceptions for non-200 statuses

        if r.status == 200
            json_obj = JSON3.read(r.body)
            if isempty(json_obj)
                return (taxid, missing, missing)
            end
            
            # Extract and clean up the primary environment term from the list of forms
            envs = [item.obt_forms[1] for item in json_obj] # Take the first, most representative term
            obtids = [item.obtid for item in json_obj]
            return (taxid, envs, obtids)
        else
            @warn "Request failed for taxid $taxid with status $(r.status)"
            return (taxid, missing, missing)
        end
    catch e
        @error "Exception for taxid $taxid: $e"
        return (taxid, missing, missing)
    end
end

#-----------------------------------------------------------------
# Argument Parsing
#-----------------------------------------------------------------
function parse_commandline()
    s = ArgParseSettings(description="Fetch environment data from OmniMicrobe for a list of taxon IDs.")
    @add_arg_table! s begin
        "--input", "-i"
            help = "Input CSV file containing a 'taxId' column"
            required = true
        "--output-dir", "-o"
            help = "Path for the output directory. Two files will be created here."
            required = true
    end
    return parse_args(s)
end


#-----------------------------------------------------------------
# Main Logic
#-----------------------------------------------------------------

function main()
    args = parse_commandline()

    if !isfile(args["input"])
        @error "Input file not found: $(args["input"])"
        return
    end

    # Define output paths
    output_dir = args["output-dir"]
    mkpath(output_dir) # Ensure directory exists
    
    # Construct output filenames based on the input filename
    base_name = first(split(basename(args["input"]), '_'))
    genome_report_path = joinpath(output_dir, base_name * "_genomes_report.csv")
    environments_path = joinpath(output_dir, base_name * "_genome_environments.csv")

    println("Reading input file: $(args["input"])")
    df = CSV.File(args["input"]) |> DataFrame

    # API endpoint
    base_url = "https://omnicrobe.migale.inrae.fr/api/search/relations?taxid=ncbi%3A"

    unique_taxids = unique(skipmissing(df.taxId))
    n_taxids = length(unique_taxids)

    println("Querying OmniMicrobe for $n_taxids unique taxon IDs...")
    
    tasks = Vector{Task}()
    # Create an asynchronous task for each API call
    for taxid in unique_taxids
        task = @async fetch_omnicrobe_env(taxid, base_url)
        push!(tasks, task)
    end
    
    # Wait for all tasks to complete and collect the results
    results = []
    for (i, t) in enumerate(tasks)
        push!(results, fetch(t))
        # Update and print a simple progress bar
        percentage = round(Int, (i / n_taxids) * 100)
        progress_bar = "[" * repeat("=", percentage) * repeat(" ", 100 - percentage) * "]"
        print("\rProgress: $progress_bar $percentage% (TaxID $i of $n_taxids)")
    end
    println("\nProcessing complete.")

    # Create a mapping from taxId to accession for the environment table
    taxid_to_accession = Dict(zip(df.taxId, df.accession))

    # --- Create the Tidy Environments DataFrame ---
    env_df = DataFrame(accession=String[], environment=String[], obtId=String[])
    
    for (taxid, envs, obtids) in results
        if !ismissing(envs)
            accession = get(taxid_to_accession, taxid, "Unknown_Accession")
            for (env, obtid) in zip(envs, obtids)
                push!(env_df, (accession, env, obtid))
            end
        end
    end

    println("Writing tidy environment data to: $environments_path")
    CSV.write(environments_path, env_df, delim='\t')
    
    # --- Create and Write the Main Genomes Report ---
    # Remove the old environment columns if they exist
    if "environments" in names(df)
        select!(df, Not(:environments))
    end
    if "obtId" in names(df)
        select!(df, Not(:obtId))
    end
    
    println("Writing main genome report to: $genome_report_path")
    CSV.write(genome_report_path, df, delim='\t')

    println("Done.")
end

main()