#-----------------------------------------------------------------
#   Description:    Retrieve environment info from omnicrobe.
#   Author:         SHP
#   Date:           2025
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
async function fetch_omnicrobe_env(taxid::Int, base_url::String)
    full_url = string(base_url, taxid)
    try
        # Perform non-blocking GET request
        r = HTTP.get(full_url, status_exception=false) # Prevent exceptions for non-200 statuses

        if r.status == 200
            json_obj = JSON3.read(r.body)
            if isempty(json_obj)
                return (taxid, missing, missing)
            end
            
            envs = [item.obt_forms for item in json_obj]
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
        "--output", "-o"
            help = "Path for the output CSV file"
            required = true
    end
    return parse_args(s)
end


#-----------------------------------------------------------------
# Main Logic
#-----------------------------------------------------------------

function taxid_to_environment()
    args = parse_commandline()

    if !isfile(args["input"])
        @error "Input file not found: $(args["input"])"
        return
    end

    println("Reading input file: $(args["input"])")
    df = CSV.File(args["input"]) |> DataFrame

    # API endpoint
    base_url = "https://omnicrobe.migale.inrae.fr/api/search/relations?taxid=ncbi%3A"

    tasks = []
    unique_taxids = unique(skipmissing(df.taxId))
    
    println("Querying OmniMicrobe for $(length(unique_taxids)) unique taxon IDs...")
    
    # Create an asynchronous task for each API call
    for taxid in unique_taxids
        task = @async fetch_omnicrobe_env(taxid, base_url)
        push!(tasks, task)
    end

    # Wait for all tasks to complete and collect the results
    results = [fetch(t) for t in tasks]
    
    # Convert the results into a DataFrame for easy merging
    results_df = DataFrame(
        taxId = [r[1] for r in results],
        environments = [r[2] for r in results],
        obtId = [r[3] for r in results]
    )

    println("Merging results...")
    # Join the fetched data back to the original dataframe
    final_df = leftjoin(df, results_df, on=:taxId)
    
    println("Writing results to: $(args["output"])")
    CSV.write(args["output"], final_df, delim='\t')
    println("Done.")
end

taxid_to_environment()
