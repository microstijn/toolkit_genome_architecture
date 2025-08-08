#!/usr/bin/env julia

#=====================================================
# Description:   Controls anc combines the separate scripts.
# Author:        SHP
# Date:          2025
# Revised:       2025-08-07
=====================================================#
using Pkg
Pkg.activate(".")
using ArgParse

# --- Script Definitions ---
const SCRIPT_DIR = @__DIR__
const TAXONOMY_SCRIPT = joinpath(SCRIPT_DIR, "taxid_from_ncbi_JSON.jl")
const ENVIRONMENT_SCRIPT = joinpath(SCRIPT_DIR, "environment_from_omnicrobe_by_taxid.jl")
const ARCHITECTURE_SCRIPT = joinpath(SCRIPT_DIR, "calc_gene_overlap_on_GFF.jl")
"""
    run_command(command, step_name)

Executes a command, prints status messages, and handles errors.
"""
function run_command(command::Cmd, step_name::String)
    println("--- Starting Step: $(step_name) ---")
    println("Executing: ", command)
    try
        # The `run` command will throw an error on failure, halting the script.
        run(command)
        println("--- Completed Step: $(step_name) ---\n")
    catch e
        @error "ERROR: Step '$step_name' failed." exception=(e, catch_backtrace())
        exit(1)
    end
end

"""
    parse_commandline()

Parses command-line arguments for the pipeline.
"""
function parse_commandline()
    s = ArgParseSettings(
        description = "Master Pipeline Controller for Genome Architecture Analysis.",
        autofix_names = true
    )
    @add_arg_table! s begin
        "--jsonl-files"
            help = "One or more paths to the NCBI 'assembly_data_report.jsonl' files."
            nargs = '+'
            required = true
        "--gff-dir"
            help = "Path to the input directory containing all GFF files."
            required = true
        "--taxdump-dir"
            help = "Path to the directory containing 'nodes.dmp' and 'names.dmp'."
            required = true
        "--output-dir"
            help = "Path to the directory where all output files will be saved."
            required = true
        "--julia-executable"
            help = "Path to the Julia executable (if not in system PATH)."
            default = "julia"
    end
    return parse_args(s)
end

"""
    main()

Main function to orchestrate the pipeline execution.
"""
function main()
    args = parse_commandline()

    # --- Setup ---
    # Create the main output directory if it doesn't exist.
    mkpath(args["output_dir"])

    # Define paths for taxdump files, normalizing them to prevent path errors.
    nodes_dmp = normpath(joinpath(args["taxdump_dir"], "nodes.dmp"))
    names_dmp = normpath(joinpath(args["taxdump_dir"], "names.dmp"))

    # --- Step 1: Extract Taxonomy from NCBI JSONL ---
    # Construct the command as an array of strings to avoid errors with
    # interpolating argument lists. The `...` is used to "splat" the
    # elements of the jsonl_files array into the command.
    taxonomy_cmd = Cmd([
        args["julia_executable"],
        TAXONOMY_SCRIPT,
        "--jsonl", args["jsonl_files"]...,
        "--nodes", nodes_dmp,
        "--names", names_dmp,
        "--output", normpath(args["output_dir"])
    ])
    run_command(taxonomy_cmd, "Extract Taxonomy and Assembly Info")

    # --- Step 2: Fetch Environment Data from OmniMicrobe ---
    step1_outputs = filter(f -> occursin("TaxId.csv", f), readdir(args["output_dir"]))
    for tax_file in step1_outputs
        # Normalize paths inside the loop as well.
        input_path = normpath(joinpath(args["output_dir"], tax_file))
        output_path = normpath(joinpath(args["output_dir"], replace(tax_file, "TaxId" => "TaxId_Omni")))
        
        # Construct the command as an array for robustness.
        environment_cmd = Cmd([
            args["julia_executable"],
            ENVIRONMENT_SCRIPT,
            "--input", input_path,
            "--output", output_path
        ])
        run_command(environment_cmd, "Fetch Environment Data for $(tax_file)")
    end

    # --- Step 3: Calculate Genome Architecture from GFF files ---
    # Normalize paths for this step as well.
    architecture_output = normpath(joinpath(args["output_dir"], "genome_architecture_metrics.csv"))
    architecture_cmd = Cmd([
        args["julia_executable"],
        ARCHITECTURE_SCRIPT,
        "--inputdir", normpath(args["gff_dir"]),
        "--output", architecture_output
    ])
    run_command(architecture_cmd, "Calculate Genome Architecture")

    println("Pipeline finished successfully!")
end

# --- Run the pipeline ---
main()
