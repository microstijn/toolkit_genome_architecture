#!/usr/bin/env julia

#=====================================================
# Description:   Controls and combines the separate scripts.
# Author:        SHP
# Date:          2025
# Revised:       2025-09-05
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

Parses command-line arguments for the main pipeline controller.
"""
function parse_commandline()
    s = ArgParseSettings(description="Main pipeline controller for genome architecture and taxonomy analysis.")
    @add_arg_table! s begin
        "--jsonl-files", "-j"
            help = "Paths to one or more NCBI assembly_data_report.jsonl files"
            nargs = '+'
            required = true
        "--gff-dir", "-g"
            help = "Path to the directory containing all GFF genome files"
            required = true
        "--taxdump-dir", "-t"
            help = "Path to the directory containing nodes.dmp and names.dmp"
            required = true
        "--output-dir", "-o"
            help = "Directory where all results will be saved"
            required = true
        "--julia-executable"
            help = "Path to the Julia executable"
            default = Sys.iswindows() ? "julia.exe" : "julia"
    end
    return parse_args(s)
end

function main()
    args = parse_commandline()
    mkpath(args["output-dir"])

    # --- Step 1: Extract Taxonomy and Assembly Info ---
    nodes_dmp = normpath(joinpath(args["taxdump-dir"], "nodes.dmp"))
    names_dmp = normpath(joinpath(args["taxdump-dir"], "names.dmp"))

    # Construct the command as an array for robustness against spaces in paths.
    # CORRECTED KEY: Use "julia-executable" with a hyphen
    taxonomy_cmd = Cmd([
        args["julia-executable"],
        TAXONOMY_SCRIPT,
        "--jsonl", args["jsonl-files"]...,
        "--nodes", nodes_dmp,
        "--names", names_dmp,
        "--output", normpath(args["output-dir"])
    ])
    run_command(taxonomy_cmd, "Extract Taxonomy and Assembly Info")

    # --- Step 2: Fetch Environment Data from OmniMicrobe ---
    step1_outputs = filter(f -> occursin("_TaxId.csv", f), readdir(args["output-dir"]))
    for tax_file in step1_outputs
        input_path = normpath(joinpath(args["output-dir"], tax_file))
        
        # CORRECTED KEY: Use "julia-executable" with a hyphen
        environment_cmd = Cmd([
            args["julia-executable"],
            ENVIRONMENT_SCRIPT,
            "--input", input_path,
            "--output-dir", normpath(args["output-dir"])
        ])
        run_command(environment_cmd, "Fetch Environment Data for $(tax_file)")
    end

    # --- Step 3: Calculate Genome Architecture from GFF files ---
    architecture_output = normpath(joinpath(args["output-dir"], "genome_architecture_metrics.csv"))
    # CORRECTED KEY: Use "julia-executable" with a hyphen
    architecture_cmd = Cmd([
        args["julia-executable"],
        ARCHITECTURE_SCRIPT,
        "--inputdir", normpath(args["gff-dir"]),
        "--output", architecture_output
    ])
    run_command(architecture_cmd, "Calculate Genome Architecture")

    println("Pipeline finished successfully!")
end

main()