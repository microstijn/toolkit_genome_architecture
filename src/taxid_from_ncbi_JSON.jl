
#-----------------------------------------------------------------
#   Description:    From JSONL from NCBI retrieve INFO. 
#   Author:         SHP
#   Date:           2025
#-----------------------------------------------------------------
#-----------------------------------------------------------------

#-----------------------------------------------------------------
#= 
Notes. This data comes from downloading data from the NCBI
using the command line toolkit. 
I used the following command to download all refence gffs. 


datasets download genome taxon 2 --assembly-source refseq --reference --include gff3,gtf,seq-report --dehydrated --filename bacteria_reference.zip
datasets download genome taxon 2157 --assembly-source refseq --reference --include gff3,gtf,seq-report --dehydrated --filename archaea_reference.zip
datasets rehydrate --directory archaea_reference/ 
datasets rehydrate --directory bacteria_reference/ 

The command download all gff3, gts and seq reports. I am to query the accompagnying JSON file for info on the assmblies etc. 

the taxdump is obtained from ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
=#
#-----------------------------------------------------------------

#-----------------------------------------------------------------
# preamble
# activate env
#-----------------------------------------------------------------

begin
    using Pkg
    cd("d:")
    cd("programming/")
    Pkg.activate("iodine")
end

#-----------------------------------------------------------------
# load packages
#-----------------------------------------------------------------

using CSV # to write CSV.
using JSON3 # to read JSON.
using DataFrames # to be able to use dataframes in general. 

#-----------------------------------------------------------------
# Identify data to load
#-----------------------------------------------------------------

datB = "D:/ncbi_downloads/bactera_reference/ncbi_dataset/data/assembly_data_report.jsonl"
datA = "D:/ncbi_downloads/archaea_reference/ncbi_dataset/data/assembly_data_report.jsonl"
datAB = [datA, datB]

#-----------------------------------------------------------------
# Create function for full taxonomic identification based on NCBI taxdump
#-----------------------------------------------------------------

using Taxonomy

taxNodes = "../ncbi_downloads/taxdump/nodes.dmp"
taxNames = "../ncbi_downloads/taxdump/names.dmp"



function getFullTaxonomy(
    taxId,
    taxNodesLoc::String,
    taxNamesLoc::String,
    )

    errNode = isfile(taxNodesLoc)
    errNode == true || throw("$taxNodesLoc is not a file")
    errName = isfile(taxNamesLoc)
    errName == true || throw("$taxNamesLoc is not a file")

    db =  Taxonomy.DB(taxNodesLoc,taxNamesLoc)

    taxid_ordered_dataframe = DataFrame(
        superkingdom = String[],
        phylum = String[],
        class = String[],
        order = String[],
        family = String[],
        genus_tax = String[],
        species = String[]
    )

    for id in taxId

        superkingdom = "NA"
        phylum = "NA"
        class = "NA"
        order = "NA"
        family = "NA"
        genus_tax = "NA"
        species = "NA"
        lineage_tax = "NA"

        try 
            tax = Taxon(id, db)
            lineage_tax = Lineage(tax)
        catch;
            lineage_tax = "NA"
        end

        try 
            superkingdom = string(lineage_tax[:superkingdom])
        catch;
            superkingdom = "NA"
        end

        try 
            phylum = string(lineage_tax[:phylum])
        catch;
            phylum = "NA"
        end

        try 
            class = string(lineage_tax[:class])
        catch;
            class = "NA"
        end

        try 
            order = string(lineage_tax[:order])
        catch;
            order = "NA"
        end

        try 
            family = string(lineage_tax[:family])
        catch;
            family = "NA"
        end
    
        try 
            genus_tax = string(lineage_tax[:genus])
        catch;
            genus_tax = "NA"
        end

        try 
            species = string(lineage_tax[:species])
        catch;
            species = "NA"
        end

        push!(taxid_ordered_dataframe, [
            superkingdom,
            phylum,
            class,
            order,
            family,
            genus_tax,
            species
        ])
    end

    return taxid_ordered_dataframe

end

#-----------------------------------------------------------------
# Identify data to load
#-----------------------------------------------------------------

for (dat, n) in zip(datAB, ["Archaea", "Bacteria"])
    output_folder = "iodine/output/"
    s = JSON3.read(dat, jsonlines=true)

    fd = DataFrame(
        organismName        = Union{Missing, String}[],
        taxId               = Union{Missing, Int64}[],
        checkmSpeciesTaxId  = Union{Missing, Int}[],
        accession           = Union{Missing, String}[],
        completeness        = Union{Missing, Float64}[],
        gcPercent           = Union{Missing, Float64}[]
    
    )
    
    for acc in s
    
        taxId               = Union{Missing, Int64}[]
        organismName        = Union{Missing, String}[]
        checkmSpeciesTaxId  = Union{Missing, Int}[]
        completeness        = Union{Missing, Float64}[]
        accession           = Union{Missing, String}[]
        gcPercent           = Union{Missing, Float64}[]
    
        try
            push!(taxId, acc.organism.taxId)
        catch
            push!(taxId, missing)
        end
    
        try
            push!(organismName, acc.organism.organismName)
        catch
            push!(organismName, missing)
        end
    
        try
            push!(checkmSpeciesTaxId, acc.checkmInfo.checkmSpeciesTaxId)
        catch
            push!(checkmSpeciesTaxId, missing)
        end
    
        try
            push!(completeness, acc.checkmInfo.completeness)
        catch
            push!(completeness, missing)
        end
    
        try
            push!(accession, acc.accession)
        catch
            push!(accession, missing)
        end
    
        try
            push!(gcPercent, acc.assemblyStats.gcPercent)
        catch
            push!(gcPercent, missing)
        end
        
        push!(
            fd,
            [
                organismName[1],
                taxId[1],
                checkmSpeciesTaxId[1],
                accession[1],
                completeness[1],
                gcPercent[1]
            ]
        );
    
    end
    
    #-----------------------------------------------------------------
    # Write data to file
    #-----------------------------------------------------------------
    
    fullTax = getFullTaxonomy(
        fd.taxId,
        taxNodes,
        taxNames
    )
    
    fdTax = hcat(fullTax, fd)

    CSV.write(
        string(output_folder, n, "TaxId.csv"),
        fdTax
    )

end






