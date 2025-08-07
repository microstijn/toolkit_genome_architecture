#-----------------------------------------------------------------
#   Description:    retrieve info from omnicrobe
#   Author:         SHP
#   Date:           2025    
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
#   preamble : packages    
#-----------------------------------------------------------------

begin
    using JSON3
    using DataFrames
    using CSV
    using HTTP
end

#-----------------------------------------------------------------
#   do the deed
#-----------------------------------------------------------------

# load file with taxonomic ids
files = readdir("iodine/output/",join = true)
files = files[occursin.(r".*TaxId\.csv", files)]

file = CSV.File(files[2]) |> DataFrame

# preacclocate columns
file.environments = missings(Vector{Array{String}}, size(file, 1))
file.obtId = missings(Vector{String}, size(file, 1))

sort!(file, :superkingdom)

url = "https://omnicrobe.migale.inrae.fr/api/search/relations?taxid=ncbi%3A"

z = 0
for (i, taxid) in enumerate(file.taxId)

    # skip missing and Viruses
    if  ismissing(file.superkingdom[i])
        continue
    end 

    # loop counter
    if mod(i, 500) == 0
        @info(string("file ", i, " of ", size(file, 1)))
    end

    # query omnicrobe
    r = HTTP.request("GET", string(url, taxid))
    b_obj = JSON3.read(r.body)

    envs = Array{String}[]
    obtids = String[]

    for depth_1 in b_obj
        env = [i for i in depth_1.obt_forms]
        obtid = depth_1.obtid
        push!(envs, env)
        push!(obtids, obtid)
    end

    if isempty(envs)
        file.environments[i] = missing
        file.obtId[i] = missing
    else
        file.environments[i] = envs
        file.obtId[i] = obtids
    end

end



CSV.write("iodine/output/BacteriaTaxIdOmni.csv", file, delim = "\t")

