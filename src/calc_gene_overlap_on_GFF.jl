#-----------------------------------------------------------------
#   Description:    Determine overlaps in any gff file
#   Author:         SHP                    Date: 2022-2023
#   Notes:         
#-----------------------------------------------------------------

# load environment

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
    using CSV
    using DataFrames
    using GFF3
    using IntervalTrees
    using GenomicFeatures
    using BenchmarkTools
end

#-----------------------------------------------------------------
#   locating data 
#-----------------------------------------------------------------

output_folder = "./output/"
data_folder = ".../ncbi_downloads/bactera_reference/"

# R like list files functions
function list_files(dir::T) where T<:String
    err = isdir(dir)
    err == true || throw("$dir is not a directory")
	contents = String[]
	for (root, dirs, files) in walkdir(dir)
	    push!.(Ref(contents), joinpath.(root, files))
	end
	return contents
end

data_folder = "../ncbi_downloads/bactera_reference/"
data_folder = "../ncbi_downloads/archaea_reference/"
files = list_files(data_folder)

files = files[contains.(files, r"gff$")]

#-----------------------------------------------------------------
#   functions  
#-----------------------------------------------------------------

# length of interval
function length_interval(array_start::AbstractArray{T}, array_end::AbstractArray{T}) where T<:Int
    array_start .- array_end
end

# space between genes 
function space_between(array_start::AbstractArray{T}, array_end::AbstractArray{T}) where T<:Int
    array_start[2:length(array_start)] .- array_end[1:(length(array_end).-1)]
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
    for overlap ∈ eachoverlap(set1, set2)
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

function get_genome_structure(file_name::AbstractArray{T}) where T<:String

    storage_dataframe = DataFrame(
        # gene lengths
        genome_name = String[],
        contig_size = Int64[],
        p_gene_nr = Int64[],
        p_gene_length_sum = Int64[],
        n_gene_nr = Int64[],
        n_gene_length_sum = Int64[],
    
        # overlap number and total size
        p_U_overlap_nr = Int64[],
        p_U_overlap_length_sum = Int64[],
        n_U_overlap_nr = Int64[],
        n_U_overlap_length_sum = Int64[],
    
        # overlap number and total size
        p_gap_length_sum = Int64[],
        n_gap_length_sum = Int64[],
    
        C_overlap_nr = Int64[],
        D_overlap_nr = Int64[],
        C_lengt_sum = Int64[],
        D_lengt_sum = Int64[]
    )

    n = length(file_name)
    
    for (enum, file) in enumerate(file_name)

        if mod(enum, 500) == 0
            @info(string("file ", enum, " of ", n))
        end
    
        genome_name = file[55:69]
        output_df = DataFrame(
            contig_nr = Int64[],
            origin = String[],
            strand = GenomicFeatures.Strand[],
            seqend = Int64[],
            seqstart = Int64[]
        );
        
        reader = open(GFF3.Reader, file);
    
        # Iterate over records.
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
        end;
    
        close(reader);
    
        # redo this part taking care to consider nr of contigs
    
        for g in groupby(output_df, :contig_nr)
    
            p_strand = @view g[[i == STRAND_POS for i in g.strand] .&& g.origin .== "gene" .&& g.origin .!= "region", :];
            n_strand = @view g[[i == STRAND_NEG for i in g.strand] .&& g.origin .== "gene" .&& g.origin .!= "region", :];
            
            p_nr_genes = size(p_strand, 1)
            n_nr_genes = size(n_strand, 1)

            genome_length = g.seqend[g.origin .== "region"][1]
        
            # gene lengths
            p_gene_lengths = length_interval(p_strand.seqend, p_strand.seqstart)
            n_gene_lengths = length_interval(n_strand.seqend, n_strand.seqstart)
        
            p_gene_length_sum = total_interval_length(p_gene_lengths)
            n_gene_length_sum = total_interval_length(n_gene_lengths)
        
            # gap lenggths
            p_gap_lengths = space_between(p_strand.seqstart, p_strand.seqend)
            n_gap_lengths = space_between(n_strand.seqstart, n_strand.seqend)
        
            # overlap number and total size
            p_U_overlap_nr, p_U_overlap_length_sum = count_overlaps(p_gap_lengths)
            n_U_overlap_nr, n_U_overlap_length_sum = count_overlaps(n_gap_lengths)
        
            # overlap number and total size
            p_gap_length_sum = total_interval_length(p_gap_lengths)
            n_gap_length_sum = total_interval_length(n_gap_lengths)
        
            # comvergent and divergent overlaps
            # first generate overlap intervals
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
    
    end

    return storage_dataframe

end

#-----------------------------------------------------------------
#   performative action
#-----------------------------------------------------------------

struct_genome = get_genome_structure(files)

output_folder = "iodine/output/"

CSV.write(
    string(output_folder, "structure_reference_archaea_2025_02_04.csv"),
    struct_genome
)




