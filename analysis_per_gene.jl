using CSV, DataFrames, Pipe
using ProgressMeter
using BioAlignments
using GenomicAnnotations
using BioSequences
using FASTX
using Dates

control_label="A2C"
ref_acc = "KJ642617"
meta_df = CSV.read("meta_data.csv", DataFrame)

post2016_df = filter(row -> row[:Date] >= Date("2016-01-1"), meta_df)
acc = post2016_df[1,:Accession]
date = post2016_df[1,:Date]
df = CSV.read("./data/mutations_per_gene_counts_$(control_label)_$(acc)_$(date).csv", DataFrame)
gene_ids = unique(df.gene_id) |> sort

function get_mutation_vec(acc,date)
    df = CSV.read("./data/mutations_per_gene_counts_$(control_label)_$(acc)_$(date).csv", DataFrame)
    # vₐₛ = # of APOBEC-induced synonymous mutations
    # vₐₙ = # of APOBEC-induced non-synonymous mutations
    vₐₛ = zeros(Int, length(gene_ids))
    vₐₙ = zeros(Int, length(gene_ids))
    for (i,gene_id) in enumerate(gene_ids)
        vₐₛ[i] = sum(df[df.gene_id .== gene_id, :oₐₚₒₛᵧₘ])
        vₐₙ[i] = sum(df[df.gene_id .== gene_id, :noₐₚₒₛᵧₘ])
    end
    return vₐₛ, vₐₙ
end

V = get_mutation_vec.(post2016_df.Accession, post2016_df.Date)
vₐₛ = [v[1] for v in V]
vₐₙ = [v[2] for v in V]

# reshape to create #=length(gene_ids) random variables
# each random variable is represented by a vector of length #=length(post2016_df.Accession)
xₐₛ = [
    [vₐₛ[n][id] for n in eachindex(vₐₛ)] 
        for id in eachindex(gene_ids)]
xₐₙ = [
    [vₐₙ[n][id] for n in eachindex(vₐₙ)] 
        for id in eachindex(gene_ids)]

# calculate the correlation corr(xₐₛ[i], xₐₙ[j]) for each i,j 
# and store the result in a matrix
corr_mat = zeros(length(gene_ids), length(gene_ids))
using Statistics
for i in eachindex(gene_ids)
    for j in eachindex(gene_ids)
            corr_mat[i,j] = cor(xₐₛ[i], xₐₛ[j])
    end
end
# replace NaN with 0
corr_mat = replace(corr_mat, NaN=>0)
@show corr_mat
# save the correlation matrix
CSV.write("corr_mat_syn.csv", DataFrame(corr_mat,:auto), headers=false)

# convert to a pgfplots compatible table where the first two columns records i and j and the third column records the correlation
pgfplots_table = DataFrame(i=Int[], j=Int[], corr=Float64[])
for i in eachindex(gene_ids)
    for j in eachindex(gene_ids)
        push!(pgfplots_table, (i=i, j=j, corr=corr_mat[i,j]))
    end
end
CSV.write("corr_mat_syn_pgfplots.csv", pgfplots_table, headers=false)

corr_mat = zeros(length(gene_ids), length(gene_ids))
using Statistics
for i in eachindex(gene_ids)
    for j in eachindex(gene_ids)
            corr_mat[i,j] = cor(xₐₛ[i], xₐₙ[j])
    end
end
# replace NaN with 0
corr_mat = replace(corr_mat, NaN=>0)
@show corr_mat
# save the correlation matrix
CSV.write("corr_mat_syn_non.csv", DataFrame(corr_mat,:auto), headers=false)

# convert to a pgfplots compatible table where the first two columns records i and j and the third column records the correlation
pgfplots_table = DataFrame(i=Int[], j=Int[], corr=Float64[])
for i in eachindex(gene_ids)
    for j in eachindex(gene_ids)
        push!(pgfplots_table, (i=i, j=j, corr=corr_mat[i,j]))
    end
end
CSV.write("corr_mat_syn_non_pgfplots.csv", pgfplots_table, headers=false)