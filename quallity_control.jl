# read fasta sequences, embedding them into a high dimensional space
# embedding based original annotations, instead of annotations based on alignment results to NC_063383
using Distributed
addprocs(3)

@everywhere begin
using CSV, DataFrames, Pipe
using ProgressMeter
using BioAlignments
using GenomicAnnotations
using BioSequences
using FASTX
using Dates
end
# read meta data
@everywhere begin
    meta_df = CSV.read("metadata.csv", DataFrame)

    records = []
    dates = Date[]
    accs = []
    # acquire sequences that are only older than 2021-01-01
    for (acc,date) ∈ zip(meta_df.Accession, meta_df.Date)
        if date > Date(1888,1,1)
            # read the fasta file
            record = FASTX.FASTA.Reader(open("data/$(acc).fasta")) |> collect
            # append the record to the list
            push!(records, record[1])
            # append the date to the list
            push!(dates, date)
            # append the accession to the list
            push!(accs, acc)
        end
        # record = FASTX.FASTA.Reader(open("data/$(acc).fasta")) |> collect
        # push!(records, record[1])
    end
    # records
end


ref_acc = "KJ642617"
ref_acc = "NC_063383"
ref_genome = readgbk("data/$(ref_acc).gbk")[1]
# sort according feature(gene) = {CDS, gene, source}
genes = [gene for gene in ref_genome.genes if feature(gene) == :gene]
aligned_genes = [gene for gene in ref_genome.genes if feature(gene) == :aligned_gene]
if length(genes) < 160 # QC control
    @info "too many gene annotations lost"
#     genes = [gene for gene in ref_genome.genes if feature(gene) == :gene]
# end
# if length(genes) < 160
#     @info "Too few genes of $ref_acc annotated" length(genes)
#     return nothing
end
# length(genes) == length(aligned_genes)
# @info "starting alignment"
diffs = []
record = records[10]
acc = accs[10]
diff = Array{Any}(undef, length(genes))
    # genomeᵢ = readgbk("data/$(acc).gbk")[1]
genome_length = length(FASTX.sequence(LongDNA{4}, record))
mismatch_counter = 0
for i ∈ 1:length(genes)
    gene = genes[i]
    gene_pos = locus(gene).position |> collect
    gene_start = gene_pos[1] |> (x -> minimum([50000 * floor(Int, x/50000) + 1, genome_length]))
    gene_end = gene_pos[end] |> (x -> minimum([50000 * ceil(Int, x/50000), genome_length]))
    gene_seq = GenomicAnnotations.sequence(gene)
    aligned_gene = aligned_genes[i]
    aligned_gene_pos = locus(aligned_gene).position |> collect
    aligned_gene_start = aligned_gene_pos[1] |> (x -> minimum([50000 * floor(Int, x/50000) + 1, genome_length]))
    aligned_gene_end = aligned_gene_pos[end] |> (x -> minimum([50000 * ceil(Int, x/50000), genome_length]))
    aligned_gene_seq = GenomicAnnotations.sequence(aligned_gene)
    if locus(gene).strand != locus(aligned_gene).strand
        @warn "strand not match" locus(gene).strand locus(aligned_gene).strand
        @info "gene locus" gene
        @info "aligned_gene locus" aligned_gene
        global mismatch_counter += 1
    end
end
@info "number of mismatched genes" mismatch_counter

