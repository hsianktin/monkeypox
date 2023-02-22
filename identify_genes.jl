# Based on reference genome and gene information
# and local sequence alignment 
# find the desired sequence.

# read fasta sequences, embedding them into a high dimensional space
using CSV, DataFrames, Pipe
using ProgressMeter
using BioAlignments # 3.0.0
using GenomicAnnotations # 0.3.2
using BioSequences # 3.1.0
using FASTX # 2.0.0
using Dates
using BioFetch # 0.2.1

# read meta data
# meta_df = CSV.read("metadata.csv", DataFrame)
meta_df = CSV.read("meta_data.csv", DataFrame)

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
end

# define reference genome by accession number
ref_acc = "NC_063383"
ref_genome = readgbk("data/$(ref_acc).gbk")[1]
features = [feature(entry) for entry in ref_genome.genes] |> unique
genes = [gene for gene in ref_genome.genes if feature(gene) == :gene]


@info "start alignment with respect to $(ref_acc)"
diffs = []
progress = Progress(length(records), 1)
for (record, acc) in zip(records, accs)
    diff = Array{Any}(undef, length(genes))
    genome_length = length(FASTX.sequence(LongDNA{4}, record))
    Threads.@threads for i ∈ eachindex(genes)
        gene = genes[i]
        gene_pos = locus(gene).position |> collect
        gene_start = gene_pos[1] |> (x -> minimum([50000 * floor(Int, x/50000) + 1, genome_length]))
        gene_end = gene_pos[end] |> (x -> minimum([50000 * ceil(Int, x/50000), genome_length]))
        gene_seq = GenomicAnnotations.sequence(gene)
        if locus(gene).strand == '+'
            genome_seq = @pipe record |> FASTX.sequence(LongDNA{4}, _) |> (x->x[gene_start:gene_end])
        else
            genome_seq = @pipe record |> FASTX.sequence(LongDNA{4}, _) |> (x->x[gene_start:gene_end]) |> reverse_complement
        end
        alignment = pairalign(
                        LocalAlignment(),
                        gene_seq,
                        genome_seq ,
                        AffineGapScoreModel(EDNAFULL, gap_open=-10,gap_extend=-2)
                    )
        diff[i] = alignment
    end
    for i in 1:length(diff)
        if !isassigned(diff,i)
            @warn "alignment over gene $i is not assigned"
        end
    end
    push!(diffs, diff)
    next!(progress; showvalues = [(:acc, acc)])
end

# bug fix for GenomicAnnotations.jl
function Base.delete!(genes::Vector{Gene})
    indices = index.(genes)
    DataFrames.delete!(parent(genes).genedata, Int.(indices))
    lastindices = length(parent(genes).genes) - length(indices) + 1 : length(parent(genes).genes)
    deleteat!(parent(genes).genes, lastindices)
    nothing
end
# construct the genomic annotations
@info "constructing the genomic annotations"
progress = Progress(length(records), 1)
bad_accs = []
for (acc, record, date, alignment_result) ∈ zip(accs, records, dates, diffs)
    # @show acc
    if isfile("data/$(acc).gbk")
        chr = readgbk("data/$(acc).gbk")[1]
    else 
        chr = fetchseq(acc, format=gb)[1]
    end
    genome_length = length(FASTX.sequence(LongDNA{4}, record))
    # end
    chr_annotations = [feature(entry) for entry in chr.genes] |> unique
    # delete all the previous :aligned_genes
    if :aligned_gene ∈ chr_annotations
        aligned_genes = [gene for gene in chr.genes if feature(gene) == :aligned_gene]
        delete!(aligned_genes)
    end
    
    # construct the annotations
    # @info "genetic annotations of $acc not found, constructing"
    bad_flag = false
    for j ∈ 1:length(genes)
        gene = genes[j]
        gene_pos = locus(gene).position |> collect
        gene_start = gene_pos[1] |> (x -> minimum([50000 * floor(Int, x/50000) + 1, genome_length]))
        gene_end = gene_pos[end] |> (x -> minimum([50000 * ceil(Int, x/50000), genome_length]))    
        strand = locus(gene).strand # whether the gene is on the positive strand
        alnⱼ = alignment_result[j]
        pairaln = alnⱼ.aln.a.aln
        if strand == '+'
            (start_nt, end_nt) = (pairaln.firstref,
                                    pairaln.lastref)
            (start_nt, end_nt) = (start_nt + gene_start - 1,
                                    end_nt + gene_start - 1)
        else
            (start_nt, end_nt) = (pairaln.firstref,
                                    pairaln.lastref)
            (start_nt, end_nt) = (gene_end - end_nt + 1,
                                    gene_end - start_nt + 1)
        end
        seq = [x[2] for x in alnⱼ.aln if x[2] ≠ DNA_Gap] |> LongDNA{4}
        # quality check 
        if strand == '+'
            focal_seq = @pipe record |> FASTX.sequence(LongDNA{4}, _) |>
                            (x->x[start_nt:end_nt])
        else 
            focal_seq = @pipe record |> FASTX.sequence(LongDNA{4}, _) |> 
                            (x->x[start_nt:end_nt]) |> reverse_complement
        end
        if seq ≠ focal_seq
            # @warn "alignment of gene $j is not correct" seq
            bad_flag = true
        end
        append_flag = (seq[1:3] == dna"ATG") # check if starts with start codon
        if append_flag
            newgene = addgene!(chr, :aligned_gene,
                    Locus(start_nt:end_nt,strand))
            newgene.gene = gene.gene # also stores the name of the gene.
            newgene.locus_tag = gene.locus_tag
        else
            # @warn "gene $j of $acc is not appended" seq
            bad_flag = true
        # end
        end
    end
    if bad_flag
        push!(bad_accs, acc)
    end
    # write to local file
    open(GenBank.Writer, "data/$(acc).gbk") do io
        write(io, chr)
    end
    next!(progress; showvalues = [(:acc, acc)])
end

# save a list of bad accessions
open("data/bad_accs.txt", "w") do io
    for acc in bad_accs
        println(io, acc)
    end
end