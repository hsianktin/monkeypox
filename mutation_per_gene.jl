# read fasta sequences, count the number of mutations per gene
using Distributed
addprocs([("rog-ld03",2), ("ryzen",2)])

@everywhere begin
    using CSV, DataFrames, Pipe
    using ProgressMeter
    using BioAlignments
    using GenomicAnnotations
    using BioSequences
    using FASTX
    using Dates
end
@everywhere using ProgressLogging
@everywhere using TerminalLoggers, ObservablePmap
@everywhere using Logging: global_logger
@everywhere global_logger(TerminalLogger())

# read meta data
# @everywhere CONTROL = [([DNA_A, DNA_C], [DNA_A, DNA_T]), ([DNA_G, DNA_T], [DNA_A, DNA_T])]
# @everywhere control_label = "AC2AT"
@everywhere CONTROL = [
    ([DNA_A], [DNA_C]), 
    ([DNA_T], [DNA_G]),
    # ([DNA_A], [DNA_T]), 
    # ([DNA_T], [DNA_A]),
    # ([DNA_A], [DNA_G]), 
    # ([DNA_T], [DNA_C]),
    ]
@everywhere control_label = "A2C"



ref_accs = [
    # "KP849470",
    "KJ642617",
    # "KJ642616",
    # "MK783029"
    # "NC_003310",
    # "JX878428",
    # "ON563414.3",
    # "ON676708",
    # "OP123040",
    # "OP257247"
]

@everywhere begin
    meta_df = CSV.read("meta_data.csv", DataFrame)
    # bad_df = CSV.read("data/bad_accs.txt", DataFrame)
    # meta_df = meta_df[[i for i ∈ eachindex(meta_df.Accession) if meta_df.Accession[i] ∉ bad_df.Accession], :]
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
    @everywhere  function QC(alg)
        return score(alg)/length(alignment(alg))
    end


    @everywhere function alg2snp(alg)
        alg = alg |> alignment |> collect
        seqₓ = [x[1] for x in alg if x[1] ≠ DNA_Gap]
        seqᵧ = [x[2] for x in alg if x[2] ≠ DNA_Gap]
        snp_frame = DataFrame(
            APOBEC = Int[],
            APOBEC⁻¹ = Int[], # reverse APOBEC
            CONTROL = Int[], # non-APOBEC relevant mutations AC -> AT, GT -> AT
            SYNONYMOUS = Int[],
        )
        if seqₓ[1:3] ≠ [DNA_A, DNA_T, DNA_G]
            # the query is not the reference
            # @warn "No start codon" seqₓ
            return snp_frame # return empty frame
        end
        # here x is the reference genome and y is the query

        posₓ = 1
        posᵧ = 1
        for i ∈ 2:length(alg)-1
            APOBEC = [([DNA_T,DNA_C], [DNA_T,DNA_T]), ([DNA_G, DNA_A], [DNA_A, DNA_A])]
            # CONTROL = [([DNA_A, DNA_C], [DNA_A, DNA_T]), ([DNA_G, DNA_T], [DNA_A, DNA_T])]
            xᵢ, yᵢ = alg[i]
            if xᵢ != DNA_Gap
                posₓ += 1
            end
            if yᵢ != DNA_Gap
                posᵧ += 1
            end
            xᵢ₋₁, yᵢ₋₁ = alg[i-1]
            xᵢ₊₁, yᵢ₊₁ = alg[i+1]
            # find the codon 
            ## initialize the codon
            codonₓ = AA_Term
            codonᵧ = AA_Term
            if posₓ % 3 == 0
                codonₓ = [seqₓ[posₓ-2], seqₓ[posₓ-1], seqₓ[posₓ]] |> LongDNA{4} |> BioSequences.translate
            elseif posₓ % 3 == 1 && posₓ + 2 ≤ length(seqₓ)
                codonₓ = [seqₓ[posₓ], seqₓ[posₓ+1], seqₓ[posₓ+2]] |> LongDNA{4} |> BioSequences.translate
            elseif posₓ % 3 == 2 && posₓ + 1 ≤ length(seqₓ)
                codonₓ = [seqₓ[posₓ-1], seqₓ[posₓ], seqₓ[posₓ+1]] |> LongDNA{4} |> BioSequences.translate
            end
            if posᵧ % 3 == 0
                codonᵧ = [seqᵧ[posᵧ-2], seqᵧ[posᵧ-1], seqᵧ[posᵧ]] |> LongDNA{4} |> BioSequences.translate
            elseif posᵧ % 3 == 1 && posᵧ + 2 ≤ length(seqᵧ)
                codonᵧ = [seqᵧ[posᵧ], seqᵧ[posᵧ+1], seqᵧ[posᵧ+2]] |> LongDNA{4} |> BioSequences.translate
            elseif posᵧ % 3 == 2 && posᵧ + 1 ≤ length(seqᵧ)
                codonᵧ = [seqᵧ[posᵧ-1], seqᵧ[posᵧ], seqᵧ[posᵧ+1]] |> LongDNA{4} |> BioSequences.translate
            end
            snp_flag = false
            APOBEC_flag = false
            SYNONYMOUS_flag = false
            CONTROL_flag = false
            APOBEC⁻¹_flag = false
            # @show x,y  # every pair of aligned characters, both x==y or x!=y
            if xᵢ ≠ yᵢ
                if xᵢ ∈ [DNA_A, DNA_T, DNA_G, DNA_C] && yᵢ ∈ [DNA_A, DNA_T, DNA_G, DNA_C]
                    snp_flag = true
                end
                if snp_flag
                    if codonₓ == codonᵧ
                        SYNONYMOUS_flag = true
                    end
                    if (
                        ([xᵢ₋₁, xᵢ], [yᵢ₋₁, yᵢ]) ∈ APOBEC ||
                        ([xᵢ, xᵢ₊₁], [yᵢ, yᵢ₊₁]) ∈ APOBEC
                    )
                        APOBEC_flag = true
                    elseif (
                            ([yᵢ₋₁, yᵢ], [xᵢ₋₁, xᵢ]) ∈ APOBEC ||
                            ([yᵢ, yᵢ₊₁], [xᵢ, xᵢ₊₁]) ∈ APOBEC
                    )
                        APOBEC⁻¹_flag = true
                    elseif (
                            # dinucleotide change
                            ([xᵢ₋₁, xᵢ], [yᵢ₋₁, yᵢ]) ∈ CONTROL ||
                            ([xᵢ, xᵢ₊₁], [yᵢ, yᵢ₊₁]) ∈ CONTROL ||
                            # single nucleotide change 
                            ([xᵢ], [yᵢ]) ∈ CONTROL
                    )
                        CONTROL_flag = true
                    end
                    # note: we will not count the same mutation twice by looking at (i-1) and (i+1) positions directly, because only if xᵢ ≠ yᵢ, we will count it
                    push!(snp_frame, [APOBEC_flag, APOBEC⁻¹_flag , CONTROL_flag ,SYNONYMOUS_flag])
                end
            end
        end
        return snp_frame
    end


  @everywhere  function pos2codon(nt, pos, seq)
        if pos % 3 == 0
            codon = [seq[pos-2], seq[pos-1], nt] |> LongDNA{4} |> BioSequences.translate
        elseif pos % 3 == 1 && pos + 2 ≤ length(seq)
            codon = [nt, seq[pos+1], seq[pos+2]] |> LongDNA{4} |> BioSequences.translate
        elseif pos % 3 == 2 && pos + 1 ≤ length(seq)
            codon = [seq[pos-1], nt, seq[pos+1]] |> LongDNA{4} |> BioSequences.translate
        else
            codon = AA_Term
        end
        return codon
    end

    @everywhere function alg2sym(alg)
        # identify all the possible synonymous mutations
        alg = alg |> alignment |> collect
        seqₓ = [x[1] for x in alg if x[1] ≠ DNA_Gap]
        seqᵧ = [x[2] for x in alg if x[2] ≠ DNA_Gap]
        if seqₓ[1:3] ≠ [DNA_A, DNA_T, DNA_G]
            # the query is not the reference
            # @warn "No start codon" seqₓ
            return (0, 0, 0, 0)
        end
        snp_frame = DataFrame(
            APOBEC = Int[],
            APOBEC⁻¹ = Int[], # reverse APOBEC
            CONTROL = Int[],
        )
        posₓ = 1
        posᵧ = 1
        for i ∈ 2:length(alg)-1
            APOBEC = [([DNA_T,DNA_C], [DNA_T,DNA_T]), ([DNA_G, DNA_A], [DNA_A, DNA_A])]
            
            xᵢ, yᵢ = alg[i]        
            xᵢ₋₁, yᵢ₋₁ = alg[i-1]
            xᵢ₊₁, yᵢ₊₁ = alg[i+1]
            if xᵢ != DNA_Gap
                posₓ += 1
                condonₓ = pos2codon(xᵢ, posₓ, seqₓ)
                for y ∈ [DNA_A, DNA_T, DNA_C, DNA_G]
                    if y ≠ xᵢ # (consider the xᵢ → y mutation)
                        condonᵧ = pos2codon(y, posₓ, seqₓ)
                        if condonₓ == condonᵧ
                            # synonymous
                            APOBEC_flag = false
                            APOBEC⁻¹_flag = false
                            CONTROL_flag = false
                            if (
                                ([xᵢ₋₁, xᵢ], [xᵢ₋₁, y]) ∈ APOBEC ||
                                ([y, xᵢ₊₁], [y, xᵢ₊₁]) ∈ APOBEC
                            )
                                APOBEC_flag = true
                            elseif (
                                ([xᵢ₋₁, y], [xᵢ₋₁, xᵢ]) ∈ APOBEC ||
                                ([xᵢ, xᵢ₊₁], [y, xᵢ₊₁]) ∈ APOBEC
                            )
                                APOBEC⁻¹_flag = true
                            elseif (
                                ([xᵢ₋₁, xᵢ], [xᵢ₋₁, y]) ∈ CONTROL ||
                                ([y, xᵢ₊₁], [y, xᵢ₊₁]) ∈ CONTROL ||
                                ([xᵢ], [y]) ∈ CONTROL
                            )
                                CONTROL_flag = true
                            end
                            push!(snp_frame, [APOBEC_flag, APOBEC⁻¹_flag, CONTROL_flag])
                        end
                    end
                end
            end
            if yᵢ != DNA_Gap
                posᵧ += 1
            end
            # # find the codon 
            # @show posᵧ
            # @show seqᵧ
        end
        return (nrow(snp_frame), sum(snp_frame.APOBEC), sum(snp_frame.APOBEC⁻¹), sum(snp_frame.CONTROL))
    end


    @everywhere function alg2sym_by_sites(alg)
        # identify all the possible synonymous mutations
        alg = alg |> alignment |> collect
        seqₓ = [x[1] for x in alg if x[1] ≠ DNA_Gap]
        seqᵧ = [x[2] for x in alg if x[2] ≠ DNA_Gap]
        if seqₓ[1:3] ≠ [DNA_A, DNA_T, DNA_G]
            # the query is not the reference
            # @warn "No start codon" seqₓ
            return (0, 0, 0)
        end
        snp_frame = DataFrame(
            APOBEC = Int[],
            APOBEC⁻¹ = Int[], # reverse APOBEC
            CONTROL = Int[],
        )
        posₓ = 1
        posᵧ = 1
        for i ∈ 2:length(alg)-1
            APOBEC = [([DNA_T,DNA_C], [DNA_T,DNA_T]), ([DNA_G, DNA_A], [DNA_A, DNA_A])]
            
            xᵢ, yᵢ = alg[i]        
            xᵢ₋₁, yᵢ₋₁ = alg[i-1]
            xᵢ₊₁, yᵢ₊₁ = alg[i+1]
            if xᵢ != DNA_Gap
                posₓ += 1
                condonₓ = pos2codon(xᵢ, posₓ, seqₓ)
                sym_flag = false
                APOBEC_flag = false
                APOBEC⁻¹_flag = false
                CONTROL_flag = false
                for y ∈ [DNA_A, DNA_T, DNA_C, DNA_G]
                    if y ≠ xᵢ # (consider the xᵢ → y mutation)
                        condonᵧ = pos2codon(y, posₓ, seqₓ)
                        if condonₓ == condonᵧ
                            sym_flag = true
                            # synonymous
                            if (
                                ([xᵢ₋₁, xᵢ], [xᵢ₋₁, y]) ∈ APOBEC ||
                                ([y, xᵢ₊₁], [y, xᵢ₊₁]) ∈ APOBEC
                            )
                                APOBEC_flag = true
                            elseif (
                                ([xᵢ₋₁, y], [xᵢ₋₁, xᵢ]) ∈ APOBEC ||
                                ([xᵢ, xᵢ₊₁], [y, xᵢ₊₁]) ∈ APOBEC
                            )
                                APOBEC⁻¹_flag = true
                            elseif (
                                ([xᵢ₋₁, xᵢ], [xᵢ₋₁, y]) ∈ CONTROL ||
                                ([y, xᵢ₊₁], [y, xᵢ₊₁]) ∈ CONTROL ||
                                ([xᵢ], [y]) ∈ CONTROL
                            )
                                CONTROL_flag = true
                            end
                        end
                    end
                end
                if sym_flag
                    push!(snp_frame, [APOBEC_flag, APOBEC⁻¹_flag, CONTROL_flag])
                end
            end
            if yᵢ != DNA_Gap
                posᵧ += 1
            end
            # # find the codon 
            # @show posᵧ
            # @show seqᵧ
        end
        return (nrow(snp_frame), sum(snp_frame.APOBEC), sum(snp_frame.APOBEC⁻¹), sum(snp_frame.CONTROL))
    end

    @everywhere function analysis(ref_acc)
        ref_genome = readgbk("data/$(ref_acc).gbk")[1]
        # sort according feature(gene) = {CDS, gene, source}
        genes = [gene for gene in ref_genome.genes if feature(gene) == :aligned_gene]
        if length(genes) < 160 # QC control
            @info "too many gene annotations lost"
        #     genes = [gene for gene in ref_genome.genes if feature(gene) == :gene]
        # end
        # if length(genes) < 160
        #     @info "Too few genes of $ref_acc annotated" length(genes)
        #     return nothing
        end
        @info "starting alignment"
        diffs = []
        @withprogress name="alignment with reference $(ref_acc)" for (record,acc,i) in zip(records,accs, eachindex(records))
            diff = Array{Any}(undef, length(genes))
            # genomeᵢ = readgbk("data/$(acc).gbk")[1]
            genome_length = length(FASTX.sequence(LongDNA{4}, record))
            Threads.@threads for i ∈ eachindex(genes)
                gene = genes[i]
                gene_pos = locus(gene).position |> collect
                gene_start = gene_pos[1] |> (x -> minimum([20000 * floor(Int, x/20000) + 1, genome_length]))
                gene_end = gene_pos[end] |> (x -> minimum([20000 * ceil(Int, x/20000), genome_length]))
                gene_seq = GenomicAnnotations.sequence(gene)  
                if locus(gene).strand == '+'
                    genome_seq = @pipe record |> FASTX.sequence(LongDNA{4}, _) |> 
                                    (x->x[gene_start:gene_end])
                else
                    genome_seq = @pipe record |> FASTX.sequence(LongDNA{4}, _) |> 
                                    (x->x[gene_start:gene_end]) |> reverse_complement
                end
                alignment = pairalign(
                    LocalAlignment(),
                    gene_seq,
                    genome_seq,
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
            @logprogress i/length(records)
        end
        #     save("tmp/local_alignment_$(ref_acc).jld", "diffs",diffs)
        # else
        #     @info "alignment records found, load alignment"
        #     diffs = load("tmp/local_alignment_$(ref_acc).jld", "diffs")
        # end
        @info "alignment finished"

        sym_mutations_per_genes = []
        nonsym_mutations_per_genes = []
        sym_mutation_sites_per_genes = []

        @info "counting SNPs"
        @withprogress name="SNP counting with reference $(ref_acc)" for i ∈ 1:length(records)
            print("$i out of $(length(records))      \r")
            diff = diffs[i]
            acc = accs[i]
            date = dates[i]
            sym_mutations_per_gene = Array{Any}(undef, length(genes), 4)
            sym_mutation_sites_per_gene = Array{Any}(undef, length(genes), 4)
            nonsym_mutations_per_gene = Array{Any}(undef, length(genes), 4)
            for j ∈ eachindex(diff) # each gene (alignment result) is represented by a diff[j]
                if QC(diff[j]) > 4
                    # filter out mismatches with low score
                    snp_frame = alg2snp(diff[j])
                    oₛᵧₘ  = sum(snp_frame.SYNONYMOUS)
                    oₐₚₒₛᵧₘ = sum(snp_frame.APOBEC .* snp_frame.SYNONYMOUS)
                    oₐₚₒ⁻¹ₛᵧₘ = sum(snp_frame.APOBEC⁻¹ .* snp_frame.SYNONYMOUS)
                    oᵪ = sum(snp_frame.CONTROL .* snp_frame.SYNONYMOUS)
                    sym_mutations_per_gene[j,:] = [oₛᵧₘ, oₐₚₒₛᵧₘ, oₐₚₒ⁻¹ₛᵧₘ, oᵪ]
                    noₐₚₒₛᵧₘ = sum(snp_frame.APOBEC .* (1 .- snp_frame.SYNONYMOUS))
                    noₐₚₒ⁻¹ₛᵧₘ = sum(snp_frame.APOBEC⁻¹ .* (1 .- snp_frame.SYNONYMOUS))
                    noᵪ = sum(snp_frame.CONTROL .* (1 .- snp_frame.SYNONYMOUS))
                    no = nrow(snp_frame) - oₛᵧₘ
                    nonsym_mutations_per_gene[j,:] = [no, noₐₚₒₛᵧₘ, noₐₚₒ⁻¹ₛᵧₘ, noᵪ]
                    nₛᵧₘ, nₐₚₒₛᵧₘ, nₐₚₒ⁻¹ₛᵧₘ, nᵪ = alg2sym(diff[j])
                    sym_mutation_sites_per_gene[j,:] = [nₛᵧₘ, nₐₚₒₛᵧₘ, nₐₚₒ⁻¹ₛᵧₘ, nᵪ]
                end
            end
            push!(sym_mutations_per_genes, sym_mutations_per_gene)
            push!(nonsym_mutations_per_genes, nonsym_mutations_per_gene)
            push!(sym_mutation_sites_per_genes, sym_mutation_sites_per_gene)
            # for each of gene, we have 12 numbers to describe the mutation
            # we can vertically stack them to get a matrix of length(genes) x 12
            # save the matrix in the form of a csv table. Each column is one of the 12 numbers
            # and each row is a gene
            record_df = DataFrame(
                gene_id = [], # index of the gene
                acc = [], # acc of the record
                date = [], # date of the record
                oₛᵧₘ = [], # synonymous mutations observed 
                oₐₚₒₛᵧₘ = [], # synonymous APOBEC mutations observed
                oₐₚₒ⁻¹ₛᵧₘ = [], # synonymous APOBEC⁻¹ mutations observed
                oᵪ = [], # synonymous control mutations observed
                noₐₚₒₛᵧₘ = [], # nonsynonymous APOBEC mutations observed
                noₐₚₒ⁻¹ₛᵧₘ = [], # nonsynonymous APOBEC⁻¹ mutations observed
                noᵪ = [], # nonsynonymous control mutations observed
                no = [], # nonsynonymous mutations observed
                nₛᵧₘ = [], # synonymous mutation sites
                nₐₚₒₛᵧₘ = [], # synonymous APOBEC mutation sites
                nₐₚₒ⁻¹ₛᵧₘ = [], # synonymous APOBEC⁻¹ mutation sites
                nᵪ = [], # synonymous control mutation sites
            )
            for j ∈ eachindex(diff)
                if QC(diff[j]) > 4
                    push!(record_df,
                    [
                        genes[j].gene,
                        acc,
                        date,
                        sym_mutations_per_gene[j,1],
                        sym_mutations_per_gene[j,2],
                        sym_mutations_per_gene[j,3],
                        sym_mutations_per_gene[j,4],
                        nonsym_mutations_per_gene[j,2],
                        nonsym_mutations_per_gene[j,3],
                        nonsym_mutations_per_gene[j,4],
                        nonsym_mutations_per_gene[j,1],
                        sym_mutation_sites_per_gene[j,1],
                        sym_mutation_sites_per_gene[j,2],
                        sym_mutation_sites_per_gene[j,3],
                        sym_mutation_sites_per_gene[j,4],
                    ])
                end
            end
            CSV.write("./data/mutations_per_gene_counts_$(control_label)_$(acc)_$(date).csv", record_df)
            @logprogress i/length(records)
        end
        @info "finished"
    end



# @showprogress pmap(analysis, accs)

summ, task = ologpmap( 1:length(ref_accs); logger_f=TerminalLogger) do id
    analysis(ref_accs[id])  
end

# Define a function that takes a string as an argument
function erase_output(str)
    # Count the number of "\n"s in the string
    n = count(isequal('\n'), str)
    # Print "\033[F" control sequences n times
    for i in 1:n
      print("\033[F")
    end
end

function replace_tab(str::AbstractString)
    # Use the replace function with a regular expression pattern to match \t and a replacement string of 4 spaces
    return replace(str, r"\t" => "    ")
end



# Define a function that takes a string as an argument
function pad_output(str)
    # Split the string by "\n" and store the lines in an array
    lines = split(str, "\n")
    # Find the maximum length of the lines
    max_len = maximum(length, lines)
    # Loop through the lines and insert spaces between "\n" to make them equal length
    for i in eachindex(lines)
      # Calculate the number of spaces needed
      n_spaces = max_len - length(lines[i])
      # Insert the spaces after the line
      lines[i] = lines[i] * " " ^ n_spaces
    end
    # Join the lines with "\n" and return the result
    return join(lines, "\n")
  end
# clc()
# erase_output(summ.val)

while !istaskdone(task)
    str = summ.val |> pad_output |> replace_tab
    print(str)
    sleep(1)
    erase_output(str)
end

println(summ.val |> pad_output)
@info "Done"

