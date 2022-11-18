# read fasta sequences, embedding them into a high dimensional space
using Distributed
addprocs(2)

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
                    end
                    push!(snp_frame, [APOBEC_flag, APOBEC⁻¹_flag ,SYNONYMOUS_flag])
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
            return (0, 0, 0)
        end
        snp_frame = DataFrame(
            APOBEC = Int[],
            APOBEC⁻¹ = Int[], # reverse APOBEC
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
                            end
                            push!(snp_frame, [APOBEC_flag, APOBEC⁻¹_flag])
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
        return (length(snp_frame.APOBEC), sum(snp_frame.APOBEC), sum(snp_frame.APOBEC⁻¹))
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
                            end
                        end
                    end
                end
                if sym_flag
                    push!(snp_frame, [APOBEC_flag, APOBEC⁻¹_flag])
                end
            end
            if yᵢ != DNA_Gap
                posᵧ += 1
            end
            # # find the codon 
            # @show posᵧ
            # @show seqᵧ
        end
        return (length(snp_frame.APOBEC), sum(snp_frame.APOBEC), sum(snp_frame.APOBEC⁻¹))
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
        @showprogress 1 for (record,acc) in zip(records,accs)
            diff = Array{Any}(undef, length(genes))
            # genomeᵢ = readgbk("data/$(acc).gbk")[1]
            genome_length = length(FASTX.sequence(LongDNA{4}, record))
            Threads.@threads for i ∈ 1:length(genes)
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
        end
        #     save("tmp/local_alignment_$(ref_acc).jld", "diffs",diffs)
        # else
        #     @info "alignment records found, load alignment"
        #     diffs = load("tmp/local_alignment_$(ref_acc).jld", "diffs")
        # end


        info_df = DataFrame(
            acc = String[],
            ref_acc = String[],
            date = Date[],
            year = Int[],
            Nₛᵧₘ = Float32[],
            Nₐₚₒₛᵧₘ = Int[],
            Nₐₚₒₛᵧₘ⁻¹ = Int[],
            Nₐ₊₋ = Int[],
            Oₛᵧₘ = Int[],
            Oₐₚₒₛᵧₘ = Int[],
            Oₐₚₒₛᵧₘ⁻¹ = Int[],
            Oₐ₊₋ = Int[],
            pₛᵧₘ = Float64[],
            pₐₚₒₛᵧₘ = Float64[],
            pₐₚₒₛᵧₘ⁻¹ = Float64[],
            pₐ₊₋ = Float64[],
        )
        @info "Start"
        for i ∈ 1:length(records)
            print("$i out of $(length(records))      \r")
            diff = diffs[i]
            acc = accs[i]
            date = dates[i]
            Nₛᵧₘ, Nₐₚₒₛᵧₘ, Nₐₚₒ⁻¹ₛᵧₘ = (0,0,0)
            SNP_FRAME = DataFrame(
                APOBEC = Int[],
                APOBEC⁻¹ = Int[],
                SYNONYMOUS = Int[],
            )
            for j ∈ 1:length(diff)
                if QC(diff[j]) > 4
                    # filter out mismatches with low score
                    snp_frame = alg2snp(diff[j])
                    nₛᵧₘ, nₐₚₒₛᵧₘ,nₐₚₒ⁻¹ₛᵧₘ = alg2sym(diff[j]) 
                    SNP_FRAME = vcat(SNP_FRAME, snp_frame)
                    Nₛᵧₘ += nₛᵧₘ
                    Nₐₚₒₛᵧₘ += nₐₚₒₛᵧₘ
                    Nₐₚₒ⁻¹ₛᵧₘ += nₐₚₒ⁻¹ₛᵧₘ
                end
            end
            Oₛᵧₘ = sum(SNP_FRAME.SYNONYMOUS)
            Oₐₚₒₛᵧₘ = sum(SNP_FRAME.APOBEC .* SNP_FRAME.SYNONYMOUS)
            Oₐₚₒₛᵧₘ⁻¹ = sum(SNP_FRAME.APOBEC⁻¹ .* SNP_FRAME.SYNONYMOUS)
            pₛᵧₘ = Oₛᵧₘ / Nₛᵧₘ
            pₐₚₒₛᵧₘ = Oₐₚₒₛᵧₘ / Nₐₚₒₛᵧₘ
            pₐₚₒₛᵧₘ⁻¹ = Oₐₚₒₛᵧₘ⁻¹ / Nₐₚₒ⁻¹ₛᵧₘ
            push!(info_df, [acc,
            ref_acc,
            date,
            Dates.year(date),
            Nₛᵧₘ,
            Nₐₚₒₛᵧₘ,
            Nₐₚₒ⁻¹ₛᵧₘ,
            Nₐₚₒₛᵧₘ + Nₐₚₒ⁻¹ₛᵧₘ,
            Oₛᵧₘ, 
            Oₐₚₒₛᵧₘ, 
            Oₐₚₒₛᵧₘ⁻¹, 
            Oₐₚₒₛᵧₘ + Oₐₚₒₛᵧₘ⁻¹,
            pₛᵧₘ, 
            pₐₚₒₛᵧₘ, 
            pₐₚₒₛᵧₘ⁻¹,
            pₐₚₒₛᵧₘ + pₐₚₒₛᵧₘ⁻¹,
            ])
        end
        CSV.write("./data/APOBEC_synonymous_$(ref_acc).csv", info_df)
    end

# ref_accs = [
#     "KP849470.1",
#     "MN648051.1",
#     "KJ642617",
#     # "MT903337.1",
#     # "ON563414.3",
#     # "ON676708",
#     # "OP123040",
#     # "OP257247"
# ]

@showprogress pmap(analysis, accs)


