using StatsBase
using LinearAlgebra

using DataFrames,CSV
if isdir("code")
    cd("code")
end

# accept id as argument
if length(ARGS) < 1
    @warn "Please provide an id"
    @warn "setting id = 1"
    id = 1
elseif length(ARGS) == 1
    id = parse(Int, ARGS[1])
else
    ids = parse.(Int, ARGS)
end


df = CSV.read("data/gtr_simulated_sequence_$(id).csv.gz", DataFrame)

# step 1: get the original sequence
df₀ = df[df.time .== 0.0, :]
seq₀ = [n for n in df₀.nucleotide]

# step 2: get sequence at other times 
sample_times = [i for i in 1:1:101] * 5e-4
function get_seq_by_t(df, t)
    seq = similar(seq₀)
    df_t = df[df.time .≤ t, :]
    for i in 1:nrow(df_t)
        seq[df_t.id[i]] = df_t.nucleotide[i]
    end
    return seq
end
sample_seqs = [get_seq_by_t(df, t) for t in sample_times]

# step 3: compare the sequences to the original sequence, and count the number of differences, including the binucleotide context
# return a 16×16 matrix, (AA, AT, ..., GG) ⊗ (AA, AT, ..., GG) -> Int = # of differences from X -> Y
# Scan the sequence, aware of the boundary treatment 
# Define a function to count the number of differences between two sequences
function count_differences(seq1, seq2)
    # Define a dictionary to store the counts of each binucleotide context
    counts = Dict{String, Int}()
    for i in 1:length(seq1)-1
        # Get the binucleotide context for each sequence
        context1 = seq1[i:i+1] |> join
        context2 = seq2[i:i+1] |> join
        # Increment the count for the corresponding binucleotide context
        if context1 != context2
            counts["$context1->$context2"] = get(counts, "$context1->$context2", 0) + 1
        end
    end
    # Return the counts dictionary
    return counts
end
letters = ["A", "C", "G", "T"]
dinucleotides = [l1*l2 for l1 in letters, l2 in letters] |> vec

# Define a function to compare all sample sequences to the original sequence
function compare_sequences(sample_seqs, seq₀)
    matrices = []
    # Loop over each sample sequence
    for sample_seq in sample_seqs
        # Count the differences between the sample sequence and the original sequence
        counts = count_differences(seq₀, sample_seq)
        # Increment the counts dictionary with the sample counts
        for (context, count) in counts
            counts[context] = get(counts, context, 0) + count
        end
        # Convert the counts dictionary to a 16x16 matrix
        matrix = zeros(Int, 16, 16)
        for (i, row_nuc) in enumerate(dinucleotides)
            for (j, col_nuc) in enumerate(dinucleotides)
                context = "$row_nuc->$col_nuc"
                matrix[i, j] = get(counts, context, 0)
            end
        end
        # Append the matrix to the list of matrices
        push!(matrices, matrix)
    end
    # Return the matrix
    return matrices
end

# Call the compare_sequences function with the sample sequences and the original sequence
matrices = compare_sequences(sample_seqs, seq₀)

# step 4: enumerate the number of "*A -> *T" mutations for each matrix in the list
#    1) get the index of the dinucleotide context
_A = findall(x -> x[2] == 'A', dinucleotides)
_C = findall(x -> x[2] == 'C', dinucleotides)
_G = findall(x -> x[2] == 'G', dinucleotides)
_T = findall(x -> x[2] == 'T', dinucleotides)
#    2) get the number of mutations where input is (matrix, start_end_pairs)
function count_mutations_by_matrix(matrix, start_end_pairs)
    # Initialize the count to 0
    count = 0
    # Loop over the start-end pairs
    for (start, stop) in start_end_pairs
        # Increment the count by the number of mutations in the matrix
        count += matrix[start, stop]
    end
    # Return the count
    return count
end

TT = findall(x -> x[1] == 'T' && x[2] == 'T', dinucleotides)
TC = findall(x -> x[1] == 'T' && x[2] == 'C', dinucleotides)
A_to_T = [count_mutations_by_matrix(
                matrix, 
                [(i,j) for i in _A, j in _T]
            ) 
            for matrix in matrices]

A_to_X = [count_mutations_by_matrix(
                matrix, 
                [(i,j) for i in _A, j in union(_C, _G, _T)]
            ) 
            for matrix in matrices]
TC_to_TT = [count_mutations_by_matrix(
                matrix, 
                [(i,j) for i in TC, j in TT]
            ) 
            for matrix in matrices]

x_y = union(
    [(a,y) for a in _A, y in union(_C, _G, _T)],
    [(t,y) for t in _T, y in union(_A, _C, _G)],
    [(c,y) for c in _C, y in union(_A, _G, _T)],
    [(g,y) for g in _G, y in union(_A, _C, _T)]
)
            
X_to_Y = [count_mutations_by_matrix(
                matrix, 
                x_y
            ) 
            for matrix in matrices]
# # Plotting
# using Plots
# gr()
# using LaTeXStrings

# plot(
#     X_to_Y, A_to_X,
#     xlabel = L"\# X \rightarrow Y", 
#     ylabel = L"\# A \rightarrow X", 
#     legend = false,
#     size = (400, 400),
#     marker = :circle,
# )

# plot(
#     X_to_Y, TC_to_TT,
#     xlabel = L"\# X \rightarrow Y",
#     ylabel = L"\# TT \rightarrow TC",
#     size = (400, 400),
#     marker = :circle,
# )

CSV.write("data/gtr_analysis_$(id).csv", DataFrame(
    X_to_Y = X_to_Y,
    A_to_X = A_to_X,
    TC_to_TT = TC_to_TT
))




