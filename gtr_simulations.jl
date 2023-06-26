using StatsBase
using LinearAlgebra

# accept id as argument
if length(ARGS) < 1
    println("Please provide an id")
    exit(1)
else
    id = parse(Int, ARGS[1])
end

# Define the equilibrium base frequencies.
pi_A = 0.25
pi_G = 0.25
pi_C = 0.25
pi_T = 0.25

# Define the transition rate parameters.
alpha = 1.0  # A <-> G
beta = 1.0  # A <-> C
gamma = 1.0  # A <-> T
delta = 1.0  # G <-> C
epsilon = 1.0  # G <-> T
eta = 1.0  # C <-> T

# Define the rate matrix.
Q = [
    -(alpha*pi_G + beta*pi_C + gamma*pi_T) alpha*pi_G beta*pi_C gamma*pi_T;
    alpha*pi_A -(alpha*pi_A + delta*pi_C + epsilon*pi_T) delta*pi_C epsilon*pi_T;
    beta*pi_A delta*pi_G -(beta*pi_A + delta*pi_G + eta*pi_T) eta*pi_T;
    gamma*pi_A epsilon*pi_G eta*pi_C -(gamma*pi_A + epsilon*pi_G + eta*pi_C);
]

# Normalize the rate matrix so that the average substitution rate is 1.
# Q = Q / sum(Q)

# Define the equilibrium distribution.
pi = [pi_A, pi_G, pi_C, pi_T]


# Define the DNA bases.
bases = ['A', 'G', 'C', 'T']

# Function to sample a base from the equilibrium distribution.
function sample_base()
    base = sample(bases, Weights(pi))
    return base
end

using Random

function simulate_transition(base)
    # Find the index of the current base.
    i = findfirst(isequal(base), bases)
    
    # Compute the transition probabilities.
    p = [Q[i, j] ./ sum([Q[i, j] for j in 1:4 if j != i]) for j in 1:4 if j != i]
    
    sample_bases = [bases[j] for j in 1:4 if j != i]

    # Sample a new base from the transition probabilities.
    new_base = sample(sample_bases, Weights(p))
    
    # Compute the rate of the exponential distribution.
    rate = sum([Q[i, j] for j in 1:4 if j != i])
    
    # Generate an exponentially-distributed random number.
    random_number = randexp() / rate
    
    return new_base, random_number
end

function simulate_sequence(length)
    # Initialize the sequence with a random base.
    sequence = [sample_base()]
    
    # Initialize the list of accumulated waiting times with the first entry set to 0.
    accumulated_times = [0.0]
    
    # Simulate the rest of the sequence.
    for i in 2:length
        new_base, waiting_time = simulate_transition(sequence[end])
        push!(sequence, new_base)
        
        # Add the waiting time to the last accumulated time to get the new accumulated time.
        new_accumulated_time = accumulated_times[end] + waiting_time
        push!(accumulated_times, new_accumulated_time)
    end
    
    return sequence, accumulated_times
end

using ProgressBars
using DataFrames, CSV

# Define the number of sequences to simulate.
N = 100000

# Preallocate a vector of DataFrames.
dfs = [DataFrame(id = Int[], nucleotide = String[], time = Float64[]) for _ in 1:Threads.nthreads()]

# Use multithreading to simulate the sequences.
Threads.@threads for i in ProgressBar(1:N)
    # Get the thread id.
    thread_id = Threads.threadid()
    
    # Simulate a sequence.
    sequence, times = simulate_sequence(100)
    
    # Append the sequence data to the DataFrame.
    for (nucleotide, time) in zip(sequence, times)
        push!(dfs[thread_id], (i, string(nucleotide), time))
    end
    
end

# Vertically concatenate the DataFrames.
df = vcat(dfs...)

sort!(df, [:time, :id])
# Print the DataFrame.
# print(df)

if isdir("code")
    cd("code")
end
# Write the DataFrame to a CSV file, compressed with gzip.
CSV.write("data/gtr_simulated_sequence_$(id).csv.gz", df, gzip=true)