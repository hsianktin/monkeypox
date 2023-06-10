# a minimal simulated program of DNA evolution
# input: number of time, 
#    number of total synonymous mutations, number of specific APOBEC3-relevant mutations
#    number of reverse APOBEC3 mutations 
#    mutations probability of each type of mutations
#    selection coefficients for APOBEC3-relevant mutations
#    replication rate
# tracking: we track the population of each genotype/genome/DNA over time.
# output: populations of each genotype/genome/DNA over time.

# load package
using Random
using StatsBase
∑ = sum
# state variables
# a dictionary of genotype id to population size and actual genotype in the form of number of total synonymous mutations, number of specific APOBEC3-relevant mutations, number of reverse APOBEC3 mutations
mutable struct Genotype
    id::Int
    total_synonymous_mutations::Int
    total_APOBEC3_mutations::Int
    total_reverse_APOBEC3_mutations::Int
    synonymous_mutations::Int
    APOBEC3_mutations::Int
    reverse_APOBEC3_mutations::Int
end

species = Dict{Genotype, Int}() # map genotype to population size

# birth rate
function β(β₀::Float64)
    return β₀ 
end

# selection coefficient
function selection(g::Genotype, s::Float64)
    return 1 + s * (g.APOBEC3_mutations)
end

# death rate with selection; more APOBEC3 mutations, lower death rate
function μ(μ₀::Float64, g::Genotype, population::Int, s::Float64, K, β₀::Float64)
    return  (μ₀ + β₀ * population / K) / selection(g, s)
end

# mutation probability


# parameters 
struct Parameters
    β₀::Float64
    μ₀::Float64
    s::Float64
    δ::Array{Float64, 1}
    K::Function
end
# time evolving operator 
function 𝚍species╱𝚍t!(species, t, Par::Parameters)
    # birth rate 
    total_birth_rate = ∑([β(Par.β₀) * population for (g, population) in species])
    # death
    total_death_rate = ∑([μ(Par.μ₀, g, population, Par.s, Par.K(t), Par.β₀) * population for (g, population) in species])
    # time forwarding 
    δt = randexp() / (total_birth_rate + total_death_rate)
    # determine which event happens
    # using ternary operator
    event = rand() < total_birth_rate / (total_birth_rate + total_death_rate) ? "birth" : "death"
    if event == "birth"
        # determine which genotype is born
        # using multinomial distribution
        # birth rate of each genotype
        𝛃 = [β(Par.β₀) * population for (g, population) in species]
        # normalize to obtain birth probability
        Pᵦ = 𝛃 ./ ∑(𝛃)
        # sample
        g = sample([g for (g, population) in species], Weights(Pᵦ))
        # update population
        # check mutation probability for each site 
        P𝛅 = vcat([Par.δ[1] for _ in 1:g.total_synonymous_mutations], 
            [Par.δ[2] for _ in 1:g.total_APOBEC3_mutations], 
            [Par.δ[3] for _ in 1:g.total_reverse_APOBEC3_mutations])
        rand_numbers = rand(length(P𝛅))
        # determine if at least one mutation happens
        mut_flag = false
        
        g′ = deepcopy(g)
        # The following implementation is very inefficient when total sites are large 
        # because rng is called for each site
        # for i in 1:g.total_synonymous_mutations
        #     if rand_numbers[i] < P𝛅[i]
        #         mut_flag = true
        #         if rand() > g.synonymous_mutations / g.total_synonymous_mutations
        #             g′.synonymous_mutations += 1
        #         else # back mutation
        #             g′.synonymous_mutations -= 1
        #         end
        #     end
        # end
        # for i in g.total_synonymous_mutations + 1:g.total_synonymous_mutations + g.total_APOBEC3_mutations
        #     if rand_numbers[i] < P𝛅[i]
        #         mut_flag = true
        #         if rand() > g.APOBEC3_mutations / g.total_APOBEC3_mutations
        #             g′.APOBEC3_mutations += 1
        #         else # back mutation
        #             g′.APOBEC3_mutations -= 1
        #         end
        #     end
        # end
        # for i in g.total_synonymous_mutations + g.total_APOBEC3_mutations + 1:g.total_synonymous_mutations + g.total_APOBEC3_mutations + g.total_reverse_APOBEC3_mutations
        #     if rand_numbers[i] < P𝛅[i]
        #         mut_flag = true
        #         if rand() > g.reverse_APOBEC3_mutations / g.total_reverse_APOBEC3_mutations
        #             g′.reverse_APOBEC3_mutations += 1
        #         else # back mutation
        #             g′.reverse_APOBEC3_mutations -= 1
        #         end
        #     end
        # end

        # we assume that total mutation probability is small, so we use the following implementation
        # assuming that at most one mutation happens
        ∑P𝛅 = ∑(P𝛅)
        if rand() < ∑P𝛅
            mut_flag = true
            rn = rand()
            if rn < P𝛅[1] / ∑P𝛅
                if rand() > g.synonymous_mutations / g.total_synonymous_mutations
                    g′.synonymous_mutations += 1
                else # back mutation
                    g′.synonymous_mutations -= 1
                end
            elseif rn < (P𝛅[1] + P𝛅[2]) / ∑P𝛅
                if rand() > g.APOBEC3_mutations / g.total_APOBEC3_mutations
                    g′.APOBEC3_mutations += 1
                else # back mutation
                    g′.APOBEC3_mutations -= 1
                end
            else
                if rand() > g.reverse_APOBEC3_mutations / g.total_reverse_APOBEC3_mutations
                    g′.reverse_APOBEC3_mutations += 1
                else # back mutation
                    g′.reverse_APOBEC3_mutations -= 1
                end
            end
        end
        if mut_flag 
            # add g′ to species
            g′.id = length(species) + 1
            species[g′] = 1
        else
            species[g] += 1
        end
    else
        # determine which genotype dies
        # using multinomial distribution
        # death rate of each genotype
        𝛍 = [μ(Par.μ₀, g, population, Par.s, Par.K(t), Par.β₀) * population for (g, population) in species]
        # normalize
        P𝛍 = 𝛍 / sum(𝛍)
        # sample
        g = sample([g for (g, population) in species], Weights(P𝛍))
        # update population
        species[g] -= 1
    end
    return δt
end

# initial condition
species = Dict{Genotype, Int}()
g₀ = Genotype(
     1, # id
     Int(1e6),  # total synonymous mutations
     Int(1e4),  # total APOBEC3 mutations
     Int(1e4),  # total reverse APOBEC3 mutations
     0,  # synonymous mutations
     0,  # APOBEC3 mutations
     0  # reverse APOBEC3 mutations
)

if length(ARGS) > 0
    profile = ARGS[1]
    include(profile)
else
    Par = Parameters(
    1.0,  # β₀
    0.5,  # μ₀
    1.0,  # s
    [1e-8, 1e-5, 1e-9],  # δ
    t -> minimum([1e7 ,1e5 * (1+exp(t))])  # K
)
end


species[g₀] = 100
using ProgressBars
using Printf
# simulation
t = 0.0
T = 100.0
species_history = []
time_history = []
prog = ProgressBar(total = round(Int,100T))
while t < T
    # @show t
    # ProgressBars.update(prog, round(Int, t))
    push!(species_history, deepcopy(species))
    push!(time_history, t)
    δt = 𝚍species╱𝚍t!(species, t, Par)
    if round(100*(t+ δt)) > round(100*t)
        ProgressBars.update(prog)
        set_postfix(prog, t = @sprintf("%.2f", t))
    end
    t += δt
end
push!(species_history, deepcopy(species))
push!(time_history, t)

# randomly sample 10 genomes from each time point in species history 
# and save them in a file
sampled_DNA = []
time_points = []
n_sample_per_time = 1
for (species, time) in zip(species_history, time_history)
    populations = [population for (g, population) in species]
    genotypes = [g for (g, population) in species]
    sampled_genotypes = sample(genotypes, Weights(populations), n_sample_per_time)
    for g in sampled_genotypes
        push!(sampled_DNA, g)
        push!(time_points, time)
    end
end

using BSON 
# save data
label = "s_$(Par.s)_μ₀_$(Par.μ₀)_β₀_$(Par.β₀)_Kt_$(Par.K(0.0))_$(Par.K(100.0))"
BSON.@save "data/simulated_$(label)_$(rand(1:100)).bson" sampled_DNA time_points Par
