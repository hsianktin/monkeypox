# a minimal simulated program of DNA evolution
# input: number of time, 
#    number of total synonymous mutations, number of specific APOBEC3-relevant mutations
#    number of reverse APOBEC3 mutations 
#    mutations probability of each type of mutations
#    selection coefficients for APOBEC3-relevant mutations
#    replication rate
# tracking: we track the population of each genotype/genome/DNA over time.
# output: populations of each genotype/genome/DNA over time.
# check directory
if isdir("code")
    cd("code")
end
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
n_species = 0 # number of species
# birth rate
function β(β₀::Float64, σ::Float64, smax::Float64, g::Genotype)
    return β₀ * selection(g,σ) / (smax)
end

# selection coefficient
function selection(g::Genotype, σ::Float64)
    return 1 + σ * (g.APOBEC3_mutations)
end

# death rate with selection; more APOBEC3 mutations, lower death rate
function μ(μ₀::Float64, g::Genotype, population::Int, σ::Float64, K, β₀::Float64, smax::Float64)
    return  (μ₀ + β₀* selection(g,σ) / (smax) * population / K)
end

# mutation probability


# parameters 
struct Parameters
    β₀::Float64
    μ₀::Float64
    σ::Float64
    δ::Array{Float64, 1}
    K::Function
end
# time evolving operator 
function 𝚍species╱𝚍t!(species, t, Par::Parameters)
    # birth rate 
    smax = maximum([selection(g, Par.σ) for (g, population) in species])
    total_population = ∑([population for (g, population) in species])
    total_birth_rate = ∑([β(Par.β₀, Par.σ, smax, g) * population for (g, population) in species])
    # s_max 
    # death
    total_death_rate = ∑([μ(Par.μ₀, g, total_population, Par.σ, Par.K(t), Par.β₀, smax) * population for (g, population) in species])
    # time forwarding 
    δt = randexp() / (total_birth_rate + total_death_rate)
    # determine which event happens
    # using ternary operator
    event = rand() < total_birth_rate / (total_birth_rate + total_death_rate) ? "birth" : "death"
    if event == "birth"
        # determine which genotype is born
        # using multinomial distribution
        # birth rate of each genotype
        𝛃 = [β(Par.β₀, Par.σ, smax, g) * population for (g, population) in species]
        # normalize to obtain birth probability
        Pᵦ = 𝛃 ./ ∑(𝛃)
        # sample
        g = sample([g for (g, population) in species], Weights(Pᵦ))
        # update population
        # check mutation probability for each site 
        # P𝛅 = vcat([Par.δ[1] for _ in 1:g.total_synonymous_mutations], 
        #     [Par.δ[2] for _ in 1:g.total_APOBEC3_mutations], 
        #     [Par.δ[3] for _ in 1:g.total_reverse_APOBEC3_mutations])
        # determine if at least one mutation happens
        mut_flag = false
        P𝛅 = [Par.δ[1] * g.total_synonymous_mutations, 
                Par.δ[2] * g.total_APOBEC3_mutations, 
                Par.δ[3] * g.total_reverse_APOBEC3_mutations]
        
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
                if rand() > g.reverse_APOBEC3_mutations / g.total_APOBEC3_mutations
                    g′.APOBEC3_mutations += 1
                    g′.total_APOBEC3_mutations -= 1
                    g′.total_reverse_APOBEC3_mutations += 1
                else # back mutation
                    g′.reverse_APOBEC3_mutations -= 1
                    g′.total_APOBEC3_mutations -= 1
                    g′.total_reverse_APOBEC3_mutations += 1
                end
            elseif rn < (P𝛅[1] + P𝛅[2] + P𝛅[3]) / ∑P𝛅
                if rand() > g.APOBEC3_mutations / g.total_reverse_APOBEC3_mutations
                    g′.reverse_APOBEC3_mutations += 1
                    g′.total_APOBEC3_mutations += 1
                    g′.total_reverse_APOBEC3_mutations -= 1
                else # back mutation
                    g′.APOBEC3_mutations -= 1
                    g′.total_APOBEC3_mutations += 1
                    g′.total_reverse_APOBEC3_mutations -= 1
                end
            end
        end
        if mut_flag 
            # add g′ to species
            global n_species += 1
            g′.id = n_species
            species[g′] = 1
        else
            species[g] += 1
        end
    else
        # determine which genotype dies
        # using multinomial distribution
        # death rate of each genotype
        𝛍 = [μ(Par.μ₀, g, total_population, Par.σ, Par.K(t), Par.β₀, smax) * population for (g, population) in species]
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


function normalize!(species)
    # remove extinct species whose population is 0
    for (g, population) in species
        if population <= 0
            delete!(species, g)
        end
    end
    return species
end


if length(ARGS) > 0
    profile = ARGS[1]
    include(profile)
else
    g₀ = Genotype(
     1, # id
     Int(1e6),  # total synonymous mutations
     Int(1e4),  # total APOBEC3 mutations
     Int(1e4),  # total reverse APOBEC3 mutations
     0,  # synonymous mutations
     0,  # APOBEC3 mutations
     0  # reverse APOBEC3 mutations
    )
    Par = Parameters(
    1.0,  # β₀
    0.5,  # μ₀
    1.0,  # σ
    [1e-7, 1e-5, 1e-6],  # δ
    t -> minimum([1e4,1e2 * (1+exp(0.005t))])  # K
)
end


species[g₀] = round(Int,Par.K(0.0))
n_species = 1
using ProgressBars
using Printf
# simulation
t = 0.0
T = 1000.0
sampled_DNA = []
N_t = []
time_points = []
n_sample_per_time = 10
# prog = ProgressBar(total = round(Int,100T))
# prog_count = 0
while t < T
    global t
    # @show t
    # ProgressBars.update(prog, round(Int, t))
    δt = 𝚍species╱𝚍t!(species, t, Par)
    if round(100*(t+ δt)) > round(100*t)
        # while prog_count < round(100*(t+ δt))
        #     global prog_count += 1
        #     ProgressBars.update(prog)
        # end
        smax = maximum([selection(g, Par.σ) for (g, population) in species])
        total_population = ∑([population for (g, population) in species])
        total_birth_rate = ∑([β(Par.β₀, Par.σ, smax, g) * population for (g, population) in species])
        # death
        total_death_rate = ∑([μ(Par.μ₀, g, total_population, Par.σ, Par.K(t), Par.β₀, smax) * population for (g, population) in species])    
        normalize!(species)
        # set_multiline_postfix(prog, "t = $(@sprintf("%.2f", t))\nN = $(sum([population for (g, population) in species]))\nK = $(Par.K(t))\nβ = $total_birth_rate\nμ = $total_death_rate\nN_species = $(length(species))\n")
        populations = [population for (g, population) in species]
        
        if length(populations) == 0
            break
        end
    
    end
    if round(1000*(t+δt)/T) > round(1000*t/T)
        populations = [population for (g, population) in species]
        if length(populations) == 0
            break
        end
        genotypes = [g for (g, population) in species]
        sampled_genotypes = sample(genotypes, Weights(populations), n_sample_per_time)
        for g in sampled_genotypes
            push!(sampled_DNA, g)
            push!(time_points, t+ δt)
            push!(N_t, ∑([population for (g, population) in species]))
        end
    end
    t += δt
end

# randomly sample 10 genomes from each time point in species history 
# and save them in a file
rand_label = rand(1:100)
using BSON 
# save data
label = "s_$(Par.σ)_μ₀_$(Par.μ₀)_β₀_$(Par.β₀)_Kt_$(Par.K(0.0))_$(Par.K(100.0))"
BSON.@save "data/simulated_$(label)_$(rand_label).bson" sampled_DNA time_points Par N_t

using Plots 
using LaTeXStrings
# plot number of synonymous mutations vs number of APOBEC3 mutations
# headless GKS 
ENV["GKSwstype"] = "nul"
N_synonymous = [g.synonymous_mutations for g in sampled_DNA]
N_APOBEC3 = [g.APOBEC3_mutations for g in sampled_DNA]
population = ∑([population for (g, population) in species])
p = scatter(N_synonymous, 
                N_APOBEC3, 
                xlabel = L"N_{\mathrm{ synonymous}}", 
                ylabel = L"N_{\mathrm{APOBEC3}}", 
                title = "σ = $(Par.σ), N = $population",
                legend = false, 
                markersize = 0.5, 
                dpi = 300, 
                size=(300,300))

p2 = scatter(
    N_synonymous,
    N_t,
    xlabel = L"N_{\mathrm{ synonymous}}",
    ylabel = L"\rm{population}",
    legend=false,
    markersize = 0.5,
    dpi = 300,
    size=(300,300)
)
# stack p and p2 side by side 
P = plot(p, p2, layout = (1,2), dpi = 300, size=(600,300))
savefig(P, "figures/simulated_$(label)_$(rand_label).png")
savefig(P, "figures/simulated_$(label)_$(rand_label).pdf")