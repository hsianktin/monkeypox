# analysis for the simulated dna mutation data.
# check directory
if isdir("code")
    cd("code")
end
# load package
using Random
using StatsBase
using BSON
using DataFrames, GLM, Statistics, HypothesisTests

∑ = sum
T = 1000.0
# state variables
# a dictionary of genotype id to population size and actual genotype in the form of number of total synonymous mutations, number of specific APOBEC3-relevant mutations, number of reverse APOBEC3 mutations
mutable struct Genotype
    id::Int
    total_synonymous_mutations::Int
    total_APOBEC3_mutations::Int
    total_reverse_APOBEC3_mutations::Int
    total_hidden_mutations::Int
    synonymous_mutations::Int
    APOBEC3_mutations::Int
    reverse_APOBEC3_mutations::Int
    hidden_mutations::Int
end

# parameters 
struct Parameters
    β₀::Float64
    μ₀::Float64
    s::Float64
    sₕ::Float64
    δ::Array{Float64, 1}
    T::Float64
    K::Function
end

function get_files(label::String, dir::String)
    # find all files of the pattern: "data/simulated_$(label)_$(rand_label).bson"
    # where rand_label is a random string of integers
    # Prepare regex pattern
    pattern = Regex("^extended_simulated_$(label)_[0-9]+.bson")
    
    # Get all files in the directory
    files = readdir(dir)
    
    # Filter the files based on the pattern
    files = filter(file -> occursin(pattern, file), files)
    
    # Prepend the directory to the file names
    files = [joinpath(dir, file) for file in files]
    
    return files
end
data_dict = Dict{Parameters, Array{Any, 1}}()
    
profiles = readdir("profiles") |> # keep only terms ending with .jl 
            (x -> filter(y -> endswith(y, ".jl"), x)) |> # filter by starting with "N"
            (x -> filter(y -> startswith(y, "extended"), x)) 

### testing 
for profile in profiles
    # load profile
    include("profiles/$profile")
    @show label = "s_$(Par.s)_sₕ_$(Par.sₕ)_μ₀_$(Par.μ₀)_β₀_$(Par.β₀)_Kt_$(Par.K(0.0))_$(Par.K(100.0))"

    dir = "data"
    files = get_files(label, dir)
    # for each of the file, parse the random id 
    DNA_samples = []
    population_sizes = []
    time_points = []
    random_ids = []
    # handle each file
    for file in files
        # load data using BSON
        try
            global data = BSON.load(file)
        catch
            rm(file)
            continue
        end
        # parse the random id
        random_id = split(file, "_")[end][1:end-5]
        if length(data[:sampled_DNA]) < 5000 # filter out early extinction
            continue
        end
        push!(DNA_samples, data[:sampled_DNA])
        push!(population_sizes, data[:N_t])
        push!(time_points, data[:time_points])
        push!(random_ids, random_id)
    end
    # build a dictionary from Par to [DNA_samples, population_sizes, time_points]
    # ADD to the dictionary
    @show [length(DNA_sample) for DNA_sample in DNA_samples]
    data_dict[Par] = [DNA_samples, population_sizes, time_points, random_ids]
end
Pars = collect(keys(data_dict))

function summarize(Pars)
    result_df = DataFrame(
        β₀ = Float64[],
        μ₀ = Float64[],
        s = Float64[],
        sₕ = Float64[],
        δ₁ = Float64[],
        δ₂ = Float64[],
        δ₃ = Float64[],
        K₀ = Float64[],
        Kₜ = Float64[],
        T = Float64[],
        coef_intercept = Float64[],
        coef_slope = Float64[],
        se_intercept = Float64[],
        se_slope = Float64[],
        p_intercept = Float64[],
        p_slope = Float64[],
        rmse = Float64[],
        r² = Float64[],
        within_group_std = Float64[],
        random_id = Int64[]
    )


    for Par_i in Pars
        T = Par_i.T
        DNA_samples_col = data_dict[Par_i][1]
        random_ids = data_dict[Par_i][4]
        DNA_samples = DNA_samples_col[1]

        for (DNA_samples, random_id) in zip(DNA_samples_col, random_ids)
            x = [gene.synonymous_mutations for gene in DNA_samples]
            y = [gene.APOBEC3_mutations for gene in DNA_samples]
            data = DataFrame(X = x, Y = y)
            # group each data point by X
            # find std of y | X 
            # then mean of std(y | X)

            # group the data by X
            grouped_data = groupby(data, :X)
            # find the std of y | X
            std_y_given_X = [std(group[!,:Y]) for group in grouped_data if nrow(group) > 1]
            # find the mean of std(y | X)
            𝔼std_yₓ = mean(std_y_given_X)
            
            # Fit the linear regression model
            lm1 = lm(@formula(Y ~ 1 + X), data)

            # Getting the coefficients
            coefficients = coef(lm1)
            coef_intercept = coefficients[1]
            coef_slope = coefficients[2]
            se_intercept = stderror(lm1)[1]
            se_slope = stderror(lm1)[2]
            p_intercept = coeftable(lm1).cols[4][1]
            p_slope = coeftable(lm1).cols[4][2]
            dof_residuals = dof_residual(lm1)
            mse = deviance(lm1) / dof_residuals
            rmse = sqrt(mse)
            r² = r2(lm1)
            # Add the results to the DataFrame
            push!(result_df, [
                Par_i.β₀,
                Par_i.μ₀,
                Par_i.s,
                Par_i.sₕ,
                Par_i.δ[1],
                Par_i.δ[2],
                Par_i.δ[3],
                Par_i.K(0.0),
                Par_i.K(T),
                T,
                coef_intercept,
                coef_slope,
                se_intercept,
                se_slope,
                p_intercept,
                p_slope,
                rmse,
                r²,
                𝔼std_yₓ,
                parse(Int64, random_id)
            ])
            end
    end
    return result_df
end

result_df = summarize(Pars)
# filter out result_df with s == 0
result_df = filter(row -> row.s != 0, result_df)
using CSV
CSV.write("data/extended_dna_simulation_result_df.csv", result_df)
# filter result_df by selecting those K₀ == Kₜ
constant_K_result_df = filter(row -> row.K₀ == row.Kₜ, result_df)
CSV.write("data/extended_dna_simulation_result_df_constant_K.csv", constant_K_result_df)
# non-constant K
non_constant_K_result_df = filter(row -> row.K₀ != row.Kₜ, result_df)
CSV.write("data/extended_dna_simulation_result_df_non_constant_K.csv", non_constant_K_result_df)

# find mean and stdd of coef_slope and coef_intercept conditioned on the same Par 
# group the result_df by Par
grouped_result_df = groupby(result_df, [:β₀, :μ₀, :s, :sₕ, :δ₁, :δ₂, :δ₃, :K₀, :Kₜ, :T])

# find the mean and std of coef_slope and coef_intercept conditioned on the same Par
mean_coef_slope = [mean(group.coef_slope) for group in grouped_result_df]
std_coef_slope = [std(group.coef_slope) for group in grouped_result_df]
mean_coef_intercept = [mean(group.coef_intercept) for group in grouped_result_df]
std_coef_intercept = [std(group.coef_intercept) for group in grouped_result_df]

# new dataframe
stat_df = DataFrame(
    β₀ = [group.β₀[1] for group in grouped_result_df],
    μ₀ = [group.μ₀[1] for group in grouped_result_df],
    s = [group.s[1] for group in grouped_result_df],
    sₕ = [group.sₕ[1] for group in grouped_result_df],
    δ₁ = [group.δ₁[1] for group in grouped_result_df],
    δ₂ = [group.δ₂[1] for group in grouped_result_df],
    δ₃ = [group.δ₃[1] for group in grouped_result_df],
    K₀ = [group.K₀[1] for group in grouped_result_df],
    Kₜ = [group.Kₜ[1] for group in grouped_result_df],
    T = [group.T[1] for group in grouped_result_df],
    mean_coef_intercept = mean_coef_intercept,
    std_coef_intercept = std_coef_intercept,
    mean_coef_slope = mean_coef_slope,
    std_coef_slope = std_coef_slope
)

# write the stat_df to csv, all, constant K, non-constant K
CSV.write("data/extended_dna_simulation_stat_df.csv", stat_df)
constant_K_stat_df = filter(row -> row.K₀ == row.Kₜ, stat_df)
CSV.write("data/extended_dna_simulation_stat_df_constant_K.csv", constant_K_stat_df)
non_constant_K_stat_df = filter(row -> row.K₀ != row.Kₜ, stat_df)
CSV.write("data/extended_dna_simulation_stat_df_non_constant_K.csv", non_constant_K_stat_df)


# # PGFPlotsX plot scatter plot of result_df.coef_slope vs result_df.s, 
# # with color result_df.Kₜ,

# using PGFPlotsX
# using LaTeXStrings
# using Statistics
# using StatsBase

# # plot the result_df
# @pgf Axis(
#     {
#         xlabel = L"s",
#         ylabel = L"\beta",
#         colorbar,
#         colormap = "viridis",
#     },
#     Plot(
#         {
#             only_marks,
#             scatter,
#             scatter_src = "explicit",
#             mark = "*",
#         },
#         Table(result_df, x = "s", y = "coef_slope", meta = "Kₜ")
#     )
# )

