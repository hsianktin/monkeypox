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
using CSV

∑ = sum
T = 1000.0
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

struct Parameters
    β₀::Float64
    μ₀::Float64
    s::Float64
    δ::Array{Float64, 1}
    K::Function
end

function get_files(label::String, dir::String)
    # find all files of the pattern: "data/simulated_$(label)_$(rand_label).bson"
    # where rand_label is a random string of integers
    # Prepare regex pattern
    pattern = Regex("^simulated_$(label)_[0-9]+.bson")
    
    # Get all files in the directory
    files = readdir(dir)
    
    # Filter the files based on the pattern
    files = filter(file -> occursin(pattern, file), files)
    
    # Prepend the directory to the file names
    files = [joinpath(dir, file) for file in files]
    
    return files
end

# get file by label, dir, and id
function get_file(label::String, dir::String, id::Int)
    # find all files of the pattern: "data/simulated_$(label)_$(rand_label).bson"
    # where rand_label is a random string of integers
    # Prepare regex pattern
    pattern = Regex("^simulated_$(label)_$(id).bson")
    
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
            (x -> filter(y -> startswith(y, "N"), x)) 

function get_profile(N, s; label="")
    return "profiles/N$(N)_s$(s)$(label).jl"
end

profiles = [get_profile(100, 0.001; label="_constant_K"), 
            get_profile(100, 0.1; label="_constant_K"),
            get_profile(1000, 0.1, label="_constant_K")]
Pars = []
labels = []
for profile in profiles
    include(profile)
    push!(labels, "s_$(Par.s)_μ₀_$(Par.μ₀)_β₀_$(Par.β₀)_Kt_$(Par.K(0.0))_$(Par.K(100.0))")
    push!(Pars, deepcopy(Par))
end

ids = [33,27,26]

using PGFPlotsX
# replace \\ with /
pgfplots_str_norm(str) = replace("$(pwd())/$(str)", "\\" => "/")

function read_sample_fit(label, id)
    file = get_file(label, "data", id)[1]

    data = BSON.load(file)
    DNA_samples = data[:sampled_DNA]
    population_size = data[:N_t]
    time_points = data[:time_points]

    DNA_samples 
    n_syn = [dna.synonymous_mutations for dna in DNA_samples ]
    n_APOBEC3 = [dna.APOBEC3_mutations for dna in DNA_samples ]
    data = DataFrame(n_syn=n_syn, n_APOBEC3=n_APOBEC3)
    # group by n_syn 
    grouped = groupby(data, :n_syn)
    # within each group, randomly pick 1 sample
    # sampled = combine(grouped, :n_APOBEC3 => (x -> rand(x)) => :n_APOBEC3)
    # within each group, randomly pick 10 samples with repeat 
    sampled = combine(grouped, :n_APOBEC3 => (x -> sample(x, 10, replace=true)) => :n_APOBEC3)
    # save sampled data to external csv file for PGFPlots to Call 
    CSV.write("data/sampled_data_$(label).csv", sampled)
    lm1 = lm(@formula(n_APOBEC3 ~ 0+ n_syn), data)
    # no intercept fitting 
    return ("data/sampled_data_$(label).csv", coef(lm1)[1])
end


plot_info = [read_sample_fit(label, id) for (label, id) in zip(labels, ids)]
plot_info[1][1]
p = @pgf Axis(
    {
        xlabel = "number of synonymous mutations \$n - n_{\\rm a}\$",
        ylabel = "number of APOBEC3-induced mutations \$n_{\\rm a}\$",
        width = "10cm",
        height = "10cm",
        # set legend pos to north west
        legend_pos = "north west",
        legend_style = {draw = "none", fill = "none", align = "left"},
        legend_cell_align = "left",
        xmin = 0,
        ymin = 0,
    },
    # Plot 1
    Plot(
        {
            only_marks,
            # set mark size to 1.5pt
            mark_size = "1pt",
            # let the mark to be semi transparent
            mark_options = {fill_opacity = 0.1},
            color = "red!50!black"
        },
        Table(
            {
                x = "n_syn",
                y = "n_APOBEC3",
                col_sep = "comma"
            },
            plot_info[1][1]  |> pgfplots_str_norm
        )
        # use "data/extended_sampled_data.csv"
    ),
    # add a linear fit based on lm1 data using slope and intercept
    # rather than table data points 
    Plot(
        {
            no_marks,
            domain = "0:100",
            color = "red",
            ultra_thick,
            # no legeng by forget plot
            forget_plot,
        },
        Expression(
            "x*$(plot_info[1][2]) "
        ),

    ),
    # legend entry, s = $(s), sₕ = $(sₕ), K(t) = $(K(0.0)); as well as linear fit 
    LegendEntry(
        "s = $(Pars[1].s), \$s_{\\rm h}\$ = 0, K(t) = $(Pars[1].K(0.0))"
    ),
    # Plot 2
    Plot(
        {
            only_marks,
            # set mark size to 1.5pt
            mark_size = "1pt",
            # let the mark to be semi transparent
            mark_options = {fill_opacity = 0.1},
            color = "green!50!black"
        },
        Table(
            {
                x = "n_syn",
                y = "n_APOBEC3",
                col_sep = "comma"
            },
            plot_info[2][1] |> pgfplots_str_norm
        )
        # use "data/extended_sampled_data.csv"
    ),
    # add a linear fit based on lm1 data using slope and intercept
    # rather than table data points 
    Plot(
        {
            no_marks,
            domain = "0:100",
            color = "green",
            ultra_thick,
            # no legeng by forget plot
            forget_plot,
        },
        Expression(
            "x*$(plot_info[2][2]) "
        ),

    ),
    # legend entry, s = $(s), sₕ = $(sₕ), K(t) = $(K(0.0)); as well as linear fit 
    LegendEntry(
        "s = $(Pars[2].s), \$s_{\\rm h}\$ = 0, K(t) = $(Pars[2].K(0.0))"
    ),
    # Plot 3
    Plot(
        {
            only_marks,
            # set mark size to 1.5pt
            mark_size = "1pt",
            # let the mark to be semi transparent
            mark_options = {fill_opacity = 0.1},
            color = "blue!50!black"
        },
        Table(
            {
                x = "n_syn",
                y = "n_APOBEC3",
                col_sep = "comma"
            },
            plot_info[3][1] |> pgfplots_str_norm
        )
        # use "data/extended_sampled_data.csv"
    ),
    # add a linear fit based on lm1 data using slope and intercept
    # rather than table data points 
    Plot(
        {
            no_marks,
            domain = "0:100",
            color = "blue",
            ultra_thick,
            # no legeng by forget plot
            forget_plot,
        },
        Expression(
            "x*$(plot_info[3][2]) "
        ),

    ),
    # legend entry, s = $(s), sₕ = $(sₕ), K(t) = $(K(0.0)); as well as linear fit 
    LegendEntry(
        "s = $(Pars[3].s), \$s_{\\rm h}\$ = 0, K(t) = $(Pars[3].K(0.0))"
    ),
) 
pgfsave("dna_plot.tex", p)