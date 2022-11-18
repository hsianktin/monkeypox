# Molecular Clock Analysis
∑(x) = sum(x)
using Statistics
using LinearAlgebra
using DataFrames,CSV
using Dates 
using Pipe

function linear(accs, ref_acc)
    df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)
    Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
    Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]
    y = Oₛᵧₘs
    X = Oₐₚₒₛᵧₘs
    β = X\y 
    ŷ = X*β
    ȳ = mean(y)
    SST = sum((y .- ȳ).^2)
    SSE = sum((ŷ - y).^2)
    R² = 1-SSE/SST
    return (β, R²)
end


function affine(accs, ref_acc)
    df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)
    Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
    Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]
    X = [ones(length(Oₐₚₒₛᵧₘs)) Oₐₚₒₛᵧₘs]
    y = Oₛᵧₘs
    β = X\y 
    ŷ = X*β
    ȳ = mean(y)
    SST = sum((y .- ȳ).^2)
    SSE = sum((ŷ - y).^2)
    R² = 1-SSE/SST
    return (β[1], β[2], R²)
end

include("grouping.jl")
ref_acc = "KJ642617"
df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)

# subdataframe containing only accsₐ
dfₐ = @pipe df |> filter(row -> row.acc ∈ accsₐ, _)
CSV.write("./data/APOBEC_synonymous_$(ref_acc)_APOBEC_group.csv", dfₐ)


include("grouping.jl")
ref_acc = "KJ642617"
β₀, β₁, R² = affine(accsₐ, ref_acc)
β′, R²′ = linear(accsᶜ, ref_acc)

# find intersection of y = β₀ + β₁x and y = β′₁x
x₀ = (β₀)/(β′-β₁)
@show y₀ = x₀ * β′

data = @pipe CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame) |>
        filter(row -> row.acc ∈ accsₐ, _)
x = data.Oₛᵧₘ

X = [ones(length(x)) x]
y = data.year

function lm(X,Y)
    β = X\Y
    ŷ = X*β
    ȳ = mean(Y)
    SST = sum((Y .- ȳ).^2)
    SSE = sum((ŷ - Y).^2)
    R² = 1-SSE/SST
    return β, R²
end

@show β, R² = lm(X,y)
@show year = floor(Int, β[1] + β[2]*y₀)
@show month = round(Int, (β[1] + β[2]*y₀ - year)*12)

