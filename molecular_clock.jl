# Molecular Clock Analysis
∑(x) = sum(x)
using Statistics
using LinearAlgebra
using DataFrames,CSV
using Dates 
using Pipe
ref_acc = "KJ642617"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
post2016_accs = [x for x in post2016_accs if x != ref_acc]

function slope(X, y)
    β = (X'X)\(X'y)
    return β
end

## Figure 2A
plot_df = @pipe df |> filter(row -> row.acc ∈ pre2016_accs, _)
a₀ = slope(plot_df.Oₛᵧₘ, plot_df.Oₐ₊₋)

plot_df = @pipe df |> filter(row -> row.acc ∈ post2016_accs, _)
a₁,b = slope([plot_df.Oₛᵧₘ ones(nrow(plot_df))], plot_df.Oₐ₊₋)

### the intersection of the two lines
y⁺ = b / (a₀ - a₁)

df = CSV.read("./data/lineage_based_clock.csv", DataFrame)

function lm(X,Y)
    β = X\Y
    ŷ = X*β
    ȳ = mean(Y)
    SST = sum((Y .- ȳ).^2)
    SSE = sum((ŷ - Y).^2)
    R² = 1-SSE/SST
    return β, R²
end



x = df.date .- df.date[1]
# convert to year
x = x .|> (x -> (Dates.value.(x))/365.25)
X = [ones(length(x)) x]
y = df.Oₛᵧₘ

@show β, R² = lm(X,y)
# find [1 x⁺] such that [1 x⁺]β = y⁺
# i.e. 1 * β[1] + x⁺ * β[2] = y⁺
# i.e. x⁺ = (y⁺ - β[1]) / β[2]
x⁺ = (y⁺ - β[1]) / β[2]
Δx = x⁺ * 365.25 |> round |> Dates.Day
date⁺ = df.date[1] + Δx
dates = Date.([
    i for i in 2016:2023])
xs = [Dates.value(i - df.date[1])/365.25 for i in dates]
ys = [β[1] + β[2]*i for i in xs]
y⁺s = [β[1] + β[2]*x⁺ for i in xs]

CSV.write("./data/lineage_based_clock_model.csv", DataFrame(
    date = dates,
    Oₛᵧₘ = ys,
    Oₛᵧₘ⁺ = y⁺s,
    ))
