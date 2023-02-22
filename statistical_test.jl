∑(x) = sum(x)
using Statistics
using LinearAlgebra
using DataFrames,CSV
using Dates 
using Pipe
# - Assume that in the animal reservoir, there is no selective pressure on APOBEC-relevant mutations.
# - Assume that all synonymous mutations are neutral and independent.

# approximation for binomial coefficients
# using Stirling's approximation
# test for homogeneity
function pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs)
    Pₐₚₒ = ∑(Oₐₚₒₛᵧₘs)/∑(Oₛᵧₘs)
    p = 1 # p-value
    for (Oₐₚₒₛᵧₘ, Oₛᵧₘ) ∈ zip(Oₐₚₒₛᵧₘs, Oₛᵧₘs)
        Ôₐₚₒₛᵧₘ = Oₛᵧₘ * Pₐₚₒ
        ∑([binom(Oₛᵧₘ, k) * Pₐₚₒ^k * (1 - Pₐₚₒ)^(Oₛᵧₘ-k) for k in 0:1:Oₛᵧₘ if abs(k - Ôₐₚₒₛᵧₘ) < abs(Oₐₚₒₛᵧₘ - Ôₐₚₒₛᵧₘ) ])
        p *= 1-∑([binom(Oₛᵧₘ, k) * Pₐₚₒ^k * (1 - Pₐₚₒ)^(Oₛᵧₘ-k) for k in 0:1:Oₛᵧₘ if abs(k - Ôₐₚₒₛᵧₘ) < abs(Oₐₚₒₛᵧₘ - Ôₐₚₒₛᵧₘ) ])
    end
    return p
end
# modified binomial function
# when n and k are large, we adopt Stirling's approximation.
function binom(n, k)
    if k > n
        return 0
    end
    if n >10 && k > 10
        return √(n/(2π*(k*(n-k)))) * exp(k*log(n/k) + (n-k)*log(n/(n-k)))
    else 
        return binomial(n, k)
    end
end

function slope(X, y)
    β = (X'X)\(X'y)
    return β
end

function random_acc_by_year(df, year)
    cur_acc = df.acc[rand([i for i in 1:length(df.acc) if Date(year+1,1,1) > df.date[i] ≥ Date(year,1,1)])]
    return cur_acc
end
function random_acc_by_year(accsₐ, times, year)
    accs = []
    for (acc, date) ∈ zip(accsₐ, times)
        if Date(year+1,1,1) > date ≥ Date(year,1,1)
            push!(accs, acc)
        end
    end
    return accs[rand(1:length(accs))]
end
years = [2017,2018,2021,2022, 2023]


# Figure 1
ref_acc = "KP849470"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
df_AC2AT = CSV.read("./data/extended_APOBEC_AC2AT_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
post2016_accs = [x for x in post2016_accs if x != ref_acc]

## Figure 1A
plot_df = @pipe df |> filter(row -> row.acc ∈ pre2016_accs, _)
slope(plot_df.Oₛᵧₘ, plot_df.Oₐ₊₋)
plot_df = @pipe df |> filter(row -> row.acc ∈ post2016_accs, _)
slope([plot_df.Oₛᵧₘ ones(nrow(plot_df))], plot_df.Oₐ₊₋)

## Figure 1B
plot_df = df_AC2AT
slope(plot_df.Oₛᵧₘ .- plot_df.Oₐ₊₋, plot_df.Oᵪ)

## Figure 1C
plot_df = df
slope(plot_df.Oₛᵧₘ .- plot_df.Oₐ₊₋, plot_df.Oᵪ)

# Figure 2
ref_acc = "KJ642617"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
post2016_accs = [x for x in post2016_accs if x != ref_acc]

## Figure 2A
plot_df = @pipe df |> filter(row -> row.acc ∈ pre2016_accs, _)
a₀ = slope(plot_df.Oₛᵧₘ, plot_df.Oₐ₊₋)

plot_df = @pipe df |> filter(row -> row.acc ∈ post2016_accs, _)
a₁,b = slope([plot_df.Oₛᵧₘ ones(nrow(plot_df))], plot_df.Oₐ₊₋)

### the intersection of the two lines
####  y = a₀ * x, y = a₁ * x + b
####  a₀ * x = a₁ * x + b
####  x = b / (a₀ - a₁)
x = b / (a₀ - a₁)
y = a₀ * x

## Figure 2B
ref_acc = "MK783029"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
post2016_accs = [x for x in post2016_accs if x != ref_acc]

plot_df = @pipe df |> filter(row -> row.acc ∈ pre2016_accs, _)
a₀, b = slope([plot_df.Oₛᵧₘ ones(nrow(plot_df))], plot_df.Oₐ₊₋)

plot_df = @pipe df |> filter(row -> row.acc ∈ post2016_accs, _)
a₁ = slope(plot_df.Oₛᵧₘ, plot_df.Oₐ₊₋)

## Test in Figure 2
ref_acc = "KJ642617"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
post2016_accs = [x for x in post2016_accs if x != ref_acc]
ps = []
for i in 1:1000
    compare_accs = [random_acc_by_year(df, year) for year in years]
    Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in compare_accs]
    Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in compare_accs]
    if isnan(pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs))
        @show compare_accs
    end
    push!(ps, pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs))
end
@show mean(ps)
#mean(ps) = 0.0009197975997642998

animal_accs = pre2016_accs ∪ ["MK783029"]
test_df = @pipe df |> filter(row -> row.acc ∈ animal_accs, _)
@show pₕ(test_df.Oₐ₊₋, test_df.Oₛᵧₘ)

ref_acc = "MK783029"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
post2016_accs = [x for x in post2016_accs if x != ref_acc]
Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in pre2016_accs]
Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in pre2016_accs]
@show pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs)

## Figure 3C
ref_acc = "KJ642617"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)   
post2016_accs = [x for x in post2016_accs if x != ref_acc]

plot_df = @pipe df |> filter(row -> row.acc ∈ pre2016_accs, _)
@show r₊₋ = slope(plot_df.Oₐₚₒₛᵧₘ, plot_df.Oₐₚₒₛᵧₘ⁻¹)
p₊ = 1 / (1 + r₊₋)
@show p₊^5
