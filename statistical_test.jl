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


include("grouping.jl")


# test for binomial distribution
function pᵦ(accs, ref_acc)
    df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)
    Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
    Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]
    return pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs)
end

function pᵦ₊(accs, ref_acc)
    df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)
    Oₐₚₒₛᵧₘs = [df.Oₐₚₒₛᵧₘ[df.acc .== acc][1] for acc in accs]
    Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]
    return pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs)
end


## Test 1: animal reservoir test
ref_acc = "KJ642617"
df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)

accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
accs = [x for x in accs if x != ref_acc]
Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]
@info "p-value for the reference group having the same evolutionary environment" pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs)

# use linear regression on Oₐₚₒₛᵧₘs and Oₛᵧₘs
X = Oₛᵧₘs
y = Oₐₚₒₛᵧₘs
β = X\y 
@show β
# @show β * df.Nₛᵧₘ[1]/df.Nₐₚₒₛᵧₘ[1]
# @show mean(y - X*β)
# @show std(y - X*β)
# calculate R²
ŷ = X*β
ȳ = mean(y)
SST = sum((y .- ȳ).^2)
SSE = sum((ŷ - y).^2)
# @show 1-SSE/SST
@info "R² for the reference group" 1-SSE/SST


## Test 2: human reservoir test
ref_acc = "MN648051"
# Linear regression
df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)
# Now we construct the test for homogeneity
# we pick KJ642617.1, KJ642617.1, JX878428.1
old_accs = accs
accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
accs = [x for x in accs if x != ref_acc]
Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]
years = [2017,2018,2021,2022]
ps = []
for i in 1:100
    compare_accs = [random_acc_by_year(df, year) for year in years]
    Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in compare_accs]
    Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in compare_accs]
    if isnan(pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs))
        @show compare_accs
    end
    push!(ps, pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs))
end
@info "mean p-value for the APOBEC3-group having the same evolutionary environment" mean(ps)

# @show pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs)
# the p-test fail because we did not take into account the correlation

# use linear regression on Oₐₚₒₛᵧₘs and Oₛᵧₘs
# X = Oₛᵧₘs
# y = Oₐₚₒₛᵧₘs
# β = X\y
# @show β
# @show β * df.Nₛᵧₘ[1]/df.Nₐₚₒₛᵧₘ[1]
# @show mean(y - X*β)
# @show std(y - X*β)
# calculate R²
ŷ = X*β
ȳ = mean(y)
SST = sum((y .- ȳ).^2)
SSE = sum((ŷ - y).^2)
@info "R² for the linear model of the APOBEC3 group"  R² = 1-SSE/SST
# Statistical Test by Sampling
df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)
# because the sequences in 2022 are not independent, but highly correlated, we just pick one sequence from each cluster
# randomly pick one sequence whose date is after 2022-01-01

# function random_acc_by_year(df, year)
#     cur_acc = df.acc[rand([i for i in 1:length(df.acc) if Date(year+1,1,1) > df.date[i] ≥ Date(year,1,1)])]
#     return cur_acc
# end

# Test 2.5 : test if all the genommes are from the same evolutionary environment
ref_acc = "KJ642617"
df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)

accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)

years = [2017,2018,2021,2022]
accs = [x for x in accs if x != ref_acc]
accs =[accs; [random_acc_by_year(df, year) for year in years]]
Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]
@info "p-value for the all genomes having the same evolutionary environment" pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs)

Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accsₐ]
Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accsₐ]
X = [ones(length(Oₐₚₒₛᵧₘs)) Oₛᵧₘs]
y = Oₐₚₒₛᵧₘs
β = X\y
@info "β for the linear model of the APOBEC3 group" β

# R\^2  
ŷ = X*β
ȳ = mean(y)
SST = sum((y .- ȳ).^2)
SSE = sum((ŷ - y).^2)
@info "R² for the linear model of the APOBEC3 group" 1-SSE/SST


# Test 3: test if KJ642617 is common ancestor of human reservoir
ref_acc = "KJ642617"
years = [2017,2018,2021,2022]
ps = []
for i in 1:100
    accs_rand = [random_acc_by_year(accsₐ, dates, year) for year in years]
    push!(ps,  pᵦ₊(accs_rand, ref_acc))
end
@info "mean p-value for the null hypothesis KJ642617 is the ancestor of APOBEC3-group" mean(ps)
# should reject the null hypothesis
df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)
accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
accs = [x for x in accs if x != ref_acc]
Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]

# use linear regression on Oₐₚₒₛᵧₘs and Oₛᵧₘs
X = Oₛᵧₘs
y = Oₐₚₒₛᵧₘs
β = X\y 
@info "Linear Regression of Reference Group with Respect to KJ642617" β

# @show β * df.Nₛᵧₘ[1]/df.Nₐₚₒₛᵧₘ[1]
# @show mean(y - X*β)
# @show std(y - X*β)
# calculate R²
ŷ = X*β
ȳ = mean(y)
SST = sum((y .- ȳ).^2)
SSE = sum((ŷ - y).^2)
# @show 1-SSE/SST
@info "R² for the reference group" 1-SSE/SST

accs = accsₐ
Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]

# use linear regression on Oₐₚₒₛᵧₘs and Oₛᵧₘs
X = [ones(length(Oₛᵧₘs)) Oₛᵧₘs]
y = Oₐₚₒₛᵧₘs
β = X\y 
@info "Linear Regression of APOBEC3 Group with Respect to KJ642617" β
# @show β * df.Nₛᵧₘ[1]/df.Nₐₚₒₛᵧₘ[1]
# @show mean(y - X*β)
# @show std(y - X*β)
# calculate R²
ŷ = X*β
ȳ = mean(y)
SST = sum((y .- ȳ).^2)
SSE = sum((ŷ - y).^2)
# @show 1-SSE/SST
@info "R² for the reference group" 1-SSE/SST

# Test 4: test if MK783028 undergoes the same evolutionary environment as the reference group.
ref_acc = "KJ642617"
df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)

accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
accs = [x for x in accs if x != ref_acc]
accs = [accs; "MK783028"]
Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]
@info "p-value for MK783028 can be joined in reference group using $(ref_acc) as reference" pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs)

ref_acc = "KJ642617"
df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)

accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
accs = [x for x in accs if x != ref_acc]
accs = [accs; "MK783028"]
Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]
@info "p-value for MK783028 can be joined in reference group using $(ref_acc) as reference" pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs)

ref_acc = "MK783028"
df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)
accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
accs = [accs; "MK783028"]
accs = [x for x in accs if x != ref_acc]
Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]
@info "p-value for MK783028 can be joined in reference group using $(ref_acc) as reference" pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs)

accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
years = [2017,2018,2021,2022]
ps = []
for i in 1:100
    accs_rand = [random_acc_by_year(accsₐ, dates, year) for year in years]
    push!(ps,  pᵦ₊(accs_rand, ref_acc))
end
@info "p-value for MK783028 is the common ancestor of the APOBEC3 group" mean(ps) minimum(ps)



