using CSV, DataFrames, Pipe
include("grouping.jl")

ref_acc = "KJ642617"
df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)
# calculate O₀ = Oₛᵧₘ - Oₐ₊₋
df.O₀ = df.Oₛᵧₘ .- df.Oₐ₊₋
CSV.write("./data/APOBEC_synonymous_$(ref_acc).csv", df)

# accs = [accsₐ; ref_acc]
accs = accsₐ
# linear regression on O₀ and df.year 
y = [df.O₀[df.acc .== acc][1] for acc in accs]
x = [df.year[df.acc .== acc][1] for acc in accs]
X = [ones(length(x)) x]
@show β = X\y