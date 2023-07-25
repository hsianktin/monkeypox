using Statistics
∑(x) = sum(x)
using Statistics, Distributions, GLM
using LinearAlgebra
using DataFrames,CSV
using Dates 
using Pipe
function linear_regression(X, y)
    n, p = size(X)
    
    # Compute the least squares solution
    beta = X \ y
    
    # Compute the residuals
    residuals = y - X * beta
    
    # Compute the residual standard error
    rss = sum(residuals.^2) / (n - p)
    se = sqrt.(diag(rss * inv(X'X)))
    
    # compute R^2
    y_mean = mean(y)
    numerator = sum((y .- X * beta).^2)
    denominator = sum((y .- y_mean).^2)
    rse = 1 - numerator / denominator

    return beta, se, sqrt(rss)
end

function weighted_least_squares(x, y, weights)
    # calculate the weighted design matrix
    W = diagm(sqrt.(weights))
    Xw = W * x
    
    # calculate the weighted response vector
    yw = W * y
    
    # calculate the regression coefficients
    β = inv(Xw' * Xw) * (Xw' * yw)
    
    # calculate the residuals
    residuals = y - x * β
    
    # calculate the mean and variance of the residuals
    mean_residuals = mean(residuals)
    var_residuals = var(residuals)
    
    # calculate the weighted R^2
    y_mean = mean(y)
    numerator = sum(weights .* (y .- x * β).^2)
    denominator = sum(weights .* (y .- y_mean).^2)
    weighted_r2 = 1 - numerator / denominator
    
    # return the regression coefficients, residual statistics, and weighted R^2
    return β, weighted_r2
end

## FIG 1a
ref_acc = "KP849470"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
post2016_accs = [x for x in post2016_accs if x != ref_acc]

plot_df = @pipe df |> filter(row -> row.acc ∈ pre2016_accs, _)
β, se, rse = linear_regression(reshape(plot_df.Oₛᵧₘ, nrow(plot_df),1), plot_df.Oₐ₊₋)
# confidence interval
@show 1.96 .* se
# ci = β .+ [-1.96, 1.96] .* se

## FIG 1b
ref_acc = "MK783028"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
post2016_accs = [x for x in post2016_accs if x != ref_acc]

plot_df = @pipe df |> filter(row -> row.acc ∈ post2016_accs, _)
β, se, rse = linear_regression(reshape(plot_df.Oₛᵧₘ, nrow(plot_df),1), plot_df.Oₐ₊₋)
@show 1.96 .* se

## FIG 2a
ref_acc = "KJ642617"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
post2016_accs = [x for x in post2016_accs if x != ref_acc]

plot_df_post = @pipe df |> filter(row -> row.acc ∈ post2016_accs, _)
β_post, se_post, rse_post = linear_regression(hcat(ones(nrow(plot_df_post)), plot_df_post.Oₛᵧₘ), plot_df_post.Oₐ₊₋)


plot_df_pre = @pipe df |> filter(row -> row.acc ∈ pre2016_accs, _)
β_pre, se_pre, rse_pre = linear_regression(reshape(plot_df_pre.Oₛᵧₘ, nrow(plot_df_pre), 1), plot_df_pre.Oₐ₊₋)

# calculate the confidence interval
coords = []
for i in 1:10000
    β_pre_sample = β_pre .+ randn() .* se_pre
    β_post_sample = β_post .+ randn() .* se_post
    # calculate the intersection of the confidence intervals
    # where n_a = β_pre_sample[1] * n_s = β_post_sample[2] * n_s + β_post_sample[1]
    # solve for n_s
    n_s = β_post_sample[1] / (β_pre_sample[1] - β_post_sample[2])
    push!(coords, (n_s, β_pre_sample[1] * n_s))
end
n_s = β_post[1] / (β_pre[1] - β_post[2])
# compute the 95% confidence interval by taking the 2.5th and 97.5th percentiles
@show quantile([coord[1] for coord in coords], [0.025, 0.975])
# quantile([coord[1] for coord = coords], [0.025, 0.975]) = [11.610302717525197, 12.767675949924808]

## FIG 4b
ref_acc = "KJ642617"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
post2016_accs = [x for x in post2016_accs if x != ref_acc]

plot_df = @pipe df |> filter(row -> row.acc ∈ post2016_accs, _)
β, se, rse = linear_regression(hcat(ones(nrow(plot_df)), plot_df.Oₛᵧₘ), plot_df.Oₐ₊₋)

@show β[1]
# -9.269549618113123
@show 1.96 .* se[1]
# 0.33792660401062125
# compute the p-value for β₁ == 0 by t-test
@show 2 * (1 - cdf(Normal(0,1), abs(β[1] / se[1])))
# 0.0

## FIG 4c
ref_acc = "MK783029"
df = CSV.read("./data/extended_APOBEC_A2C_synonymous_$(ref_acc).csv", DataFrame)
pre2016_accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
pre2016_accs = [x for x in pre2016_accs if x != ref_acc]
post2016_accs = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
post2016_accs = [x for x in post2016_accs if x != ref_acc]


plot_df = @pipe df |> filter(row -> row.acc ∈ pre2016_accs, _)
β, se, rse = linear_regression(hcat(ones(nrow(plot_df)), plot_df.Oₛᵧₘ), plot_df.Oₐ₊₋)
@show β[1]
# 2.6628151420390633
@show 1.96 .* se[1]
# 0.5548464513510758
# compute the p-value for β₁ == 0 by t-test
@show 2 * (1 - cdf(Normal(0,1), abs(β[1] / se[1])))


plot_df = @pipe df |> filter(row -> row.acc ∈ post2016_accs, _)
β, se, rse = linear_regression(hcat(ones(nrow(plot_df)), plot_df.Oₛᵧₘ), plot_df.Oₐ₊₋)



## FIG 3b
ref_acc = "KJ642617"
df = CSV.read("./data/lineage_based_clock.csv", DataFrame)
x = df.date .- Date("2016-10-01")
# convert to year
x = x .|> (x -> (Dates.value.(x))/365.25)
df.time = x
X = [ones(length(x)) x]
y = df.Oₛᵧₘ
# weight by year 
years = unique(df.date .|> (x -> Dates.year(x)))
n_years = [sum(df.date .|> (x -> Dates.year(x)) .== year) for year in years]
weights = [ 1/(ξ) for ξ in x]
β, se, rse = linear_regression(X, y)
ci = confint(glm(@formula(Oₛᵧₘ ~ 1 + time), df, Normal(), wts=weights), 0.95)
β = mean(ci, dims=2)
se = [ci[1,2]- ci[1,1], ci[2,2]- ci[2,1]] ./ 1.96*2

yꜛ=12.67
yꜜ=11.61
# find time at which yꜛ is reached
t₀s = []
push!(t₀s, (yꜛ - ci[1,1]) / ci[2,1])
push!(t₀s, (yꜛ - ci[1,2]) / ci[2,2])
push!(t₀s, (yꜛ - ci[1,1]) / ci[2,2])
push!(t₀s, (yꜛ - ci[1,2]) / ci[2,1])
push!(t₀s, (yꜜ - ci[1,1]) / ci[2,1])
push!(t₀s, (yꜜ - ci[1,2]) / ci[2,2])
push!(t₀s, (yꜜ - ci[1,1]) / ci[2,2])
push!(t₀s, (yꜜ - ci[1,2]) / ci[2,1])
Δt = (minimum(t₀s)* 365.25 |> round |> Dates.Day, maximum(t₀s)* 365.25 |> round |> Dates.Day)
@show Δt .+ Date("2016-10-01")

# plot the 95% confidence interval of the lines
function yₘₐₓ(t, ci)
ys = []
push!(ys, ci[1,1] + ci[2,1] * t)
push!(ys, ci[1,2] + ci[2,2] * t)
push!(ys, ci[1,1] + ci[2,2] * t)
push!(ys, ci[1,2] + ci[2,1] * t)
return maximum(ys)
end

function yₘᵢₙ(t, ci)
ys = []
push!(ys, ci[1,1] + ci[2,1] * t)
push!(ys, ci[1,2] + ci[2,2] * t)
push!(ys, ci[1,1] + ci[2,2] * t)
push!(ys, ci[1,2] + ci[2,1] * t)
return minimum(ys)
end
δΔts = [Date("2015-12-26") - Date("2016-10-01"), Date("2023-01-01") - Date("2016-10-01")]
δΔts = δΔts .|> (x -> Dates.value(x) / 365.25)
Δts = collect(-1.8:0.1:6) ∪ δΔts

yₘₐₓs = [yₘₐₓ(Δt, ci) for Δt in Δts]
yₘᵢₙs = [yₘᵢₙ(Δt, ci) for Δt in Δts]
t = [Dates.Day(floor(Δt * 365.25)) + Date("2016-10-01") for Δt in Δts]
plot_df = DataFrame(t=t, yₘₐₓs=yₘₐₓs, yₘᵢₙs=yₘᵢₙs)
# filter plot_df such that Date("2016-01-01") < t < Date("2023-01-01")
plot_df = plot_df[plot_df.t .≥ Date("2015-12-01"), :]
plot_df = plot_df[plot_df.t .≤ Date("2023-01-01"), :]
# sort by t
plot_df = sort(plot_df, :t)
CSV.write("./data/lineage_based_clock_CI.csv", plot_df)


## FIG 3a
ref_acc = "KJ642617"
df = CSV.read("./data/molecular_clock.csv", DataFrame)
x = df.date .- Date("2016-10-01")
# convert to year
x = x .|> (x -> (Dates.value.(x))/365.25)
df.time = x
######

X = [ones(length(x)) x]
y = df.Oₛᵧₘ
# weight by year 
years = unique(df.date .|> (x -> Dates.year(x)))
n_years = [sum(df.date .|> (x -> Dates.year(x)) .== year) for year in years]
weights = [ 1/(ξ) for ξ in x]
β, se, rse = linear_regression(X, y)
ci = confint(glm(@formula(Oₛᵧₘ ~ 1 + time), df, Normal(), wts=weights), 0.95)
β = mean(ci, dims=2)
se = [ci[1,2]- ci[1,1], ci[2,2]- ci[2,1]] ./ 1.96*2
yꜛ=12.67
yꜜ=11.61
# find time at which yꜛ is reached
t₀s = []
push!(t₀s, (yꜛ - ci[1,1]) / ci[2,1])
push!(t₀s, (yꜛ - ci[1,2]) / ci[2,2])
push!(t₀s, (yꜛ - ci[1,1]) / ci[2,2])
push!(t₀s, (yꜛ - ci[1,2]) / ci[2,1])
push!(t₀s, (yꜜ - ci[1,1]) / ci[2,1])
push!(t₀s, (yꜜ - ci[1,2]) / ci[2,2])
push!(t₀s, (yꜜ - ci[1,1]) / ci[2,2])
push!(t₀s, (yꜜ - ci[1,2]) / ci[2,1])
Δt = (minimum(t₀s)* 365.25 |> round |> Dates.Day, maximum(t₀s)* 365.25 |> round |> Dates.Day)
@show Δt .+ Date("2016-10-01")

# plot the 95% confidence interval of the lines
function yₘₐₓ(t, ci)
    ys = []
    push!(ys, ci[1,1] + ci[2,1] * t)
    push!(ys, ci[1,2] + ci[2,2] * t)
    push!(ys, ci[1,1] + ci[2,2] * t)
    push!(ys, ci[1,2] + ci[2,1] * t)
    return maximum(ys)
end

function yₘᵢₙ(t, ci)
    ys = []
    push!(ys, ci[1,1] + ci[2,1] * t)
    push!(ys, ci[1,2] + ci[2,2] * t)
    push!(ys, ci[1,1] + ci[2,2] * t)
    push!(ys, ci[1,2] + ci[2,1] * t)
    return minimum(ys)
end

δΔts = [Date("2016-01-01") - Date("2016-10-01"), Date("2023-01-01") - Date("2016-10-01")]
δΔts = δΔts .|> (x -> Dates.value(x) / 365.25)
Δts = collect(-7.8:0.1:1) ∪ δΔts
yₘₐₓs = [yₘₐₓ(Δt, ci) for Δt in Δts]
yₘᵢₙs = [yₘᵢₙ(Δt, ci) for Δt in Δts]
ŷs = [β[1] + β[2] * Δt for Δt in Δts]
t = [Dates.Day(floor(Δt * 365.25)) + Date("2016-10-01") for Δt in Δts]
plot_df = DataFrame(t=t, yₘₐₓs=yₘₐₓs, yₘᵢₙs=yₘᵢₙs, ŷs=ŷs)
# filter plot_df such that Date("2016-01-01") < t < Date("2023-01-01")
plot_df = plot_df[plot_df.t .≥ Date("2015-12-01"), :]
plot_df = plot_df[plot_df.t .≤ Date("2023-01-01"), :]
# sort by t
plot_df = sort(plot_df, :t)
CSV.write("./data/molecular_clock_CI.csv", plot_df)


### Fig Resp 2
df = CSV.read("./data/model_data.csv", DataFrame)
function pₕ(Oₐₚₒₛᵧₘs, Oₛᵧₘs)
    Pₐₚₒ = ∑(Oₐₚₒₛᵧₘs)/∑(Oₛᵧₘs)
    ps = [] # p-value
    for (Oₐₚₒₛᵧₘ, Oₛᵧₘ) ∈ zip(Oₐₚₒₛᵧₘs, Oₛᵧₘs)
        Ôₐₚₒₛᵧₘ = Oₛᵧₘ * Pₐₚₒ
        ∑([binom(Oₛᵧₘ, k) * Pₐₚₒ^k * (1 - Pₐₚₒ)^(Oₛᵧₘ-k) for k in 0:1:Oₛᵧₘ if abs(k - Ôₐₚₒₛᵧₘ) < abs(Oₐₚₒₛᵧₘ - Ôₐₚₒₛᵧₘ) ])
        push!(ps, maximum([
            1-∑(    
                [binom(Oₛᵧₘ, k) * Pₐₚₒ^k * (1 - Pₐₚₒ)^(Oₛᵧₘ-k) for k in 0:1:Oₛᵧₘ if abs(k - Ôₐₚₒₛᵧₘ) < abs(Oₐₚₒₛᵧₘ - Ôₐₚₒₛᵧₘ) ]
            ), 0]) )
    end
    @show ps
    # Fisher's method for combining p-values
    function fishers_method(p_values)
        k = length(p_values)
        test_statistic = -2 * sum(log.(p_values))
        combined_p_value = ccdf(Chisq(2k), test_statistic)
        return combined_p_value
    end
    
    return fishers_method(ps)
end
Oₐₚₒₛᵧₘs = df.Oₐₚₒₛᵧₘ
Oₛᵧₘs = df.Oₛᵧₘ
pₕ(df.Oₐ₊₋, df.Oₛᵧₘ)

select_df = df[2:end,:]
pₕ(select_df.Oₐ₊₋, select_df.Oₛᵧₘ)
Oₐₚₒₛᵧₘs = select_df.Oₐₚₒₛᵧₘ
Oₛᵧₘs = select_df.Oₛᵧₘ

xs = 1:1:350
function prediction_interval(Oₐₚₒₛᵧₘs, Oₛᵧₘs, new_Oₛᵧₘs)
    Pₐₚₒ = ∑(Oₐₚₒₛᵧₘs)/∑(Oₛᵧₘs)
    prob(x,lower_bound, upper_bound) = ∑([binom(x, k) * Pₐₚₒ^k * (1 - Pₐₚₒ)^(x-k) for k in lower_bound:1:upper_bound  ])
    upper_bounds = []
    lower_bounds = []
    for x in new_Oₛᵧₘs
        if x < 0
            error("new_Oₛᵧₘs should be positive")
        end
        @show x
        @show ȳ = Pₐₚₒ * x 
        upper_bound = minimum([ceil(Int,ȳ + 1), x])
        lower_bound = maximum([floor(Int,ȳ - 1), 0])
        @show p = prob(x,lower_bound , upper_bound)
        while upper_bound < x && lower_bound > 0 && prob(x,lower_bound , upper_bound) < 0.95
            @show prob(x,lower_bound , upper_bound)
            upper_bound += 1
            lower_bound -= 1
        end
        push!(upper_bounds, upper_bound)
        push!(lower_bounds, lower_bound)
    end
    return upper_bounds, lower_bounds
end
x = 341
ys = prediction_interval(select_df.Oₐ₊₋, select_df.Oₛᵧₘ, xs)

# save xs, ys[1], ys[2] to csv
pi_df = DataFrame(x=xs, upper=ys[1], lower=ys[2])
CSV.write("./data/model_prediction_interval.csv", pi_df)

using PGFPlotsX
p = @pgf Axis(
    {
        xlabel = "number of synonymous mutations \$n \$",
        ylabel = "number of APOBEC3-relevant mutations \$n_{\\rm a}\$",
        width = "10cm",
        height = "10cm",
        # show grid lines
        grid = "major",
        # set legend pos to north west
        legend_pos = "north west",
    },
    Plot(
        {
            only_marks,
            # set mark size to 1.5pt
            # mark_size = "1pt",
            # let the mark to be semi transparent
            mark_options = {fill_opacity = 0.1},
            color = "black"
        },
        Table(x=df.Oₛᵧₘ, y=df.Oₐ₊₋)
    ),
    Plot(
        {
            no_marks,
            red!30,
        },
        Table(x=xs, y=ys[1])
    )
    ,
    Plot(
        {
            no_marks,
            red!30,
        },
        Table(x=xs, y=ys[2])
    ),
) 

