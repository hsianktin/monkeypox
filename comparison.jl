using DataFrames
using CSV
using Statistics
using Pipe
using Dates
∑(x) = sum(x)
include("grouping.jl")

accs = accsₐ

function D₊(acc1,acc2)
    temp_df = CSV.read("APOBEC_synonymous_$(acc2).csv", DataFrame)
    # find the row where row.acc == acc1 and return the value pₐₚₒₛᵧₘ
    Oₐₚₒₛᵧₘ = temp_df[temp_df.acc .== acc1, :pₐₚₒₛᵧₘ][1]
    return Oₐₚₒₛᵧₘ
end
acc1 = "KJ642617"
acc2 = "MK783028"

D₊(acc1,acc2)
D₊(acc2,acc1)


function D₋(acc1,acc2)
    temp_df = CSV.read("APOBEC_synonymous_$(acc2).csv", DataFrame)
    # find the row where row.acc == acc1 and return the value pₐₚₒₛᵧₘ
    Oₐₚₒₛᵧₘ⁻¹ = temp_df[temp_df.acc .== acc1, :Oₐₚₒₛᵧₘ⁻¹][1]
    return Oₐₚₒₛᵧₘ⁻¹
end


D₊(accs[10],accs[3])

# get the distance matrix
M₊ = [D₊(accs[i],accs[j]) for i in 1:length(accs), j in 1:length(accs)]
# M₋ = [D₋(accs[i],accs[j]) for i in 1:length(accs), j in 1:length(accs)]

# for robustness, first identify equivalent sets.
# acc1 and acc2 are equivalent if d(acc1,acc2) = d(acc2,acc1) = 0

eq_accs = []

for i in 1:length(accs), j in 1:length(accs)
    if i < j
        if M₊[i,j] == M₊[j,i] == 0
            push!(eq_accs, [accs[i], accs[j]])
        end
    end
end

# identify equivalent classes
eq_classes = []
for i in 1:length(eq_accs)
    acc1, acc2 = eq_accs[i]
    if length(eq_classes) == 0
        push!(eq_classes, Set([acc1, acc2]))
    else
        append_flag = true
        for j in 1:length(eq_classes)
            if acc1 ∈ eq_classes[j]
                push!(eq_classes[j], acc2)
                append_flag = false
            elseif acc2 ∈ eq_classes[j]
                push!(eq_classes[j], acc1)
                append_flag = false
            end
            # println("i,j=($i,$j)")
        end
        if append_flag
            push!(eq_classes, Set([acc1, acc2]))
        end
    end
end


for acc in accs
    in_eq_class_flag = false
    for eq_class ∈ eq_classes
        if acc ∈ eq_class 
            in_eq_class_flag = true
        end
    end
    if !in_eq_class_flag
        push!(eq_classes,Set([acc]))
    end
end

@info "equivalent classes" eq_classes

# now we have the equivalent classes, we can calculate the distance matrix
# for each class, we take the average of the distance matrix

distance_matrix = zeros(length(eq_classes), length(eq_classes))
for i in 1:length(eq_classes)
    for j in 1:length(eq_classes)
        if i ≠ j
            distance_matrix[i,j] = mean(
                [M₊[findfirst(isequal(acc1), accs), findfirst(isequal(acc2), accs)] for acc1 in eq_classes[i], acc2 in eq_classes[j]]
                )
        end
    end
end
@info "distance matrix" distance_matrix
# write the distance matrix to a csv file

# define i → j if distance_matrix[j,i] > 0, distance_matrix[i,j] = 0
# find all maximal chains of the form i1→i2→i3→...→in
# a maximal chain is a chain that is not a subchain of any other chain
function trees(distance_matrix)
    n = size(distance_matrix,1)
    i→j = if distance_matrix[j,i] > 0 && distance_matrix[i,j] == 0
        true
    else
        false
    end
    # the value of i → j is true when distance_matrix[j,i] > 0, false otherwise

    # first, start with i, find all j such that i → j
    # then, for each j, find all k such that j → k
    # repeat until no more j can be found

    # find all maximal chains
    # maximal_chains = []
    function extend_chain(chain, set)
        # chain is a list of indices such that chain[i] → chain[i+1]
        # set is the set of indices that have not been used yet
        # find all j such that chain[end] → j
        # return a set of chains
        # @show chain
        # @show set
        chains = []
        if length(set) == 0
            return [chain]
        else
            extensible_flag = false
            for j in 1:length(set)
                # @show chain[end]→set[j]
                # @show set[j] ∉ chain
                if (chain[end]→set[j]) && (set[j] ∉ chain)
                    extensible_flag = true
                    new_chain = copy(chain)
                    push!(new_chain, set[j])
                    # cat the new chain to the existing chains
                    # @show chains
                    # @info "extended chain" new_chain
                    chains = vcat(chains, extend_chain(new_chain, set[set .!= set[j]]))
                    # @info "intermediate chains" chains
                end
            end
            if !extensible_flag
                # @info "terminal chain" chain
                return [chain]
            else
                return chains
            end
            # return chains
        end
    end

    chains = [extend_chain([i], collect(1:n)) for i in 1:n]
    chains = vcat(chains...)
    # remove the chain that is a subchain of another chain
    # for each chain, check if it is a subchain of any other chain
    # if it is, remove it

    ∑(A) = if length(A) == 0
        0
    else
        sum(A)
    end

    A ⊂ B = if ∑([1 for i in 1:length(A) if A[i] ∉ B]) == 0
        true
    else
        false
    end

    function issubchain(a,b)
        return a ⊂ b
    end

    # chain₀ ⊂ chain₁

    maximal_chains = []
    for chain₀ in chains
        subchain_flag = false 
        for chain₁ in chains
            # @show chain₀
            # @show chain₁
            # @show issubchain(chain₀, chain₁)
            if chain₀ ≠ chain₁
                if issubchain(chain₀, chain₁)
                    subchain_flag = true
                end
            end
        end
        if !subchain_flag
            push!(maximal_chains, chain₀)
        end
    end
    return maximal_chains
end

maximal_chains =  trees(distance_matrix)

@info "maximal chains" maximal_chains
@show roots = unique([maximal_chain[1] for maximal_chain in maximal_chains])
# molecular clock estimation based on the roots
root = roots[1]

eq_classes_dates = []
for eq_class in eq_classes
    eq_class_date = dates[[acc in eq_class for acc in accs]]
    push!(eq_classes_dates, eq_class_date)
end
eq_class_mean_dates = []
for eq_class_dates in eq_classes_dates
    eq_class_mean_date = @pipe eq_class_dates |> Dates.year.(_) |> minimum
    push!(eq_class_mean_dates,eq_class_mean_date)
end
eq_class_mean_dates
# combine maximal chains and equivalent classes into a newick tree
using Phylo
# simpletree = parsenewick("(human_root)animal_root;")

# ref_acc = "KP849470"
# df = CSV.read("./data/APOBEC_synonymous_$(ref_acc).csv", DataFrame)
# pₐₚₒₛᵧₘs = [df.pₐ₊₋[df.acc .== acc][1] for acc in accs]
# pₛᵧₘs = [df.pₛᵧₘ[df.acc .== acc][1] for acc in accs]
# # @show β₀, β₁, R² = affine(accsₐ, ref_acc)
# # @show β′, R²′ = linear(accsᶜ, ref_acc)

# # find intersection of y = β₀ + β₁x and y = β′₁x
# x₀ = (β₀)/(β′-β₁)
# y₀ = x₀ * β′
# τs = pₐₚₒₛᵧₘs .- y₀
# eq_classes_taus = []
# eq_classes_mean_taus = []
# for eq_class in eq_classes
#     eq_class_taus = τs[[acc in eq_class for acc in accs]]
#     push!(eq_classes_taus, eq_class_taus)
#     eq_class_mean_sym = @pipe eq_class_taus |> mean
#     push!(eq_classes_mean_taus, eq_class_mean_sym)
# end
# βs = []
# using Plots
# plot(legend = false)
# for i in 1:length(maximal_chains)
#     chain = maximal_chains[i]
#     if chain[1] == roots[1] && length(chain)> 1
#         plot_df = DataFrame(
#             acc = [],
#             date = [],
#             t = Float32[],
#             year = Int[],
#         )
#         for node in chain 
#             plot_accs = eq_classes[node] |> collect
#             plot_dates = [dates[accs .== acc][1] for acc in plot_accs]
#             plot_taus = [τs[accs .== acc][1] for acc in plot_accs]
#             plot_year = Dates.year.(plot_dates)
#             plot_df = vcat(plot_df, DataFrame(
#                 acc = plot_accs,
#                 date = plot_dates,
#                 t = plot_taus,
#                 year = plot_year,
#             ))
#         end
#         CSV.write("chain_clock_$(i).csv", plot_df)
#                 # get the dates and taus of the chain
#         chain_dates = [eq_class_mean_dates[node] for node in chain]
#         chain_taus = [eq_classes_mean_taus[node] for node in chain]    
#         # linear regression 
#         plot!(plot_df.t, plot_df.year, mark=:circle, label="chain $(i)")
#         X = [ones(length(plot_df.t)) plot_df.t]
#         y = plot_df.year
#         @show X
#         @show y
#         β = X\y
#         X = [1 0;X]
#         plot!(X[:,2], X*β)
#         CSV.write("chain_clock_linear_$(i).csv", DataFrame(
#             t = X[:,2],
#             date = X*β,
#         ))
#         push!(βs, β)
#     end
# end
# xaxis!("evolutionary time")
# yaxis!("collection date")

# @show βs 

# function d(i,j)
#     return (distance_matrix[i,1] + distance_matrix[j,1] - distance_matrix[i,j] - distance_matrix[j,i])/2
# end

# D = [d(i,j) for i in 1:length(eq_classes), j in 1:length(eq_classes)]

# use linear regression on Oₐₚₒₛᵧₘs and Oₛᵧₘs
# Xh = [ones(length(Oₐₚₒₛᵧₘs)) Oₛᵧₘs]
# yh = Oₐₚₒₛᵧₘs
# βh = Xh\yh


# accs = @pipe df |> filter(row ->row.date < Date(2010,1,1),_) |> (x -> x.acc)
# accs = [x for x in accs if x != ref_acc]
# Oₐₚₒₛᵧₘs = [df.Oₐ₊₋[df.acc .== acc][1] for acc in accs]
# Oₛᵧₘs = [df.Oₛᵧₘ[df.acc .== acc][1] for acc in accs]

# Xa = Oₛᵧₘs
# ya = Oₐₚₒₛᵧₘs
# βa = Xa\ya

# # compute the intersection y = βa x = βh[1] + βh[2] x
# x = ( - βh[1])/(βh[2]- βa)

# # evolutionary time 
# τs = yh .- x*βa

# using Plots
# acch = @pipe df |> filter(row ->row.date > Date(2010,1,1),_) |> (x -> x.acc)
# acch = [x for x in acch if x != exclude_acc]

# # scatter(τs, yh, label = "human")
# # eq_classes[1]
# # # plot acch which belong to eq_classes[1]
# # scatter( ones(length(τs[[acch[i] ∈ eq_classes[1] for i in 1:length(acch)]])),
# # τs[[acch[i] ∈ eq_classes[1] for i in 1:length(acch)]],
# # color = "red",
# #  label = "human eq_class 1")

# # # plot acch which belong to eq_classes[2]
# # scatter!( ones(length(τs[[acch[i] ∈ eq_classes[2] for i in 1:length(acch)]])).*2,
# # τs[[acch[i] ∈ eq_classes[2] for i in 1:length(acch)]],
# # color = "blue",
# #  label = "human eq_class 2")

# # # plot acch which belong to eq_classes[3]
# # scatter!( ones(length(τs[[acch[i] ∈ eq_classes[3] for i in 1:length(acch)]])).*3,
# # τs[[acch[i] ∈ eq_classes[3] for i in 1:length(acch)]],
# # color = "green",
# #  label = "human eq_class 3")

# # # plot acch which belong to eq_classes[4]
# # scatter!( ones(length(τs[[acch[i] ∈ eq_classes[4] for i in 1:length(acch)]])).*4,
# # τs[[acch[i] ∈ eq_classes[4] for i in 1:length(acch)]],
# # color = "yellow",
# #  label = "human eq_class 4")

# # # plot acch which belong to eq_classes[5]
# # scatter!( ones(length(τs[[acch[i] ∈ eq_classes[5] for i in 1:length(acch)]])).*5,
# # τs[[acch[i] ∈ eq_classes[5] for i in 1:length(acch)]],
# # color = "black",
# #  label = "human eq_class 5")

# # plot acch which belong to eq_classes[6]
# # scatter!( ones(length(τs[[acch[i] ∈ eq_classes[6] for i in 1:length(acch)]])).*6,
# # τs[[acch[i] ∈ eq_classes[6] for i in 1:length(acch)]],
# # color = "orange",
# #  label = "human eq_class 6")

# # draw arrows between eq_classes
# # according to maximal_chains
# # eq_classes[1] -> eq_classes[2]
# # eq_classes[1] -> eq_classes[3] -> eq_classes[4]
# # eq_classes[1] -> eq_classes[3] -> eq_classes[5] -> eq_classes[6]

# x1 = 0
# x2 = x1 - distance_matrix[2,1]
# x3 = x1 + distance_matrix[3,1]
# x4 = x3 - distance_matrix[4,3]
# distance_matrix[3,1]
# x5 = x3 + distance_matrix[5,3]
# distance_matrix[5,1]
# # x6 = x5 + distance_matrix[6,5]
# # distance_matrix[6,1]

# x = distance_matrix[:,1]
# x = [x1 x2 x3 x4 x5]
# X = []
# Y = []
# for i in 1:length(τs)
#     # find the eq_class of acch[i]
#     @show i
#     @show   j = findfirst(x -> acch[i] ∈ eq_classes[x] , collect(1:length(eq_classes)))
#     if !isnothing(j)
#         push!(X, x[j])
#         push!(Y, τs[i])
#     end
# end
# scatter(X, Y, label = "human")




# # for each eq_class
# # associate the mean of the evolutionary time
# # of the accs in the eq_class
# t = []
# for i in 1:length(eq_classes)
#     push!(t, minimum(τs[[acch[j] ∈ eq_classes[i] for j in 1:length(acch)]]))
# end

# # for each entry in eq_classes, create time_classes
# time_classes = []
# for i in 1:length(eq_classes)
#     accsᵢ = eq_classes[i]
#     timeclass = []
#     for acc in accsᵢ
#         pos = findfirst(x -> x == acc, acch)
#         push!(timeclass, τs[pos])
#     end
#     push!(time_classes, timeclass)
# end

# Δt = [t[1], t[2] - t[1], t[3] - t[1], t[4] - t[3], t[5] - t[4]]
# # simpletree = parsenewick("(((((class6:$(t[6]-t[5]))class5:$(t[5]-t[3]),class4:$(t[4]-t[3]))class3:$(t[3]-t[1]),class2:$(t[2]-t[1]))class1:$(t[1]),class11:$(t[1]))human_root:$(t[1]))animal_root;")
# # get the complete tree by adding the names provided by eq_classes[i]
# # obtain the newick string
# newick = "animal_root;"
# function get_layer_newick(eq_classes, i)
#     layer_newick = ""
#     eq_class= collect(eq_classes[i])
#     for j in 1:length(eq_class)
#         if j == 1
#             layer_newick = "$(eq_class[j]):$(time_classes[i][j]-t[i])"
#         else
#             layer_newick ="$(eq_class[j]):$(time_classes[i][j]-t[i]),$(layer_newick)"
#         end
#     end
#     layer_newick = "($(layer_newick))eq_class$(i):0.1"
#     return layer_newick
# end
# newick_string = "animal_root;"
# newick_string = "(((((((($(get_layer_newick(eq_classes,5)))class5:$(Δt[5]))$(get_layer_newick(eq_classes,4)))class4:$(Δt[4]),$(get_layer_newick(eq_classes, 3)))class3:$(Δt[3]),$(get_layer_newick(eq_classes,2)))class2:$(Δt[2]),$(get_layer_newick(eq_classes, 1)))class1:$(Δt[1]))$(newick_string)"
# # get_layer_newick(eq_classes, 2)
# complex_tree= parsenewick(newick_string)
# default(linecolor = :black, size = (400, 400)) # looks nicer with black lines
# plot(complex_tree,
# markersize = 5, markercolor = :steelblue, markerstrokecolor = :white,)

# savefig("time_to_acch.pdf")

# # write newick_string to file 
# open("time_to_acch.newick", "w") do io
#     write(io, newick_string)
# end
