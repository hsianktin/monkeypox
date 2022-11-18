using DataFrames,CSV,Pipe
# utils to divide the entries into two groups
# the APOBEC group and the reference group
df = CSV.read("metadata.csv",DataFrame)
accsₐ = @pipe df |> filter(row ->row.Date > Date(2010,1,1),_) |> (x -> x.Accession) 
dates = @pipe df |> filter(row ->row.Date > Date(2010,1,1),_) |> (x -> x.Date)
# exclude NC_063383
# dates = dates[accsₐ .!= "NC_063383"]
# accsₐ = accsₐ[accsₐ .!= "NC_063383"]
accsᶜ = df.Accession[[df.Accession[i] ∉ accsₐ for i in 1:length(df.Accession)]]

# utils for sampling by year
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
