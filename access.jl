using CSV
using Pipe
using DataFrames
using Dates
using ProgressMeter
using BioSequences
using BioFetch
using FASTX
using GenomicAnnotations


@info "starting to download data"
acc_df = CSV.read("collection_dates.csv", DataFrame)
acc_string = String[]
# download using BioFetch

progress = Progress(length(acc_df.Accession), 1)
for acc in acc_df.Accession
    s = "downloading $(acc)"
    # gbseq = fetchseq(acc, format=gb)[1]
    push!(acc_string, acc)
    # download the fasta file
    if !isfile("data/$(acc).fasta")
        # fetch fasta and write to file 
        fastaseq = fetchseq(acc, format=fasta)[1]
        open(FASTA.Writer, "data/$(acc).fasta") do io
            write(io, fastaseq)
        end
        # sleep(0.33)
    end
    # download the gbk file
    if !isfile("data/$(acc).gbk")
        gbseq = fetchseq(acc, format=gb)[1]
        open(GenBank.Writer, "data/$(acc).gbk") do io
            write(io, gbseq)
        end
        # sleep(0.33)
    end
    # avoid being blocked by NCBI
    next!(progress; showvalues = [(:message, s)])
end

df = @pipe DataFrame(Accession = acc_string, Date = acc_df.Date) |>
    sort(_, :Date)

# save the dataframe to a csv file
CSV.write("meta_data.csv", df)
@info "Done"