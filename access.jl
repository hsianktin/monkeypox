using CSV
using Pipe
using DataFrames
using Dates
using ProgressMeter
using BioSequences
using BioFetch
using FASTX
using GenomicAnnotations

acc_df = CSV.read("./data/accession.txt", DataFrame)

# download using BioFetch
@showprogress for acc in acc_df.Accession
    @info "download $(acc)"
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
end

df = @pipe DataFrame(Accession = acc_df.Accession, Date = (Date∘Year).(acc_df.Year)) |>
    sort(_, :Date) 

# save the dataframe to a csv file
CSV.write("metadata.csv", df)