using CSV
using DataFrames
df = CSV.read("meta_data_doi.csv", DataFrame)
Reference = [
    "\\cite{$(bib_cite_entry)}" for bib_cite_entry in df.Bib_cite_entry
]

# `_` is recognized as a special character in LaTeX, so we need to replace it with `\_`
df.Accession = [
    replace(accession, "_" => "\\_") for accession in df.Accession
]

df.Reference = Reference

tex_df = df[:, [:Accession, :Date, :Clade, :Reference]]
CSV.write("meta_data.tex", tex_df)