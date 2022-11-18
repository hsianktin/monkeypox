# Introduction
The repository to extract the number of possible and observed
synonymous, APOBEC3-relevant and -irrelevant mutations.

Due to recent overhaul of API in `BioSequences`, please
make sure the package version is 3.X by using
```julia
pkg>add BioSequences@3
```



# Workflow
Data preprocessing based on accession numbers.
```
accession.txt |> access.jl  |> data, metadata.csv
```

## Alignment
Run `identify_genes.jl` to annotate the genes based on the reference genome.

Run `embedding.jl` to generate pairwise alignment and data for further analysis.

## Statistical tests
Run `statistical_test.jl` to generate the results of statistical tests.

Run `comparison.jl` to identify the equivalent classes of genomes and generate a distance matrix of the equivalent classes.

Run `molecular_clock.jl` to generate the results of molecular clock analysis.
