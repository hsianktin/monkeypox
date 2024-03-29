# N = 1000
g₀ = Genotype(
    1, # id
    Int(1e6),  # total synonymous mutations
    Int(1e4),  # total APOBEC3 mutations
    Int(1e4),  # total reverse APOBEC3 mutations
    0,  # synonymous mutations
    0,  # APOBEC3 mutations
    0  # reverse APOBEC3 mutations
   )
   Par = Parameters(
   1.0,  # β₀
   0.5,  # μ₀
   0.001,  # s
   [1e-7, 1e-5, 1e-6],  # δ
   t -> 5000  # K
)
