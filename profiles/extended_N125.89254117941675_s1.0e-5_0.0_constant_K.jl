# N = 125.89254117941675
g₀ = Genotype(
    1, # id
    Int(1e6),  # total synonymous mutations
    Int(1e4),  # total APOBEC3 mutations
    Int(1e4),  # total reverse APOBEC3 mutations
    Int(1e3),  # total hidden mutations
    0,  # synonymous mutations
    0,  # APOBEC3 mutations
    0,  # reverse APOBEC3 mutations
    0  # hidden mutations
   )
   Par = Parameters(
   1.0,  # β₀
   0.5,  # μ₀
   1.0e-5,  # s
    0.0,  # sₕ
   [1e-7, 1e-5, 1e-6, 1e-4],  # δ
   1000.0, #T
   t -> 125.89254117941675  # K
)
