# reference
# profiles/N100_s0.jl
# # standard case, no selection 
# # N ≈ 500
# g₀ = Genotype(
#     1, # id
#     Int(1e6),  # total synonymous mutations
#     Int(1e4),  # total APOBEC3 mutations
#     Int(1e4),  # total reverse APOBEC3 mutations
#     0,  # synonymous mutations
#     0,  # APOBEC3 mutations
#     0  # reverse APOBEC3 mutations
#    )
#    Par = Parameters(
#    1.0,  # β₀
#    0.5,  # μ₀
#    0.0,  # s
#    [1e-7, 1e-5, 1e-6],  # δ
#    t -> minimum([1e3 ,1e2 * (1+exp(0.01t))])  # K
# )
# check directory
if isdir("code")
    cd("code")
end
# this code generates the profiles for given N and s and writes into the profiles folder
# the profiles are generated by the function generate_profile

function generate_profile(N, s, sₕ = 0.0)
    return """
# N = $N
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
   $s,  # s
    $sₕ,  # sₕ
   [1e-7, 1e-5, 1e-6, 1e-4],  # δ
   1000.0, #T
   t -> minimum([$(10N) ,$(N) * (1+exp(0.005t))])  # K
)
"""
end

function constant_K_generate_profile(N, s,sₕ = 0.0)
    return """
# N = $N
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
   $s,  # s
    $sₕ,  # sₕ
   [1e-7, 1e-5, 1e-6, 1e-4],  # δ
   1000.0, #T
   t -> $(5N)  # K
)
"""
end

Ns = [200]
ss = [10.0^x for x in -5:1:1]
sₕs = [10.0^x for x in -7:0.5:2]
# for N in Ns
#     for s in ss
#         for sₕ in sₕs
#             open("profiles/extended_N$(N)_s$(s)_$(sₕ).jl", "w") do f
#                 write(f, generate_profile(N, s, sₕ))
#             end
#         end
#     end
# end

for N in Ns
    for s in ss
        for sₕ in sₕs
            open("profiles/extended_N$(N)_s$(s)_$(sₕ)_constant_K.jl", "w") do f
                write(f, constant_K_generate_profile(N, s, sₕ))
            end
        end
    end
end

for N in [10.0^x for x in 2:0.1:4]
    for s in ss
        sₕ = 0.0
            open("profiles/extended_N$(N)_s$(s)_$(sₕ)_constant_K.jl", "w") do f
                write(f, constant_K_generate_profile(N, s, sₕ))
            end
    end
end