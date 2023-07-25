using Distributed 
addprocs(16)
@everywhere using ProgressMeter

@everywhere function run_profile(profile)
    (run∘pipeline)(`julia extended_dna_simulation.jl $profile`)
end

if length(ARGS) == 0
    @error "No profile is specified!"
else
    # 1. if the arguments does not contain wild card, then use the arguments as the profiles
    # 2. if the arguments contain wild card, then use the arguments as the pattern to match the profiles
    if !occursin("*", ARGS[1])
        profiles = ARGS
    else
        profiles = filter(x -> occursin(ARGS[1], x), readdir("profiles"))
    end
end


@info "The following profiles will be executed:"

using Random
try
    # for profile in profiles
    #     @show profile
    #     run_profile(profile)
    # end
    # permute the profiles
    Random.shuffle!(profiles)
    progress_pmap(run_profile, profiles; progress = Progress(length(profiles), color=:white, showspeed = true,barglyphs=BarGlyphs("[=> ]")))
    # @time pmap(run_profile, profiles)
    # the content of mail_success is 
    mail_body = """
    From: "Automata" <automata@ryzen>
    To: "Xiangting Li" <xiangting.li@ucla.edu>
    Subject: monkeypox simulation task completed

    Hi Xiangting,

    The tasks scheduled in monkeypox simulation are completed.
    The profiles are:
        - $(join(profiles, "\n    - "))

    Please check the result.
    Automata.
    """
    open("mail_success.txt", "w") do f
        write(f, mail_body)
    end
    run(pipeline(`ssmtp xiangting.li@ucla.edu`; stdin = "mail_success.txt"))
    # remove the mail_success
    rm("mail_success.txt")
catch e
    mail_body = """
    From: "Automata" <automata@ryzen>
    To: "Xiangting Li" <xiangting.li@ucla.edu>
    Subject: monkeypox simulation task failed

    Hi Xiangting,

    The tasks scheduled in monkeypox simulation encountered a problem. Please review the logs and fix the problem.
    $(e)
    The profiles are:
        - $(join(profiles, "\n     -"))
    
    Automata.
    """
    open("mail_failure.txt", "w") do f
        write(f, mail_body)
    end
    run(pipeline(`ssmtp xiangting.li@ucla.edu`; stdin = "mail_failure.txt"))
    # remove the mail_failure
    rm("mail_failure.txt")
end

