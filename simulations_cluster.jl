using Distributed

using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using Dates
using DelimitedFiles

using ClusterManagers
ENV["JULIA_WORKER_TIMEOUT"] = "600"

addprocs(ClusterManagers.SlurmManager(250), N=8, topology=:master_worker, exeflags="--project=Project.toml"; W="300")
#addprocs(10)
@everywhere using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@everywhere include("abm_behavior.jl")
@everywhere const cv=abmbehavior


function run(myp::cv.ModelParameters, nsims=1000, folderprefix="./")
    
    println("starting $nsims simulations...\nsave folder set to $(folderprefix)")
    dump(myp)
        
    cdr = try
        pmap(1:nsims) do x                 
            result = cv.runsim(x, myp)
            GC.gc()
            return result
        end  
    catch e
        open("erro.txt", "a") do io
            println(io, "==== ERRO no pmap ====")
            println(io, "Data/hora: ", Dates.now())
            println(io, "Tipo: ", typeof(e))
            println(io, "Erro: ", e)
            println(io, "Stacktrace:")
            showerror(io, e, catch_backtrace())
            println(io, "\n")
        end
        nothing
    end
    
    if isnothing(cdr)
        error("Processo morto")
    end

    println("simulations finished")
    println("total size of simulation dataframes: $(Base.summarysize(cdr))")
    ## write the infectors 
    
    ## write contact numbers
    #writedlm("$(folderprefix)/ctnumbers.dat", [cdr[i].ct_numbers for i = 1:nsims])    
    ## stack the sims together
    allag = vcat([cdr[i].a  for i = 1:nsims]...)
    allv = vcat([cdr[i].allv  for i = 1:nsims]...)
    # ag1 = vcat([cdr[i].g1 for i = 1:nsims]...)
    # ag2 = vcat([cdr[i].g2 for i = 1:nsims]...)
    # ag3 = vcat([cdr[i].g3 for i = 1:nsims]...)
    # ag4 = vcat([cdr[i].g4 for i = 1:nsims]...)
    # ag5 = vcat([cdr[i].g5 for i = 1:nsims]...)
    # ag6 = vcat([cdr[i].g5 for i = 1:nsims]...)
    # ag7 = vcat([cdr[i].g5 for i = 1:nsims]...)
    mydfs = Dict("all" => allag, "allv" => allv)#, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "ag4" => ag4, "ag5" => ag5, "ag6" => ag6, "ag7" => ag7)
    #mydfs = Dict("all" => allag)
    
    ## save at the simulation and time level
    ## to ignore for now: miso, iiso, mild 
    c1 = Symbol.((:LAT, :INF, :DED), :_INC)
    #c2 = Symbol.((:LAT, :INF, :IISO, :HOS, :ICU, :DED), :_PREV)
    #c1 = Symbol.((:LAT, :HOS, :ICU, :DED), :_INC)
    c2 = Symbol.((:PRO, :LIK, :HES, :ANT, :UNDEFV), :_INC)
    c3 = Symbol.((:PRO, :LIK, :HES, :ANT, :UNDEFV), :_PREV)
   for (k, df) in mydfs

        if k != "allv"
            println("saving dataframe sim level: $k")
            # simulation level, save file per health status, per age group
            #for c in vcat(c1..., c2...)
            for c in vcat(c1...)
                udf = unstack(df, :time, :sim, c) 
                fn = string("$(folderprefix)/simlevel_", lowercase(string(c)), "_", k, ".dat")
                CSV.write(fn, udf)
            end
            println("saving dataframe time level: $k")
        else
            println("saving dataframe sim level: $k")
            # simulation level, save file per health status, per age group
            #for c in vcat(c1..., c2...)
            for c in vcat(c2...)
                udf = unstack(df, :time, :sim, c) 
                fn = string("$(folderprefix)/simlevel_", lowercase(string(c)), "_", k, ".dat")
                CSV.write(fn, udf)
            end
            for c in vcat(c3...)
                udf = unstack(df, :time, :sim, c) 
                fn = string("$(folderprefix)/simlevel_", lowercase(string(c)), "_", k, ".dat")
                CSV.write(fn, udf)
            end
            println("saving dataframe time level: $k")
        end

    end


    R01 = [cdr[i].R0 for i=1:nsims]
    vac_number = [cdr[i].vac_number for i=1:nsims]

    writedlm(string(folderprefix,"/R01.dat"),R01)
    writedlm(string(folderprefix,"/vac_number.dat"),vac_number)

    return mydfs
end



function create_folder_new(ip::cv.ModelParameters, change_l::Bool)
    #strategy = ip.apply_vac_com == true ? "S1" : "S2"
    if !change_l
        RF = string("/data/thomas/homophily/results_", ip.file_index, "_", ip.β, "_", ip.h, "_", ip.κ, "_", ip.function_homophily) ## 
        #RF = string("./outputs/results_", ip.file_index, "_", ip.β, "_", ip.h, "_", ip.b, "_", ip.b_value) ## 
            #RF = string("./outputs/results_", ip.file_index, "_", ip.β, "_", ip.h, "_", ip.b, "_", ip.b_value,"_", ip.probrec)
    else
        RF = string("/data/thomas/homophily/results_", ip.file_index, "_", ip.β, "_", ip.h, "_", ip.κ, "_", ip.function_homophily, "_", ip.L) ## 
    end
    #RF = string("./outputs/results_", ip.β, "_", ip.h, "_", ip.b, "_", ip.b_value) ## 
    if !Base.Filesystem.isdir(RF)
        Base.Filesystem.mkpath(RF)
    end
    return RF
end

function run_param(fidx::Int64, beta::Float64, hh::Float64, kapa::Float64 = 2.0, vac_rate::Int64 = 0, group_init::Union{Nothing, Int8} = nothing, nsims::Int64 = 1000, multib::Vector{Float64} = ones(Float64, 4), xi::Float64 = 1.0, law::Symbol = :power, ll::Vector{Float64} = [0.7, 2.0, 0.2], change_l::Bool = false, im_p::Int8 = Int8(0))
     GC.gc()
    @everywhere ip = cv.ModelParameters(β=$beta,h = $(hh), κ = $kapa, file_index = $fidx,
    vaccination_rate = $vac_rate, ξ = $xi,
    groupinitial = $group_init, mult_b = $multib, vac_immunity_period = $im_p, function_homophily=$(QuoteNode(law)), L = $ll[1], alphan = $ll[2], k = $ll[3])
    folder = create_folder_new(ip, change_l)

    run(ip,nsims,folder)
    
end