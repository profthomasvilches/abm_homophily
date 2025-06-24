using Distributed

using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using ClusterManagers
using Dates
using DelimitedFiles

ENV["JULIA_WORKER_TIMEOUT"] = "180"

addprocs(ClusterManagers.SlurmManager(250), N=8, topology=:master_worker, exeflags="--project=Project.toml"; W="300")
#addprocs(10)
@everywhere using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@everywhere include("abm_behavior.jl")
@everywhere const cv=abmbehavior


function run(myp::cv.ModelParameters, nsims=1000, folderprefix="./")
    println("starting $nsims simulations...\nsave folder set to $(folderprefix)")
    dump(myp)
    
    cdr = pmap(1:nsims) do x                 
            cv.runsim(x, myp)
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
    c1 = Symbol.((:LAT, :ASYMP, :INF, :DED), :_INC)
    #c2 = Symbol.((:LAT, :ASYMP, :INF, :IISO, :HOS, :ICU, :DED), :_PREV)
    #c1 = Symbol.((:LAT, :HOS, :ICU, :DED), :_INC)
    c2 = Symbol.((:PRO, :LIK, :HES, :ANT, :UNDEFV), :_INC)
    c3 = Symbol.((:PRO, :LIK, :HES, :ANT, :UNDEFV), :_PREV)
    R0::Float64 = 0.0

    if ip.calibrating 
        c = :LAT_INC
        udf = unstack(mydfs["all"], :time, :sim, c)
        R0 = mean(sum(Matrix(udf[2:end, 2:end]), dims = 1))
        println("###########################")
        println("Calibrating: the R0 is $R0")
        println("###########################")
        return 0
    end 

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


    vac_number = [cdr[i].vac_number for i=1:nsims]
   
    writedlm(string(folderprefix,"/vac_number.dat"),vac_number)

    return mydfs
end



function create_folder_new(ip::cv.ModelParameters)
    #strategy = ip.apply_vac_com == true ? "S1" : "S2"
    if isnothing(ip.probrec)
        #RF = string("/data/thomas/homophily/results_", ip.file_index, "_", ip.β, "_", ip.h, "_", ip.b, "_", ip.b_value) ## 
        RF = string("./outputs/results_", ip.file_index, "_", ip.β, "_", ip.h, "_", ip.b, "_", ip.b_value) ## 
    else
        #RF = string("/data/thomas/homophily/results_", ip.file_index, "_", ip.β, "_", ip.h, "_", ip.b, "_", ip.b_value,"_", ip.probrec)
        RF = string("./outputs/results_", ip.file_index, "_", ip.β, "_", ip.h, "_", ip.b, "_", ip.b_value,"_", ip.probrec)
    end
    #RF = string("./outputs/results_", ip.β, "_", ip.h, "_", ip.b, "_", ip.b_value) ## 
    if !Base.Filesystem.isdir(RF)
        Base.Filesystem.mkpath(RF)
    end
    return RF
end

function run_param(fidx::Int64, beta::Float64, hh::Float64, bb::Float64 = 0.5, scen::Symbol = :prob, vac_rate::Int64 = 0, ve::Float64 = 0.0, ves::Float64 = 0.0, prob_rec::Union{Nothing, Float64} = nothing, group_init::Union{Nothing, Int8} = nothing, pops::Int64 = 10000, nsims::Int64 = 1000, calib::Bool = false, mt::Int64 = 200)
    
    @everywhere ip = cv.ModelParameters(β=$beta,h = $(hh), b = $bb, b_value=$(QuoteNode(scen)), file_index = $fidx,
    vaccination_rate = $vac_rate, popsize = $pops, vaccine_eff = $ve, vaccine_eff_symp = $ves, probrec = $prob_rec,
    groupinitial = $group_init, calibrating = $calib, modeltime = $mt)
    folder = create_folder_new(ip)

    run(ip,nsims,folder)
    
end