using Distributed
using Base.Filesystem
using DataFrames
using CSV
using Query
using Statistics
using UnicodePlots
using ClusterManagers
using Dates
using DelimitedFiles

## load the packages by covid19abm

#using covid19abm

#addprocs(2, exeflags="--project=.")


#@everywhere using covid19abm

addprocs(SlurmManager(250), N=9, topology=:master_worker, exeflags="--project=.")
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
    DelimitedFiles.writedlm("$(folderprefix)/infectors.dat", [cdr[i].infectors for i = 1:nsims])    

    ## write contact numbers
    #writedlm("$(folderprefix)/ctnumbers.dat", [cdr[i].ct_numbers for i = 1:nsims])    
    ## stack the sims together
    allag = vcat([cdr[i].a  for i = 1:nsims]...)
    allv = vcat([cdr[i].allv  for i = 1:nsims]...)
    ag1 = vcat([cdr[i].g1 for i = 1:nsims]...)
    ag2 = vcat([cdr[i].g2 for i = 1:nsims]...)
    ag3 = vcat([cdr[i].g3 for i = 1:nsims]...)
    ag4 = vcat([cdr[i].g4 for i = 1:nsims]...)
    ag5 = vcat([cdr[i].g5 for i = 1:nsims]...)
    ag6 = vcat([cdr[i].g5 for i = 1:nsims]...)
    ag7 = vcat([cdr[i].g5 for i = 1:nsims]...)
    mydfs = Dict("all" => allag, "allv" => allv, "ag1" => ag1, "ag2" => ag2, "ag3" => ag3, "ag4" => ag4, "ag5" => ag5, "ag6" => ag6, "ag7" => ag7)
    #mydfs = Dict("all" => allag)
    
    ## save at the simulation and time level
    ## to ignore for now: miso, iiso, mild 
    c1 = Symbol.((:LAT, :ASYMP, :INF, :DED), :_INC)
    #c2 = Symbol.((:LAT, :ASYMP, :INF, :IISO, :HOS, :ICU, :DED), :_PREV)
    #c1 = Symbol.((:LAT, :HOS, :ICU, :DED), :_INC)
    c2 = Symbol.((:PRO, :LIK, :HES, :ANT, :UNDEF), :_INC)
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
            println("saving dataframe time level: $k")
        end
        
        



    end

    return mydfs
end



function savestr(p::cv.ModelParameters, custominsert="/", customstart="")
    datestr = (Dates.format(Dates.now(), dateformat"mmdd_HHMM"))
    ## setup folder name based on model parameters
    taustr = replace(string(p.τmild), "." => "")
    fstr = replace(string(p.fmild), "." => "")
    rstr = replace(string(p.β), "." => "")
    prov = replace(string(p.prov), "." => "")
    eldr = replace(string(p.eldq), "." => "")
    eldqag = replace(string(p.eldqag), "." => "")     
    fpreiso = replace(string(p.fpreiso), "." => "")
    tpreiso = replace(string(p.tpreiso), "." => "")
    fsev = replace(string(p.fsevere), "." => "")    
    frelasymp = replace(string(p.frelasymp), "." => "")
    strat = replace(string(p.ctstrat), "." => "")
    pct = replace(string(p.fctcapture), "." => "")
    cct = replace(string(p.fcontactst), "." => "")
    idt = replace(string(p.cidtime), "." => "") 
    tback = replace(string(p.cdaysback), "." => "")     
    fldrname = "/data/covid19abm/simresults/$(custominsert)/$(customstart)_$(prov)_strat$(strat)_pct$(pct)_cct$(cct)_idt$(idt)_tback$(tback)_fsev$(fsev)_tau$(taustr)_fmild$(fstr)_q$(eldr)_qag$(eldqag)_relasymp$(frelasymp)_tpreiso$(tpreiso)_preiso$(fpreiso)/"
    mkpath(fldrname)
end

function create_folder_new(ip::cv.ModelParameters)
    #strategy = ip.apply_vac_com == true ? "S1" : "S2"
    RF = string("/data/thomas/homophily/") ## 
    if !Base.Filesystem.isdir(RF)
        Base.Filesystem.mkpath(RF)
    end
    return RF
end

function run_param(fidx::Int16, beta::Float64, hh::Float64, bb::Float64 = 0.5; scen::Symbol = :prob, )
    
    @everywhere ip = cv.ModelParameters(β=$beta,h = $(hh), b = $bb, b_value = $scen, file_index = $fidx)
    folder = create_folder_new(ip)

    println("$h_i $isos")
    run(ip,nsims,folder)
    
end