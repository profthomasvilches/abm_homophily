module abmbehavior

using Base
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
@enum HEALTH SUS LAT PRE ASYMP INF REC DED UNDEF # 
@enum VACS PRO LIK HES ANT UNDEFV


Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUS
    health_status::HEALTH = SUS
    vac_behavior::VACS = UNDEFV
    swap::HEALTH = UNDEF
    swap_status::HEALTH = UNDEF
    sickfrom::HEALTH = UNDEF
    sickby::Int64 = -1
    nextday_meetcnt::UInt8 = 0 ## how many contacts for a single day
    age::Int16   = 0    # in years. don't really need this but left it incase needed later
    ag::Int16   = 0
    tis::Int16   = 0   # time in state 
    exp::Int16   = 0   # max statetime
    dur::NTuple{4, Int8} = (0, 0, 0, 0)   # Order: (latents, asymps, pres, infs) TURN TO NAMED TUPS LATER
    doi::Int16   = 999   # day of infection.
    iso::Bool = false  ## isolated (limited contacts)
    isovia::Symbol = :null ## isolated via quarantine (:qu), preiso (:pi), intervention measure (:im), or contact tracing (:ct)    
    contacts_vac::Vector{UInt8} = UInt8.([0;0;0;0;0;0;0;0]) # PRO NV; LIK; HES; ANTI; REC vac; DED vac; PRO VAC; VAC  
    connections::Vector{UInt32} = UInt32.([0])
    #comorbidity::Int8 = 0 ##does the individual has any comorbidity?
    
    vac_status::Int8 = 0 ##
    wentto::Int8 = 0
    incubationp::Int16 = 0

    ag_new::Int16 = -1
    
   
    recovered::Bool = false 
    daysisolation::Int64 = 999 
    daysinf::Int64 = 999 
    isofalse::Bool = false
    
    betas_vac::Vector{Float32} = Float32.([0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]) # bs, bh, bl, ba, be
   
    totaldaysiso::Int16 = 0  
end

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0345  
    h::Float64 = 0.5     
    seasonal::Bool = false ## seasonal betas or not
    prov::Symbol = :usa
    calibration::Bool = false
    start_several_inf::Bool = true
    modeltime::Int64 = 200
    initialinf::Int64 = 1
   
    fsevere::Float64 = 1.0 #
    frelasymp::Float64 = 0.26 ## relative transmission of asymptomatic
    vaccine_eff::Float16 = 0.0   ## change this to Float32 typemax(Float32) typemax(Float64)
    vaccine_eff_symp::Float16 = 0.0   ## change this to Float32 typemax(Float32) typemax(Float64)
    vaccination_rate::Int64 = 10
    
    file_index::Int16 = 0
    α::Float64 = 1.0
    b_value::Symbol = :prob
    b = 0.5
    # hosp_red::Float64 = 3.1 # Taiye: We can add this if we decide to include hospitalizations.
    isolation_days::Int64 = 5
    probrec::Union{Nothing, Float64} = nothing
    groupinitial::Union{Nothing, Int8} = nothing
    κ::Float64 = 2.0
end


Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

include("matrices.jl")
## constants 
const popsize = 50000
const humans = Array{Human}(undef, popsize) 
const p = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99]
const v_basis = @SVector [PRO, LIK, HES, ANT]
const pos_s = collect(1:popsize)
const behaviors = falses(popsize)

export popsize, ModelParameters, HEALTH, VACS, Human, humans, BETAS, v_basis, pos_s

function runsim(simnum, ip::ModelParameters)
    GC.gc()
    # function runs the `main` function, and collects the data as dataframes. 
    hmatrix, vmatrix, hh1,vac_number, initial_state = main(ip,simnum)            

    #Get the R0
    
    R01 = length(findall(k -> k.sickby in hh1, humans))/length(hh1)
    
    ###use here to create the vector of comorbidity
    # get simulation age groups
    #ags = [x.ag for x in humans] # store a vector of the age group distribution 
    #ags = [x.ag_new for x in humans] # store a vector of the age group distribution 
    
    all1 = _collectdf(hmatrix)
    
    allv = _collectdfv(vmatrix)


    age_groups = [0:14, 15:24, 25:34, 35:44, 45:54, 55:64, 65:999]
    ags = map(x->findfirst(y-> x.age in y, age_groups),humans) # store a vector of the age group distribution 
    #spl = _splitstate(hmatrix, ags)
    # ag1 = _collectdf(spl[1])
    # ag2 = _collectdf(spl[2])
    # ag3 = _collectdf(spl[3])
    # ag4 = _collectdf(spl[4])
    # ag5 = _collectdf(spl[5])
    # ag6 = _collectdf(spl[6])
    # ag7 = _collectdf(spl[7])
    insertcols!(all1, 1, :sim => simnum);insertcols!(allv, 1, :sim => simnum); # insertcols!(ag1, 1, :sim => simnum); insertcols!(ag2, 1, :sim => simnum); 
    #insertcols!(ag3, 1, :sim => simnum); insertcols!(ag4, 1, :sim => simnum); insertcols!(ag5, 1, :sim => simnum);
    #insertcols!(ag6, 1, :sim => simnum); insertcols!(ag7, 1, :sim => simnum);
    

    pos = findall(y-> y in (11,22,33),hmatrix[:,end])

    vector_ded::Vector{Int64} = zeros(Int64,100)

    for i = pos
        x = humans[i]
        vector_ded[(x.age+1)] += 1
    end

    return (a=all1, allv = allv, #g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5,g6=ag6,g7=ag7,
    vector_dead=vector_ded, R0 = R01, vac_number = vac_number)
end
export runsim

function main(ip::ModelParameters,sim::Int64)
    Random.seed!(sim*726)
    ## datacollection            
    # matrix to collect model state for every time step

    # reset the parameters for the simulation scenario
    reset_params(ip)  #logic: outside "ip" parameters are copied to internal "p" which is a global const and available everywhere. 

    popsize == 0 && error("no population size given")
    
    hmatrix = zeros(Int16, popsize, p.modeltime)
    vmatrix = zeros(Int16, popsize, p.modeltime)
    initialize() # initialize population
    
    #h_init::Int64 = 0
    # insert initial infected agents into the model
    # and setup the right swap function. 

    #create herd immunity
    # herd_immu_dist_4(sim,1) # Taiye: We are not considering herd immunity.

    # split population in agegroups 
    grps = get_ag_dist()
    
    
    Bj = map(x-> get_age_behav(x), grps)
    Bj = map(x-> [Bj[j][x] for j in 1:length(Bj)], 1:length(v_basis))

    vac_number = zeros(Int64, p.modeltime)
    #insert one infected in the latent status in age group 4
    insert_infected(LAT, p.initialinf, 4)

    # h_init1 = findall(x->x.health_status  in (LAT,MILD,INF,PRE,ASYMP), humans) # Taiye: MILD is unnecessary.
    h_init1 = findall(x->x.health_status  in (LAT,INF,PRE,ASYMP), humans)
    
    ## save the preisolation isolation parameters
    #we need the workplaces to get the next days counts
    for x in humans
        get_nextday_counts(x)
        cnts = Int(x.nextday_meetcnt)
        cnts == 0 && continue # skip person if no contacts
        #general population contact
        gpw = Int.(round.(cm[x.ag].*cnts))
         
        x.connections = return_contacts(x, gpw, Vector(Bj[Int(x.vac_behavior)+1]))
        #got back to a more random network
        # for j in eachindex(x.connections)
            
        #     if rand() < p.h
        #         x.connections[j] = rand(grps[humans[x.connections[j]].ag])
        #     end
        # end
    end
    
    initial_state = [Int(humans[i].vac_behavior) for i in eachindex(humans)]

    # start the time loop
    for st = 1:p.modeltime
        
        
        vaccination(sim)
        vac_number[st] = sum([x.vac_status for x in humans])
        #     for x in humans
        # #        if x.iso && !(x.health_status in (HOS,ICU,DED)) # Taiye: Depends on whether we are considering HOS, ICU and DED.
        #         if x.iso && !(x.health_status in (DED)) #&& !(x.health_status in (HOS,ICU,DED))
        #             x.totaldaysiso += 1
        #         end
        #     end
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        _get_model_state2(st, vmatrix) ## this datacollection needs to be at the start of the for loop
        dyntrans(st, grps,sim)
        sw = time_update(vac_number[st]) ###update the system
        
        # end of day
    end
    

    
    
    return hmatrix, vmatrix, h_init1, vac_number, initial_state## return the model state as well as the age groups. 
end
export main

function vaccination(sim::Int64)
    if p.vaccination_rate == 0
        return 0
    end
    rng = MersenneTwister(1234*sim)
    # pos = findall(x -> x.vac_status == 0 && x.vac_behavior ∈ (UNDEFV, PRO), humans)
    pos = Iterators.take(shuffle(rng, findall(x -> x, behaviors)), p.vaccination_rate)

    for i in pos
        humans[i].vac_status = 1
        behaviors[i] = false
    end
    #remaining_doses = max(0, p.vaccination_rate+remaining_doses-length(pos))
    
    return 0#remaining_doses
end

function reset_params(ip::ModelParameters)
    # the p is a global const
    # the ip is an incoming different instance of parameters 
    # copy the values from ip to p. 
    for x in propertynames(p)
        setfield!(p, x, getfield(ip, x))
    end

    # reset the contact tracing data collection structure
    #= for x in propertynames(ct_data)
        setfield!(ct_data, x, 0)
    end

    =#    # resize and update the BETAS constant array
    #init_betas()

    # resize the human array to change population size
    #resize!(humans, popsize)
    
end
export reset_params, reset_params_default


## Data Collection/ Model State functions
function _get_model_state(st, hmatrix)
    # collects the model state (i.e. agent status at time st)
    for i=1:length(humans)
        hmatrix[i, st] = Int(humans[i].health)
    end    
end
export _get_model_state

## Data Collection/ Model State functions
function _get_model_state2(st, vmatrix)
    # collects the model state (i.e. agent status at time st)
    for i=1:length(humans)
        vmatrix[i, st] = Int(humans[i].vac_behavior)
    end    
end
export _get_model_state2

function _collectdf(hmatrix)
    ## takes the output of the humans x time matrix and processes it into a dataframe
    #_names_inci = Symbol.(["lat_inc", "mild_inc", "miso_inc", "inf_inc", "iiso_inc", "hos_inc", "icu_inc", "rec_inc", "ded_inc"])    
    #_names_prev = Symbol.(["sus", "lat", "mild", "miso", "inf", "iiso", "hos", "icu", "rec", "ded"])
    mdf_inc, mdf_prev = _get_incidence_and_prev(hmatrix)
    mdf = hcat(mdf_inc, mdf_prev)    
    _names_inc = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_INC"))
    _names_prev = Symbol.(string.((Symbol.(instances(HEALTH)[1:end - 1])), "_PREV"))
    _names = vcat(_names_inc..., _names_prev...)
    datf = DataFrame(mdf, _names)
    insertcols!(datf, 1, :time => 1:p.modeltime) ## add a time column to the resulting dataframe
    return datf
end

function _collectdfv(vmatrix)
    ## takes the output of the humans x time matrix and processes it into a dataframe
    #_names_inci = Symbol.(["lat_inc", "mild_inc", "miso_inc", "inf_inc", "iiso_inc", "hos_inc", "icu_inc", "rec_inc", "ded_inc"])    
    #_names_prev = Symbol.(["sus", "lat", "mild", "miso", "inf", "iiso", "hos", "icu", "rec", "ded"])
    mdf_inc, mdf_prev = _get_incidence_and_prev_v(vmatrix)
    mdf = hcat(mdf_inc, mdf_prev)    
    _names_inc = Symbol.(string.((Symbol.(instances(VACS)[1:end])), "_INC"))
    _names_prev = Symbol.(string.((Symbol.(instances(VACS)[1:end])), "_PREV"))
    _names = vcat(_names_inc..., _names_prev...)
    datf = DataFrame(mdf, _names)
    insertcols!(datf, 1, :time => 1:p.modeltime) ## add a time column to the resulting dataframe
    return datf
end


function _splitstate(hmatrix, ags)
    #split the full hmatrix into 4 age groups based on ags (the array of age group of each agent)
    #sizes = [length(findall(x -> x == i, ags)) for i = 1:4]
    matx = []#Array{Array{Int64, 2}, 1}(undef, 4)
    for i = 1:maximum(ags)#length(agebraks)
        idx = findall(x -> x == i, ags)
        push!(matx, view(hmatrix, idx, :))
    end
    return matx
end
export _splitstate

function _get_incidence_and_prev(hmatrix)
    cols = instances(HEALTH)[1:end - 1] ## don't care about the UNDEF health status
    inc = zeros(Int64, p.modeltime, length(cols))
    pre = zeros(Int64, p.modeltime, length(cols))
    for i = 1:length(cols)
        inc[:, i] = _get_column_incidence(hmatrix, cols[i])
        pre[:, i] = _get_column_prevalence(hmatrix, cols[i])
    end
    return inc, pre
end


function _get_incidence_and_prev_v(vmatrix)
    cols = instances(VACS)[1:end] ## don't care about the UNDEF health status
    inc = zeros(Int64, p.modeltime, length(cols))
    pre = zeros(Int64, p.modeltime, length(cols))
    for i = 1:length(cols)
        inc[:, i] = _get_column_incidence(vmatrix, cols[i])
        pre[:, i] = _get_column_prevalence(vmatrix, cols[i])
    end
    return inc, pre
end
# Checkpoint

function _get_column_incidence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for r in eachrow(hmatrix)
        idx = findall(x-> r[x] == inth && r[x] != r[x-1],2:length(r))
        idx = idx .+ 1
        #idx = findfirst(x -> x == inth, r)
        if idx !== nothing
            for i in idx 
                timevec[i] += 1
            end
        end
    end
    return timevec
end

function _get_column_prevalence(hmatrix, hcol)
    inth = Int(hcol)
    timevec = zeros(Int64, p.modeltime)
    for (i, c) in enumerate(eachcol(hmatrix))
        idx = findall(x -> x == inth, c)
        if idx !== nothing
            ps = length(c[idx])    
            timevec[i] = ps    
        end
    end
    return timevec
end

export _collectdf, _get_incidence_and_prev, _get_column_incidence, _get_column_prevalence

## initialization functions 
function get_province_ag(prov) 
    ret = @match prov begin
        :usa => Distributions.Categorical(@SVector [0.059444636404977,0.188450296592341,0.396101793107413,0.189694011721906,0.166309262173363])
        # :ontario => Distributions.Categorical(@SVector [0.04807822, 0.10498712, 0.12470340, 0.14498051, 0.13137129, 0.12679091, 0.13804896, 0.10292032, 0.05484776, 0.02327152])
        # :canada => Distributions.Categorical(@SVector [0.04922255,0.10812899,0.11792442,0.13956709,0.13534216,0.12589012,0.13876094,0.10687438,0.05550450,0.02278485])
        _ => error("shame for not knowing your canadian provinces and territories")
    end       
    return ret  
end
export get_province_ag


function initialize() 
    agedist = get_province_ag(p.prov)
    #agebraksnew = [0:4,5:14,15:24,25:34,35:44,45:54,55:64,65:74,75:84,85:99]
    
    age_break_behav = [0:17, 18:24,25:49,50:64,65:99]
    prob_behav = [[0.271;0.241;0.303;0.084],[0.38;0.214;0.205;0.088],[0.271;0.241;0.303;0.084],[0.322;0.302;0.223;0.053],[0.471;0.288;0.131;0.04]]
    a = [4;19;49;64;79;999]
     
    for i = 1:popsize 
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        agn = rand(agedist)
        x.age = rand(agebraks[agn]) 
        x.ag = findfirst(y-> x.age in y, agebraks)
        
        g = findfirst(y->y>=x.age,a)
        
        
        x.ag_new = g
        
        if x.age >= 0
            g = findfirst(y->x.age ∈ y, age_break_behav)
            x.vac_behavior = sample(v_basis, Weights(prob_behav[g]))
            behaviors[x.idx] = x.vac_behavior == PRO ? true : false
        end
        x.exp = 999  ## susceptible people don't expire.
        x.betas_vac = getting_b_values(p)
        #x.dur = sample_epi_durations() # sample epi periods   

    end
end
export initialize

function getting_b_values(p::ModelParameters)

    if p.b_value == :fixed
        return [1-p.b; p.b; 1-p.b; p.b; p.b]# bs, bh, bl, ba, be
    elseif p.b_value == :prob

        b1 = round(rand(Distributions.Uniform(0.08, 0.12)), digits = 5)#p.b
        b2 = round(rand(Distributions.Uniform(0.05, 0.08)), digits = 5)#round(rand(Distributions.Beta(2, 25)), digits = 5) #bh
        b3 = round(rand(Distributions.Uniform(0.02, 0.05)), digits = 5)#round(rand(Distributions.Beta(1.5, 20)), digits = 5) #ba
        b4 = round(rand(Distributions.Uniform(0.03, 0.06)), digits = 5)#round(rand(Distributions.Beta(1, 20)), digits = 5) #bl
        b5 = isnothing(p.probrec) ? round(rand(Distributions.Uniform(0.01, 0.03)), digits = 5) : Float32(p.probrec)  #round(rand(Distributions.Beta(2, 20)), digits = 5) : Float32(p.probrec) #be
        b6 = round(rand(Distributions.Uniform(0.01, 0.03)), digits = 5) #bpl
        b7 = round(rand(Distributions.Uniform(0.06, 0.09)), digits = 5) #bpl
        return Float32.([b1; b2; b3; b4;  b5; b6; b7])
    else
        error("no strategy set for b values")
    end

end

function get_ag_dist() 
    # splits the initialized human pop into its age groups
    grps =  map(x -> findall(y -> y.ag == x, humans), 1:length(agebraks)) 
    return grps
end

function insert_infected(health, num, ag) 
    ## inserts a number of infected people in the population randomly
    ## this function should resemble move_to_inf()
    if isnothing(p.groupinitial)
        l = findall(x -> x.health == SUS && x.ag == ag, humans)
    else
        l = findall(x -> x.health == SUS && x.vac_behavior == v_basis[p.groupinitial] && x.ag == ag, humans)
    end
    if length(l) > 0 && num < length(l)
        h = sample(l, num; replace = false)
        @inbounds for i in h 
            x = humans[i]
            x.dur = sample_epi_durations(x)
            
            if health == PRE
                x.swap = health
                x.swap_status = PRE
                x.daysinf = x.dur[1]+1
                x.wentto = 1
                move_to_pre(x) ## the swap may be asymp, mild, or severe, but we can force severe in the time_update function
            elseif health == LAT
                x.swap = health
                x.swap_status = LAT
                x.daysinf = rand(1:x.dur[1])
                move_to_latent(x)
            
            # Taiye: We are not currently considering MILD**.
            #elseif health == MILD
                #   x.swap = health
                #  x.swap_status = MILD
                # x.wentto = 1
                #x.daysinf = x.dur[2]+1
                #move_to_mild(x)
            elseif health == INF
                x.swap = health
                x.swap_status = INF
                x.wentto = 1
                move_to_inf(x)
            elseif health == ASYMP
                x.swap = health
                x.swap_status = ASYMP
                x.wentto = 2
                move_to_asymp(x)
            
            
            elseif health == REC 
                x.swap = health
                x.swap_status = REC
                x.wentto = rand(1:2)
                move_to_recovered(x)

            else 
                error("can not insert human of health $(health)")
            end
            
            
            x.sickfrom = INF # this will add +1 to the INF count in _count_infectors()... keeps the logic simple in that function.    
            
        end
    end    
    return h
end
export insert_infected

function get_alphas()
    αv = p.κ*rand(Distributions.Beta(2,3))
    αsv = p.κ*rand(Distributions.Beta(3,4))
    αsl = p.κ*rand(Distributions.Beta(3,3))
    αsh = p.κ*rand(Distributions.Beta(4,3))
    αsa = p.κ*rand(Distributions.Beta(3,2))
    αrv = p.κ*rand(Distributions.Beta(3,2))

    return αv, αsv, αsl, αsh, αsa, αrv
end

function time_update(nvac::Int64)

    # counters to calculate incidence
    αv, αsv, αsl, αsh, αsa, αrv = get_alphas()
    contagens = countmap([x.vac_behavior for x in humans])
    contagem_vetor = [get(contagens, i, 0) for i in v_basis]
    
    ninfvac = length(findall(x -> x.vac_status == 1 && x.wentto == 1, humans))
    
    contagem_vetor = [nvac; contagem_vetor; ninfvac]
    contagem_vetor[2] = contagem_vetor[2]-nvac

    for x in humans 
        x.tis += 1 
        x.doi += 1 # increase day of infection. variable is garbage until person is latent
         
        x.daysinf += 1
        
        if x.tis >= x.exp             
            if x.swap_status == LAT
                move_to_latent(x) 
            elseif x.swap_status == PRE
                move_to_pre(x)
            elseif x.swap_status == ASYMP
                move_to_asymp(x)
            elseif x.swap_status == INF
                move_to_inf(x)
            elseif x.swap_status == REC
                move_to_recovered(x)
            elseif x.swap_status == DED
                move_to_dead(x)
            else
                dump(x)
                error("swap expired, but no swap set.")
            end
        end

        # behavioral change
        update_behavior(x, [αv, αsv, αsl, αsh, αsa, αrv], contagem_vetor)
        
        x.contacts_vac = UInt8.([0;0;0;0;0;0;0;0;0])
    
        #get_nextday_counts(x)
        
    end

end
export time_update


function update_behavior(x::Human, alphas::Vector{Float64}, ninds::Vector{Int64})

    #? create a vector of Pv, Pa...
    #? 

    αv, αsv, αsl, αsh, αsa, αrv = alphas
    V, Sv, L, H, A, Rv = ninds


    if x.vac_behavior == LIK
        # probability of getting PRO
        n1 = x.contacts_vac[7]+p.α*x.contacts_vac[1]
        prob1 = (1-x.betas_vac[1])^n1 # bs, bh, bl, ba, be

        # probability of getting hesitant
        n2 = x.contacts_vac[3]+x.contacts_vac[4]
        prob2 = (1-x.betas_vac[2])^n2

        # probability related to vaccinated and recovered
        n3 = x.contacts_vac[5]
        prob3 = (1-x.betas_vac[5])^n3
        #? go to pro or Hes
        #? 1-prob1*prob2*prob3

        #! global influence
         # probability of getting PRO
        n1 = (αv*V+αsv*Sv)/popsize
        prob11 = (1-x.betas_vac[1])^n1 # bs, bh, bl, ba, be

        # probability of getting hesitant
        n2 = (αsh*H+αsa*A)/popsize
        prob21 = (1-x.betas_vac[2])^n2

        # probability related to vaccinated and recovered
        n3 = αrv*Rv/popsize
        prob31 = (1-x.betas_vac[5])^n3


        if rand() < (1-prob1*prob11)/((1-prob1*prob11)+(1-prob2*prob3*prob21*prob31))
            if rand() < 1-prob1*prob11 #? Is the individual changing?
                x.vac_behavior = PRO
                behaviors[x.idx] = true
            end
        else
            if rand() < 1-prob2*prob3*prob21*prob31
                
                x.vac_behavior = HES
            end
        end
    elseif x.vac_behavior == HES
          # probability of getting LIK
          n1 = x.contacts_vac[7]+p.α*x.contacts_vac[1] 
          prob1 = (1-x.betas_vac[7])^n1 # bs, bh, bl, ba, be
  
          # probability of getting anti
          n2 = x.contacts_vac[4]
          prob2 = (1-x.betas_vac[4])^n2


          n1 = (αsl*L+αsv*Sv+αv*V)/popsize 
          prob11 = (1-x.betas_vac[7])^n1 # bs, bh, bl, ba, be
  
          # probability of getting anti
          n2 = αsa*A/popsize
          prob21 = (1-x.betas_vac[4])^n2
  
        if rand() < (1-prob1*prob11)/((1-prob1*prob11)+(1-prob2*prob21))
            if rand() < 1-prob1*prob11
                x.vac_behavior = LIK
            end
        else
            if rand() < 1-prob2*prob21
                x.vac_behavior = ANT
            end
        end 
    elseif x.vac_behavior == ANT
        # probability of getting hes
        n1 = x.contacts_vac[8] # number of vaccinated in contact
        prob1 = (1-x.betas_vac[3])^n1 # bs, bh, bl, ba, be
        n1 = (αsl*L+αsh*H)/popsize # number of vaccinated in contact
        prob11 = (1-x.betas_vac[3])^n1 # bs, bh, bl, ba, be


        if rand() < 1-prob1*prob11 # probability of getting at least one of the changes
            
            x.vac_behavior = HES
            
        end
    elseif x.vac_behavior == PRO && x.vac_status == 0
        # probability of getting likely
        # probability of getting likely
        n1 = x.contacts_vac[Int(LIK)+1]+ x.contacts_vac[Int(HES)+1]+x.contacts_vac[Int(ANT)+1]# number of vaccinated in contact
        prob1 = (1-x.betas_vac[6])^n1 # bs, bh, bl, ba, be

        n1 = (αsl*L+αsh*H)/popsize # number of vaccinated in contact
        prob11 = (1-x.betas_vac[6])^n1 # bs, bh, bl, ba, be


        if rand() < 1-prob1*prob11 # probability of getting at least one of the changes
            
            x.vac_behavior = LIK
            behaviors[x.idx] = false
            
        end

    end

end

function sample_epi_durations(y::Human)
    # when a person is sick, samples the 

    lat_dist = Distributions.truncated(LogNormal(1.434, 0.661), 4, 7) # truncated between 4 and 7
    pre_dist = Distributions.truncated(Gamma(1.058, 5/2.3), 0.8, 3)#truncated between 0.8 and 3
    asy_dist = Gamma(5, 1)
    inf_dist = Gamma((3.2)^2/3.7, 3.7/3.2)

    latents = Int.(round.(rand(lat_dist)))
    y.incubationp = latents
    pres = Int.(round.(rand(pre_dist)))
    latents = latents - pres # ofcourse substract from latents, the presymp periods
    asymps = max(Int.(ceil.(rand(asy_dist))),1)
    infs = max(Int.(ceil.(rand(inf_dist))),1)

    return (latents, asymps, pres, infs)
end

function move_to_latent(x::Human)
    ## transfers human h to the incubation period and samples the duration
    x.health = x.swap
    x.health_status = x.swap_status
    x.doi = 0 ## day of infection is reset when person becomes latent
    x.tis = 0   # reset time in state 
    x.exp = x.dur[1] # get the latent period
   
    #0-18 31 19 - 59 29 60+ 18 going to asymp
    symp_pcts = [0.7, 0.623, 0.672, 0.672, 0.812, 0.812] #[0.3 0.377 0.328 0.328 0.188 0.188]
    age_thres = [4, 19, 49, 64, 79, 999]
    g = findfirst(y-> y >= x.age, age_thres)

 
    if rand() < (symp_pcts[g])*(1-p.vaccine_eff_symp)^x.vac_status
       
        x.swap = PRE
        x.swap_status = PRE
        x.wentto = 1
        
    else
        x.swap = ASYMP
        x.swap_status = ASYMP
        x.wentto = 2
        
    end
    
    ## in calibration mode, latent people never become infectious.
    
end
export move_to_latent

function move_to_asymp(x::Human)
    ## transfers human h to the asymptomatic stage 
    x.health = x.swap  
    x.health_status = x.swap_status
    x.tis = 0 
    x.exp = x.dur[2] # get the presymptomatic period
   
    x.swap = REC
    x.swap_status = REC
    
    # x.iso property remains from either the latent or presymptomatic class
    # if x.iso is true, the asymptomatic individual has limited contacts
end
export move_to_asymp

function move_to_pre(x::Human)
    #θ = (0.95, 0.9, 0.85, 0.6, 0.2)
    
    x.health = x.swap
    x.health_status = x.swap_status
    x.tis = 0   # reset time in state 
    x.exp = x.dur[3] # get the presymptomatic period
    ##########

    x.swap = INF
    x.swap_status = INF
 
    
end
export move_to_pre




function move_to_inf(x::Human)
    ## transfers human h to the severe infection stage for γ days
    ## for swap, check if person will be hospitalized, selfiso, die, or recover
 
    groups = [0:34,35:54,55:69,70:84,85:100]
    gg = findfirst(y-> x.age in y,groups)

    mh = [0.0002; 0.0015; 0.011; 0.0802; 0.381] # death rate for cases.
 
    x.health = x.swap
    x.health_status = x.swap_status
    x.swap = UNDEF
    
    x.tis = 0 
    

       
   # else ## no hospital for this lucky (but severe) individual 
    if rand() < mh[gg]
            x.exp = x.dur[4] 
            x.swap = DED
            x.swap_status = DED
    else 
            x.exp = x.dur[4]  
            x.swap = REC
            x.swap_status = REC
    end
         
   # end
    ## before returning, check if swap is set 
    #x.swap == UNDEF && error("agent I -> ?")
end

function move_to_dead(h::Human)
    # no level of alchemy will bring someone back to life. 
    h.health = h.swap
    h.health_status = h.swap_status
    h.swap = UNDEF
    h.swap_status = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely
end

function move_to_recovered(h::Human)
    h.health = h.swap
    h.health_status = h.swap_status
    
    h.recovered = true

    h.swap = UNDEF
    h.swap_status = UNDEF
    h.tis = 0 
    h.exp = 999 ## stay recovered indefinitely

    #h.iso = false ## a recovered person has ability to meet others
    
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
end


@inline function _get_betavalue(xhealth) 
    #bf = p.β ## baseline PRE
    #length(BETAS) == 0 && return 0
    bf = p.β#BETAS[sys_time]
    # values coming from FRASER Figure 2... relative tranmissibilities of different stages.
    if xhealth == ASYMP
        bf = bf * p.frelasymp #0.11
    elseif xhealth == INF 
        bf = 0.89*bf
    end

    return bf
end
export _get_betavalue

@inline function get_nextday_counts(x::Human)
    # get all people to meet and their daily contacts to recieve
    # we can sample this at the start of the simulation to avoid everyday    
    cnt = 0
    ag = x.ag
    #if person is isolated, they can recieve only 3 maximum contacts

    cnt = min(max(1, rand(negative_binomials(ag))), 254) ##using the contact average for shelter-in
    
    x.nextday_meetcnt = cnt

    

    if x.health_status == DED
        x.nextday_meetcnt = 0
    end
   
    #x.contacts_vac = UInt8.([0;0;0;0;0;0;0;0;0])

    return cnt
end

function calculate_H(gp)
    
    # how many times it appears
    contagens = countmap([humans[i].vac_behavior for i in gp])

    contagem_vetor = [get(contagens, i, 0) for i in v_basis]
    # probabilities based on homophily
    
    Pj = zeros(Float64, length(v_basis), length(v_basis))
    if sum(contagem_vetor) > 0
        for br in v_basis
            H = p.h .* Int.(br .== v_basis)+(1-p.h) .* contagem_vetor./sum(contagem_vetor)
            Pj[Int(br)+1,:] = H ./ sum(H)
        end
    end
    return Pj
end

function get_age_behav(grpi)
    vector = Iterators.take([humans[i].vac_behavior for i in grpi], length(grpi))
    aa = [findall(y -> y == v_basis[i], collect(vector)) for i in eachindex(v_basis)]
    return map(x-> grpi[x], aa)
end

function dyntrans(sys_time, grps,sim)
    #totalmet = 0 # count the total number of contacts (total for day, for all INF contacts)
    #totalinf = 0 # count number of new infected 
    ## find all the people infectious
    #rng = MersenneTwister(246*sys_time*sim)

    shuffle!(pos_s)
    #todo create a vector of vectors: each age group with their respective behaviors
    #? this should improve performance
    
    # go through every infectious person
    for i in pos_s    
            x = humans[i] 
            xhealth = x.health
            perform_contacts(x,grps,xhealth)
            
    end
    #return totalmet, totalinf
end
export dyntrans

function return_contacts(x::Human, gpw::Vector{Int64}, Bji::Vector{Vector{Int64}})

    
    aux = [Int[] for _ in eachindex(Bji)]
    
    for (i, g) in enumerate(gpw)
        if g == 0
            aux[i] = Vector{Int}(undef, 0)
            continue
        end
        if length(Bji[i]) > 0
            aux[i] = sample(Bji[i], g, replace = true)
        end

    end

    

    aux = reduce(vcat, aux)

    return aux

end

function perform_contacts(x::Human,grp_sample::Vector{Vector{Int64}},xhealth::HEALTH)

    for j in x.connections 
        # go through edach person
         
        y = humans[j]
        # getting back to this one
        if rand() < p.h
            y = humans[rand(grp_sample[y.ag])]
        end
    
        if x.vac_status*y.vac_status == 0 && y.health != DED# only gets to it if it is necessary
            
            if y.vac_status == 0 && y.vac_behavior != UNDEFV
                x.contacts_vac[Int(y.vac_behavior)+1] += 1
            elseif y.vac_status == 1
                if y.wentto == 1
                    if  y.health_status == REC
                        x.contacts_vac[5] += 1
                    elseif y.health_status == DED
                        x.contacts_vac[6] += 1
                    end
                else
                    if y.vac_behavior == PRO # redundant
                        x.contacts_vac[7] += 1
                    end
                    x.contacts_vac[8] += 1
                end
               
                
            end

        end
        
        if xhealth == SUS && y.health in (ASYMP, INF) && x.swap == UNDEF
            
            beta = _get_betavalue(y.health)*(1-p.vaccine_eff)^x.vac_status
            
            if rand() < beta
                x.exp = x.tis   ## force the move to latent in the next time step.
                x.sickfrom = y.health ## stores the infector's status to the infectee's sickfrom
                x.sickby = x.sickby < 0 ? y.idx : x.sickby
                x.swap = LAT
                x.swap_status = LAT
                x.daysinf = 0
                x.dur = sample_epi_durations(x)
            end
        end



        
    end  
    
end
function contact_matrix()
    # regular contacts, just with 5 age groups. 
    #  0-4, 5-19, 20-49, 50-64, 65+
    CM = Array{Array{Float64, 1}, 1}(undef, 5)
    CM[1] = [0.25,0.132,0.44,0.144,0.034]
    CM[2] = [0.0264,0.43,0.404,0.108,0.0316]
    CM[3] = [0.03,0.13,0.602,0.179,0.059]
    CM[4] = [0.026,0.086,0.456,0.3,0.132]
    CM[5] = [0.012,0.052,0.303,0.266,0.367]  
   
    return CM
end

# 
# calibrate for 2.7 r0
# 20% selfisolation, tau 1 and 2.

function negative_binomials(ag) 
    ## the means/sd here are calculated using _calc_avgag
    # [0:4, 5:19, 20:49, 50:64, 65:99]
    means = [6.97;9.54;10.96;8.05;4.41]
    sd = [5.22;6.66;8.35;6.86;3.83]
    means = means
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms[ag]
end
#const nbs = negative_binomials()
const cm = contact_matrix()
#export negative_binomials, contact_matrix, nbs, cm

export negative_binomials


function negative_binomials_shelter(ag) 
    ## the means/sd here are calculated using _calc_avgag
    #72% reduction
    means = [1.95; 2.67;3.07; 2.255;1.234]
    sd = [1.461518;1.863781;2.337172;1.920688;1.12]
    means = means
    #sd = sd*mult
    totalbraks = length(means)
    nbinoms = Vector{NegativeBinomial{Float64}}(undef, totalbraks)
    for i = 1:totalbraks
        p = 1 - (sd[i]^2-means[i])/(sd[i]^2)
        r = means[i]^2/(sd[i]^2-means[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end
    return nbinoms[ag]   
end

#const vaccination_days = days_vac_f()
#const vac_rate_1 = vaccination_rate_1()
#const vac_rate_2 = vaccination_rate_2()
## references: 
# critical care capacity in Canada https://www.ncbi.nlm.nih.gov/pubmed/25888116
end