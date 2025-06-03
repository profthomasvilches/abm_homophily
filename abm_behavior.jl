module abmbehavior

using Base
using Parameters, Distributions, StatsBase, StaticArrays, Random, Match, DataFrames
# @enum HEALTH SUS LAT PRE ASYMP MILD MISO INF IISO HOS ICU REC DED UNDEF
@enum HEALTH SUS LAT PRE ASYMP INF REC DED UNDEF # 
@enum VACS PRO LIK HES ANT UNDEFV

#TODO iniciar contatos em zeros 
#!ok
#TODO iniciar status vac_behavior
#!ok
#todo implement b for each individual
#! ok
#todo implement status behavior change
#! ok
#Todo implement contact to recovered (from symptomatic) people that were vaccinated.
#! ok
#todo implement vaccination
#! ok
#todo implement transmission reduction and symptoms based on vaccination
#! ok
#todo implement matrix for behavior

Base.@kwdef mutable struct Human
    idx::Int64 = 0 
    health::HEALTH = SUS
    health_status::HEALTH = SUS
    vac_behavior::VACS = UNDEFV
    swap::HEALTH = UNDEF
    swap_status::HEALTH = UNDEF
    sickfrom::HEALTH = UNDEF
    sickby::Int64 = -1
    nextday_meetcnt::Int16 = 0 ## how many contacts for a single day
    age::Int16   = 0    # in years. don't really need this but left it incase needed later
    ag::Int16   = 0
    tis::Int16   = 0   # time in state 
    exp::Int16   = 0   # max statetime
    dur::NTuple{4, Int8} = (0, 0, 0, 0)   # Order: (latents, asymps, pres, infs) TURN TO NAMED TUPS LATER
    doi::Int16   = 999   # day of infection.
    iso::Bool = false  ## isolated (limited contacts)
    isovia::Symbol = :null ## isolated via quarantine (:qu), preiso (:pi), intervention measure (:im), or contact tracing (:ct)    
    contacts_vac::Vector{Int64} = [0;0;0;0;0;0;0;0] # PRO NV; LIK; HES; ANTI; REC vac; DED vac; PRO VAC; VAC  
    #comorbidity::Int8 = 0 ##does the individual has any comorbidity?
    # Taiye: We are not considering comorbidities at this stage.

    vac_status::Int8 = 0 ##
    wentto::Int8 = 0
    incubationp::Int16 = 0

    got_inf::Bool = false
    
    # herd_im::Bool = false
    # Taiye: We are not considering herd immunity at this stage.

    ag_new::Int16 = -1
    
   
    recovered::Bool = false 
    daysisolation::Int64 = 999 
    daysinf::Int64 = 999 
    isofalse::Bool = false
    
    betas_vac::Vector{Float64} = [0.0; 0.0; 0.0; 0.0; 0.0] # bs, bh, bl, ba, be
   
    totaldaysiso::Int32 = 0  
end

## default system parameters
@with_kw mutable struct ModelParameters @deftype Float64    ## use @with_kw from Parameters
    β = 0.0345  
    h::Float64 = 0.5     
    seasonal::Bool = false ## seasonal betas or not
    popsize::Int64 = 100000
    prov::Symbol = :usa
    calibration::Bool = false
    start_several_inf::Bool = true
    modeltime::Int64 = 200
    initialinf::Int64 = 1
    fmild::Float64 = 0.5  ## percent of people practice self-isolation
    # Taiye: Could be useful later for keeping track of the population in isolation.

    fsevere::Float64 = 1.0 #
    frelasymp::Float64 = 0.26 ## relative transmission of asymptomatic
    vaccine_eff::Float16 = 0.0   ## change this to Float32 typemax(Float32) typemax(Float64)
    vaccine_eff_symp::Float16 = 0.0   ## change this to Float32 typemax(Float32) typemax(Float64)
    vaccination_rate::Int64 = 10
    # herd::Int8 = 0 #typemax(Int32) ~ millions
    # Taiye: We are not considering herd immunity at this stage.

    file_index::Int16 = 0
    α::Float64 = 1.0
    b_value::Symbol = :prob
    b = 0.5
    # hosp_red::Float64 = 3.1 # Taiye: We can add this if we decide to include hospitalizations.
    isolation_days::Int64 = 5
end


Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

include("matrices.jl")
## constants 
const humans = Array{Human}(undef, 0) 
const p = ModelParameters()  ## setup default parameters
const agebraks = @SVector [0:4, 5:19, 20:49, 50:64, 65:99]
#const agebraks_vac = @SVector [0:0,1:4,5:14,15:24,25:44,45:64,65:74,75:100]

export ModelParameters, HEALTH, VACS, Human, humans, BETAS

function runsim(simnum, ip::ModelParameters)
    # function runs the `main` function, and collects the data as dataframes. 
    hmatrix, vmatrix, hh1,vac_number = main(ip,simnum)            

    #Get the R0
    
    R01 = length(findall(k -> k.sickby in hh1,humans))/length(hh1)
    
    ###use here to create the vector of comorbidity
    # get simulation age groups
    #ags = [x.ag for x in humans] # store a vector of the age group distribution 
    #ags = [x.ag_new for x in humans] # store a vector of the age group distribution 
    
    all1 = _collectdf(hmatrix)
    
    allv = _collectdfv(vmatrix)

    age_groups = [0:14, 15:24, 25:34, 35:44, 45:54, 55:64, 65:999]
    ags = map(x->findfirst(y-> x.age in y, age_groups),humans) # store a vector of the age group distribution 
    spl = _splitstate(hmatrix, ags)
    ag1 = _collectdf(spl[1])
    ag2 = _collectdf(spl[2])
    ag3 = _collectdf(spl[3])
    ag4 = _collectdf(spl[4])
    ag5 = _collectdf(spl[5])
    ag6 = _collectdf(spl[6])
    ag7 = _collectdf(spl[7])
    insertcols!(all1, 1, :sim => simnum);insertcols!(allv, 1, :sim => simnum);  insertcols!(ag1, 1, :sim => simnum); insertcols!(ag2, 1, :sim => simnum); 
    insertcols!(ag3, 1, :sim => simnum); insertcols!(ag4, 1, :sim => simnum); insertcols!(ag5, 1, :sim => simnum);
    insertcols!(ag6, 1, :sim => simnum); insertcols!(ag7, 1, :sim => simnum);
    

    pos = findall(y-> y in (11,22,33),hmatrix[:,end])

    vector_ded::Vector{Int64} = zeros(Int64,100)

    for i = pos
        x = humans[i]
        vector_ded[(x.age+1)] += 1
    end

    return (a=all1, g1=ag1, g2=ag2, g3=ag3, g4=ag4, g5=ag5,g6=ag6,g7=ag7,allv = allv,
    vector_dead=vector_ded, R0 = R01, vac_number = vac_number)
end
export runsim

function main(ip::ModelParameters,sim::Int64)
    Random.seed!(sim*726)
    ## datacollection            
    # matrix to collect model state for every time step

    # reset the parameters for the simulation scenario
    reset_params(ip)  #logic: outside "ip" parameters are copied to internal "p" which is a global const and available everywhere. 

    p.popsize == 0 && error("no population size given")
    
    hmatrix = zeros(Int16, p.popsize, p.modeltime)
    vmatrix = zeros(Int16, p.popsize, p.modeltime)
    initialize() # initialize population
    
    #h_init::Int64 = 0
    # insert initial infected agents into the model
    # and setup the right swap function. 

    #create herd immunity
    # herd_immu_dist_4(sim,1) # Taiye: We are not considering herd immunity.

    # split population in agegroups 
    grps = get_ag_dist()
    
    vac_number = zeros(Int64, p.modeltime)
    #insert one infected in the latent status in age group 4
    insert_infected(LAT, p.initialinf, 4)

    # h_init1 = findall(x->x.health_status  in (LAT,MILD,INF,PRE,ASYMP), humans) # Taiye: MILD is unnecessary.
    h_init1 = findall(x->x.health_status  in (LAT,INF,PRE,ASYMP), humans)
    
    ## save the preisolation isolation parameters
    #we need the workplaces to get the next days counts
    for x in humans
        get_nextday_counts(x)
    end
    
    remaining_doses::Int64 = 0
    # start the time loop
    for st = 1:p.modeltime
        
        vac_number[st] = sum([x.vac_status for x in humans])

        remaining_doses = vaccination(remaining_doses)
        for x in humans
    #        if x.iso && !(x.health_status in (HOS,ICU,DED)) # Taiye: Depends on whether we are considering HOS, ICU and DED.
            if x.iso && !(x.health_status in (DED)) #&& !(x.health_status in (HOS,ICU,DED))
                x.totaldaysiso += 1
            end
        end
        _get_model_state(st, hmatrix) ## this datacollection needs to be at the start of the for loop
        _get_model_state2(st, vmatrix) ## this datacollection needs to be at the start of the for loop
        dyntrans(st, grps,sim)
        sw = time_update() ###update the system
        
        # end of day
    end
    

    
    
    return hmatrix, vmatrix, h_init1, vac_number## return the model state as well as the age groups. 
end
export main

function vaccination(remaining_doses::Int64)
    if p.vaccination_rate == 0
        return 0
    end
    pos = findall(x -> x.vac_status == 0 && x.vac_behavior ∈ (UNDEFV, PRO), humans)

    if length(pos) >= p.vaccination_rate+remaining_doses

        pp = sample(pos,  p.vaccination_rate+remaining_doses, replace = false)

        for i in pp
            humans[i].vac_status = 1
        end
        remaining_doses = 0
    else
        
        for i in pos
            humans[i].vac_status = 1
        end

        remaining_doses = p.vaccination_rate+remaining_doses-length(pos)
    end
    return remaining_doses
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
    resize!(humans, p.popsize)
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
       #= :alberta => Distributions.Categorical(@SVector [0.0655, 0.1851, 0.4331, 0.1933, 0.1230])
        :bc => Distributions.Categorical(@SVector [0.0475, 0.1570, 0.3905, 0.2223, 0.1827])
        :manitoba => Distributions.Categorical(@SVector [0.0634, 0.1918, 0.3899, 0.1993, 0.1556])
        :newbruns => Distributions.Categorical(@SVector [0.0460, 0.1563, 0.3565, 0.2421, 0.1991])
        :newfdland => Distributions.Categorical(@SVector [0.0430, 0.1526, 0.3642, 0.2458, 0.1944])
        :nwterrito => Distributions.Categorical(@SVector [0.0747, 0.2026, 0.4511, 0.1946, 0.0770])
        :novasco => Distributions.Categorical(@SVector [0.0455, 0.1549, 0.3601, 0.2405, 0.1990])
        :nunavut => Distributions.Categorical(@SVector [0.1157, 0.2968, 0.4321, 0.1174, 0.0380])
        :pei => Distributions.Categorical(@SVector [0.0490, 0.1702, 0.3540, 0.2329, 0.1939])
        :quebec => Distributions.Categorical(@SVector [0.0545, 0.1615, 0.3782, 0.2227, 0.1831])
        :saskat => Distributions.Categorical(@SVector [0.0666, 0.1914, 0.3871, 0.1997, 0.1552])
        :yukon => Distributions.Categorical(@SVector [0.0597, 0.1694, 0.4179, 0.2343, 0.1187])
        :newyorkcity   => Distributions.Categorical(@SVector [0.064000, 0.163000, 0.448000, 0.181000, 0.144000])=#
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
    
    age_break_behav = [18:24,25:49,50:64,65:99]
    prob_behav = [[0.38;0.214;0.205;0.088],[0.271;0.241;0.303;0.084],[0.322;0.302;0.223;0.053],[0.471;0.288;0.131;0.04]]
    a = [4;19;49;64;79;999]
     
    for i = 1:p.popsize 
        humans[i] = Human()              ## create an empty human       
        x = humans[i]
        x.idx = i 
        agn = rand(agedist)
        x.age = rand(agebraks[agn]) 
        x.ag = findfirst(y-> x.age in y, agebraks)
        
        g = findfirst(y->y>=x.age,a)
        
        
        x.ag_new = g
        
        if x.age >= 18
            g = findfirst(y->x.age ∈ y, age_break_behav)
            x.vac_behavior = sample([PRO; LIK; HES; ANT], Weights(prob_behav[g]))
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

        b1 = p.b
        b2 = rand(Distributions.Beta(1.5, 20)) #ba
        b3 = rand(Distributions.Beta(2, 25)) #bh
        b4 = rand(Distributions.Beta(1, 200)) #bl
        b5 = rand(Distributions.Beta(1.2, 20)) #be
        return [b1; b3; b4; b2; b5]
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
    l = findall(x -> x.health == SUS && x.ag == ag, humans)
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

function time_update()
    # counters to calculate incidence

    for x in humans 
        x.tis += 1 
        x.doi += 1 # increase day of infection. variable is garbage until person is latent
         
        x.daysinf += 1
        
        #TODO add behavioral change here
        if x.tis >= x.exp             
            @match Symbol(x.swap_status) begin
                :LAT  => begin 
                    move_to_latent(x); 
                end
                :PRE  => begin move_to_pre(x); end
                :ASYMP => begin move_to_asymp(x);  end
                :INF  => begin move_to_inf(x); end    
                :REC  => begin move_to_recovered(x); end
                :DED  => begin move_to_dead(x); end
                _    => begin dump(x); error("swap expired, but no swap set."); end
            end
        end

        # behavioral change
        update_behavior(x)

        get_nextday_counts(x)
        
    end

end
export time_update

function update_behavior(x)

    #? create a vector of Pv, Pa...
    #? 

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
        if rand() < (1-prob1)/((1-prob1)+(1-prob2*prob3))
            if rand() < 1-prob1 #? Is the individual changing?
                x.vac_behavior = PRO
            end
        else
            if rand() < 1-prob2*prob3
                x.vac_behavior = HES
            end
        end
    elseif x.vac_behavior == HES
          # probability of getting PRO
          n1 = x.contacts_vac[7]+p.α*x.contacts_vac[1] 
          prob1 = (1-x.betas_vac[1])^n1 # bs, bh, bl, ba, be
  
          # probability of getting anti
          n2 = x.contacts_vac[4]
          prob2 = (1-x.betas_vac[4])^n2
  
  
        if rand() < (1-prob1)/((1-prob1)+(1-prob2))
            if rand() < 1-prob1
                x.vac_behavior = PRO
            end
        else
            if rand() < 1-prob2
                x.vac_behavior = ANT
            end
        end 
    elseif x.vac_behavior == ANT
        # probability of getting likely
        n1 = x.contacts_vac[8] # number of vaccinated in contact
        prob1 = (1-x.betas_vac[3])^n1 # bs, bh, bl, ba, be


        if rand() < 1-prob1 # probability of getting at least one of the changes
            
            x.vac_behavior = LIK
            
        end
    end
    
end

#@inline _set_isolation(x::Human, iso) = _set_isolation(x, iso, x.isovia)
@inline function _set_isolation(x::Human, iso, via)
    # a helper setter function to not overwrite the isovia property. 
    # a person could be isolated in susceptible/latent phase through contact tracing
    # --> in which case it will follow through the natural history of disease 
    # --> if the person remains susceptible, then iso = off
    # a person could be isolated in presymptomatic phase through fpreiso
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # a person could be isolated in mild/severe phase through fmild, fsevere
    # --> if x.iso == true from CT and x.isovia == :ct, do not overwrite
    # --> if x.iso == true from PRE and x.isovia == :pi, do not overwrite
    if x.isovia == :null || via == :sev
        x.iso = iso 
        x.isovia = via
        x.daysisolation = 0
        x.days_after_detection = 0
    elseif !iso
        x.iso = iso 
        x.isovia = via
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
    
    x.got_inf = true
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

   # h.iso = true # a dead person is isolated
   # _set_isolation(h, true)  # do not set the isovia property here.  
    # isolation property has no effect in contact dynamics anyways (unless x == SUS)
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
   
    #h.daysinf = 999

    h.got_inf = false 
    
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
    else
        bf = 0.0
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
    
    if !x.iso 
        cnt = rand(negative_binomials(ag)) ##using the contact average for shelter-in
        
        x.nextday_meetcnt = cnt
    elseif !(x.health_status  in (DED))
        cnt = rand(negative_binomials_shelter(ag))  # expensive operation, try to optimize
        x.nextday_meetcnt_w = 0
        x.nextday_meetcnt = cnt
    end
    

    if x.health_status == DED
        x.nextday_meetcnt = 0
    end
   
    x.contacts_vac = [0;0;0;0;0;0;0;0;0]

    return cnt
end

function calculate_H(gp)
    
    status = [PRO; LIK; HES; ANT]
    
    vector = [humans[i].vac_behavior for i in gp]
    
    
    # how many times it appears
    contagens = countmap(vector)

    contagem_vetor = [get(contagens, i, 0) for i in status]
    # probabilities based on homophily
    
    Pj = zeros(Float64, length(status), length(status))
    if sum(contagem_vetor) > 0
        for br in status
            H = p.h .* Int.(br .== status)+(1-p.h) .* contagem_vetor./sum(contagem_vetor)
            Pj[Int(br)+1,:] = H ./ sum(H)
        end
    end
    return Pj
end

function dyntrans(sys_time, grps,sim)
    #totalmet = 0 # count the total number of contacts (total for day, for all INF contacts)
    #totalinf = 0 # count number of new infected 
    ## find all the people infectious
    #rng = MersenneTwister(246*sys_time*sim)

    Pj = map(x->calculate_H(x), grps)

    pos = shuffle(1:length(humans))

    # go through every infectious person
    for x in humans[pos]        
           
            xhealth = x.health
            cnts = x.nextday_meetcnt
            cnts == 0 && continue # skip person if no contacts
            #general population contact
            gpw = Int.(round.(cm[x.ag]*cnts))
            
            #cnts = number of contacts on that day

            
            perform_contacts(x,gpw,grps,xhealth, Pj)
    end
    #return totalmet, totalinf
end
export dyntrans

function return_contacts(x, gp, g, Pj)

    if g == 0
        aux = []
        return aux
    end

    if x.vac_behavior == UNDEFV
        # if the individual does not require behavior, we sample it from uniform
        aux1 = sample(gp, g)
        return aux1
    end


    status = [PRO; LIK; HES; ANT]
    
    vector = [humans[i].vac_behavior for i in gp]

    if any(vector .== UNDEFV)

        pos = findall(x-> humans[x].vac_behavior == UNDEFV, gp)

        g1 = Int(round(g*length(pos)/length(gp)))

        aux1 = sample(pos, g1)
        g = g-g1

        if g == 0
            return gp[aux1]
        end
    else
        aux1 = []
    end

    Pjj = Pj[Int(x.vac_behavior)+1, :]
    # sampling the groups based on H
    gr = sample(1:length(Pjj), Weights(Pjj), g, replace = true)
    
    vector2 = [status[i] for i in gr]
    
    # counting the number of contact in each group
    contagens = countmap(vector2)
    contagem_vetor = [get(contagens, i, 0) for i in status]
    aa = [findall(y -> y == status[i], vector) for i in eachindex(contagem_vetor)]
    aux = repeat([[]], length(aa))
    for i in eachindex(contagem_vetor)
        if length(aa[i]) > 0
            aux[i] = sample(aa[i], contagem_vetor[i], replace = true)
        end
    end

    aux = vcat(aux...)

    aux = [aux;aux1]

    return gp[aux]

end

function perform_contacts(x,gpw,grp_sample,xhealth, Pj)

    for (i, g) in enumerate(gpw) 
        meet = return_contacts(x, grp_sample[i], g, Pj[i])#rand(grp_sample[i], g)   # sample the people from each group
        # go through each person
        for j in meet 
            y = humans[j]

        
            ycnt = y.nextday_meetcnt  
            y.nextday_meetcnt = y.nextday_meetcnt - min(1,ycnt) # remove a contact

            ycnt == 0 && continue
            
            if x.vac_status*y.vac_status == 0 # only gets to it if it is necessary
                if x.vac_status == 0 && x.vac_behavior != UNDEF
                    y.contacts_vac[Int(x.vac_behavior)+1] += 1
                end
    
                if y.vac_status == 0 && y.vac_behavior != UNDEF
                    x.contacts_vac[Int(y.vac_behavior)+1] += 1
                end
    
                if x.vac_status == 1
                    if x.wentto == 1
                        if  x.health_status == REC
                            y.contacts_vac[5] += 1
                        elseif x.health_status == DED
                            y.contacts_vac[6] += 1
                        end
                    end
                    if x.vac_status == PRO
                        y.contacts_vac[7] += 1
                    end
                    y.contacts_vac[8] += 1
    
                elseif y.vac_status == 1
                    if y.wentto == 1
                        if  y.health_status == REC
                            x.contacts_vac[5] += 1
                        elseif y.health_status == DED
                            x.contacts_vac[6] += 1
                        end
                    end
                    if y.vac_status == PRO
                        x.contacts_vac[7] += 1
                    end
                    x.contacts_vac[8] += 1
                end
            end
            


            aux = 0

            if y.health == SUS && xhealth ∈ (PRE, INF, ASYMP) && y.swap == UNDEF
                aux = 1
                beta = _get_betavalue(xhealth)*(1-p.vaccine_eff)^y.vac_status
            elseif xhealth == SUS && y.health ∈ (PRE, INF, ASYMP) && y.swap == UNDEF
                aux = 2
                beta = _get_betavalue(y.health)*(1-p.vaccine_eff)^x.vac_status
            else
                beta = 0.0
            end

            if rand() < beta
                if aux == 1
                    y.exp = y.tis   ## force the move to latent in the next time step.
                    y.sickfrom = xhealth ## stores the infector's status to the infectee's sickfrom
                    y.sickby = y.sickby < 0 ? x.idx : y.sickby
                           
                    y.swap = LAT
                    y.swap_status = LAT
                    y.daysinf = 0
                    y.dur = sample_epi_durations(y)
                elseif aux == 2
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