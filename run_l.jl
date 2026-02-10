
### Power

for l in 0.0:0.2:0.8

    xi = 5.0
    # L-P; L-H; AH ;HA ;LHR; PL; HL
    kapa = 1.0

    ww = 0.05
    aa = [0.5,ww*1.0,ww*1.0,0.5]#

    LV = [l, 2.0, 0.2]

    h = 0.0
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :power, LV, true) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :power, LV, true) #divided by 10 and 2

    ###### aqui

    for h in 0.05:0.05:1.0
        run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :power, LV, true) #divided by 10 and 2
        run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :power, LV, true) #divided by 10 and 2
    end



    ### :assortative

    h = 0.0
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :assortative, LV, true) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :assortative, LV, true) #divided by 10 and 2

    ###### aqui

    for h in 0.05:0.05:1.0
        run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :assortative, LV, true) #divided by 10 and 2
        run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :assortative, LV, true) #divided by 10 and 2
    end




    ### :created



    h = 0.0
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :created, LV, true) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :created, LV, true) #divided by 10 and 2

    ###### aqui

    for h in 0.05:0.05:1.0
        run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :created, LV, true) #divided by 10 and 2
        run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :created, LV, true) #divided by 10 and 2
    end

end
