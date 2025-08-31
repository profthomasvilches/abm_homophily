

# sim = 1
# ip = ModelParameters(prov = :usa, h = 0.0001, vaccination_rate = 10, β = 0.1)


# without vaccine
# who is initially infected
# adding vaccine on top

# L-P; L-H; AH ;HA ;LHR; PL; HL
kapa = 1.0
aa = [1.0;1.85;0.8;1.0;2.57;2.6;0.8]#
h = 0.0
run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
# run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(1), 2000, aa) #divided by 10 and 2
# run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(1), 2000, aa) #divided by 10 and 2
# run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(2), 2000, aa) #divided by 10 and 2
# run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(2), 2000, aa) #divided by 10 and 2
# run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(3), 2000, aa) #divided by 10 and 2
# run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(3), 2000, aa) #divided by 10 and 2
# run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(4), 2000, aa) #divided by 10 and 2
# run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(4), 2000, aa) #divided by 10 and 2

###### aqui

for h in 0.15:0.05:1.0
    run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(1), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(1), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(2), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(2), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(3), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(3), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(4), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(4), 2000, aa) #divided by 10 and 2
end



for h in 0.02:0.04:0.1
    run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(1), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(1), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(2), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(2), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(3), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(3), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(4), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(4), 2000, aa) #divided by 10 and 2
end

for h in 0.002:0.004:0.01
    run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(1), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(1), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(2), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(2), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(3), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(3), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(4), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(4), 2000, aa) #divided by 10 and 2
end


for h in 0.0001:0.0003:0.001
    run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(1), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(1), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(2), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(2), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(3), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(3), 2000, aa) #divided by 10 and 2
    # run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, Int8(4), 2000, aa) #divided by 10 and 2
    # run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, Int8(4), 2000, aa) #divided by 10 and 2
end




#? testing 0 p->l

aa = [1.0;1.85;0.8;1.0;2.57;0.0;0.8]#
#aa = [1.0;1.85;0.8;1.0;2.57;2.6;0.8]#
h = 0.0
#run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

###### aqui

for h in 0.15:0.05:1.0
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end



for h in 0.02:0.04:0.1
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end

for h in 0.002:0.004:0.01
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
   
end


for h in 0.0001:0.0003:0.001
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end

#? testing p->l very high
aa = [1.0;1.85;0.8;1.0;2.57;5.2;0.8]#
h = 0.0
#run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

###### aqui

for h in 0.15:0.05:1.0
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end



for h in 0.02:0.04:0.1
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end

for h in 0.002:0.004:0.01
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
   
end


for h in 0.0001:0.0003:0.001
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end



#! CHANGING KAPPA

# L-P; L-H; AH ;HA ;LHR; PL; HL
kapa = 0.0
aa = [1.0;1.85;0.8;1.0;2.57;2.6;0.8]#
h = 0.0
run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
# 

###### aqui

for h in 0.15:0.05:1.0
    run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end



for h in 0.02:0.04:0.1
    run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end

for h in 0.002:0.004:0.01
    run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end


for h in 0.0001:0.0003:0.001
    run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end




#? testing 0 p->l

aa = [1.0;1.85;0.8;1.0;2.57;0.0;0.8]#
#aa = [1.0;1.85;0.8;1.0;2.57;2.6;0.8]#
h = 0.0
#run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

###### aqui

for h in 0.15:0.05:1.0
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end



for h in 0.02:0.04:0.1
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end

for h in 0.002:0.004:0.01
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
   
end


for h in 0.0001:0.0003:0.001
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end

#? testing p->l very high
aa = [1.0;1.85;0.8;1.0;2.57;5.2;0.8]#
h = 0.0
#run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

###### aqui

for h in 0.15:0.05:1.0
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end



for h in 0.02:0.04:0.1
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end

for h in 0.002:0.004:0.01
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
   
end


for h in 0.0001:0.0003:0.001
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end



#! CHANGING KAPPA

# L-P; L-H; AH ;HA ;LHR; PL; HL
kapa = 2.0
aa = [1.0;1.85;0.8;1.0;2.57;2.6;0.8]#
h = 0.0
#run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
# 

###### aqui

for h in 0.15:0.05:1.0
    #run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end



for h in 0.02:0.04:0.1
    #run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end

for h in 0.002:0.004:0.01
    #run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end


for h in 0.0001:0.0003:0.001
    #run_param(1, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end




#? testing 0 p->l

aa = [1.0;1.85;0.8;1.0;2.57;0.0;0.8]#
#aa = [1.0;1.85;0.8;1.0;2.57;2.6;0.8]#
h = 0.0
#run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

###### aqui

for h in 0.15:0.05:1.0
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end



for h in 0.02:0.04:0.1
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end

for h in 0.002:0.004:0.01
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
   
end


for h in 0.0001:0.0003:0.001
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(4, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end

#? testing p->l very high
aa = [1.0;1.85;0.8;1.0;2.57;5.2;0.8]#
h = 0.0
#run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

###### aqui

for h in 0.15:0.05:1.0
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end



for h in 0.02:0.04:0.1
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end

for h in 0.002:0.004:0.01
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
   
end


for h in 0.0001:0.0003:0.001
    #run_param(2, 0.1, h, kapa, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, kapa, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    
end