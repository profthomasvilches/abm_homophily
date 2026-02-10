

# sim = 1
#ip = ModelParameters(prov = :usa, h = 0.9, vaccination_rate = 10, β = 0.1, function_homophily = :created)
# # xi = 10.0

# for xi in 5.0:1.0:15.0
# kapa = 1.0
# aa = [1.0;1.85;0.8;1.0;2.57;2.6;0.8]#
# h = 0.0
# run_param(Int(xi), 0.1, h, kapa, 200,  nothing, 2000, aa, xi)
# h = 0.5
# run_param(Int(xi), 0.1, h, kapa, 200,  nothing, 2000, aa, xi)
# h = 0.55
# run_param(Int(xi), 0.1, h, kapa, 200,  nothing, 2000, aa, xi)
# h = 1.0
# run_param(Int(xi), 0.1, h, kapa, 200,  nothing, 2000, aa, xi)
# end
# without vaccine
# who is initially infected
# adding vaccine on top


### Power



xi = 5.0
# L-P; L-H; AH ;HA ;LHR; PL; HL
kapa = 1.0


ww = 0.05
aa = [0.5,ww*1.0,ww*1.0,0.5]#
h = 0.0
run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :power) #divided by 10 and 2
run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :power) #divided by 10 and 2

###### aqui

for h in 0.15:0.05:1.0
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :power) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :power) #divided by 10 and 2
end



for h in 0.02:0.04:0.1
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :power) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :power) #divided by 10 and 2
end

for h in 0.002:0.004:0.01
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :power) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :power) #divided by 10 and 2
end


for h in 0.0005:0.0005:0.001
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :power) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :power) #divided by 10 and 2
end


### :assortative

xi = 5.0
# L-P; L-H; AH ;HA ;LHR; PL; HL
kapa = 1.0


ww = 0.05
aa = [0.5,ww*1.0,ww*1.0,0.5]#
h = 0.0
run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :assortative) #divided by 10 and 2
run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :assortative) #divided by 10 and 2

###### aqui

for h in 0.15:0.05:1.0
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :assortative) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :assortative) #divided by 10 and 2
end



for h in 0.02:0.04:0.1
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :assortative) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :assortative) #divided by 10 and 2
end

for h in 0.002:0.004:0.01
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :assortative) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :assortative) #divided by 10 and 2
end


for h in 0.0005:0.0005:0.001
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :assortative) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :assortative) #divided by 10 and 2
end


### :created



h = 0.0
run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :created) #divided by 10 and 2
run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :created) #divided by 10 and 2

###### aqui

for h in 0.15:0.05:1.0
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :created) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :created) #divided by 10 and 2
end



for h in 0.02:0.04:0.1
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :created) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :created) #divided by 10 and 2
end

for h in 0.002:0.004:0.01
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :created) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :created) #divided by 10 and 2
end


for h in 0.0005:0.0005:0.001
    run_param(1, 0.1, h, kapa, 0,  nothing, 2000, aa, xi, :created) #divided by 10 and 2
    run_param(6, 0.1, h, kapa, 200,  nothing, 2000, aa, xi, :created) #divided by 10 and 2
end