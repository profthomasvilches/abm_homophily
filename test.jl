

sim = 1
ip = ModelParameters(prov = :usa, h = 0.0001, vaccination_rate = 10, β = 0.1)


# without vaccine
# who is initially infected
# adding vaccine on top

# L-P; L-H; AH ;HA ;LHR; PL; HL

aa = [1.0;1.85;0.8;1.0;2.57;2.6;0.8]#
h = 0.0
run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 1, 2000, aa) #divided by 10 and 2
run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 1, 2000, aa) #divided by 10 and 2
run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 2, 2000, aa) #divided by 10 and 2
run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 2, 2000, aa) #divided by 10 and 2
run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 3, 2000, aa) #divided by 10 and 2
run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 3, 2000, aa) #divided by 10 and 2
run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 4, 2000, aa) #divided by 10 and 2
run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 4, 2000, aa) #divided by 10 and 2


for h in 0.15:0.05:1.0
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 1, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 1, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 2, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 2, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 3, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 3, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 4, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 4, 2000, aa) #divided by 10 and 2
end



for h in 0.02:0.04:0.1
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 1, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 1, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 2, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 2, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 3, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 3, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 4, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 5, 2000, aa) #divided by 10 and 2
end


for h in 0.002:0.004:0.01
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 1, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 1, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 2, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 2, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 3, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 3, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 4, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 4, 2000, aa) #divided by 10 and 2
end


for h in 0.0001:0.0003:0.001
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 1, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 1, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 2, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 2, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 3, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 3, 2000, aa) #divided by 10 and 2
    run_param(1, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, 4, 2000, aa) #divided by 10 and 2
    run_param(6, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, 4, 2000, aa) #divided by 10 and 2
end


aa = [1.0;6.0;1.0;1.0;1.0;6.0;1.0]
h = 0.0
run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(3, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2


for h in 0.15:0.05:1.0
    run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(3, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end



for h in 0.02:0.01:0.1
    run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(3, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end


for h in 0.002:0.001:0.01
    run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(3, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end


for h in 0.0001:0.0001:0.001
    run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(3, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end




aa = [1.0;1.0;10.0;1.0;1.0;1.0;5.0]
h = 0.0
run_param(4, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(5, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2


for h in 0.15:0.05:1.0
    run_param(4, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end



for h in 0.02:0.01:0.1
    run_param(4, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end


for h in 0.002:0.001:0.01
    run_param(4, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end


for h in 0.0001:0.0001:0.001
    run_param(4, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
    run_param(5, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
end




aa = [1.0;1.4;1.0;1.0;2.;2.;1.0]
h = 0.01
#run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(8, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

h = 0.0
#run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(8, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

h = 1.0
#run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(8, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2



# aa = [1.0;1.8;1.0;1.0;2.5;2.5;1.0]
# aa = [1.0;2.0;1.0;1.0;2.7;2.7;1.0]
# aa = [1.0;1.85;1.0;1.0;2.55;2.55;1.0]

# L-P; L-H; AH ;HA ;LHR; PL; HL
#aa = [1.0;1.87;1.0;1.0;2.57;2.57;1.0]
aa = [1.0;1.85;0.8;1.0;2.57;2.6;0.8]#

h = 0.01
#run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(9, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

h = 0.0
#run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(9, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

h = 1.0
#run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(9, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

h =  0.5
#run_param(2, 0.1, h, 0.05, :prob, 0, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2
run_param(9, 0.1, h, 0.05, :prob, 325, 0.7, 0.7, nothing, nothing, 2000, aa) #divided by 10 and 2

