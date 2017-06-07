#=

Simualtion of the forwarding graphs in

"Feedback for increased robustness of forwarding graphs in the Cloud"
by Victor Millnert, Johan Eker, and Enrico Bini.


The simulation allows for random generation of a forwarding graph of
length n. It the solves for the optimal schedule for that forwarding
graph and ends by simulating if with and without stochastic input
traffic.

=#

@everywhere include("./NFV.jl")
using PyPlot

    
# ------------------------------------------------------------
#          GENERATE THE FORWARDING GRAPHS
# ------------------------------------------------------------

m = 3 # number of functions in the chain 

recreate_example = true # set to true if you wish to use the
                         # example-chain in the paper

randomness = 0.4  # the amount of randomness in the stochastic input
stochastic = true # if the input should be stochastic
feedback = true  # if the functions should use feedback or not

dt = 1e-5 # sample time
tend = 20 # simulation end-time
t = 0:dt:tend
N = length(t)

deadline_violation = zeros(1,N) # vector to compute the probability of a
                              # missed deadlines

# ------------------------------------------------------------
#          GENERATE THE FORWARDING GRAPHS
# ------------------------------------------------------------
println("initializing forwarding graph")

# Deadline
Dmax = 0.10;

Delta_min = Dmax/(m*10) # minumum time-overhead
Delta_max = Dmax/m # maximum time-overhead
Delta_d   = 0.001 # steps for the difference
Delta = rand(Delta_min:Delta_d:Delta_max, m)



# Define costs
# buffer cost
jq_min = 0.1
jq_max = 1.0
jq_d   = 0.01
jq = rand(jq_min:jq_d:jq_max, m)

# compute cost
jc_min = 2.0
jc_max = 20.0
jc_d   = 0.01
jc = rand(jc_min:jc_d:jc_max, m)

# define input traffic rate to the function
r_min = 20.0 
r_max = 50.0
r_d   = 0.1
r = rand(r_min:r_d:r_max)

# define the service capacity for one VM
ns_min = 5.0
ns_max = 10.0
ns_d   = 0.01
ns = rand(ns_min:ns_d:ns_max, m) 


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Parameters used for the example (uncomment to run)  #!
if recreate_example
    m = 2                                                 #!
    r = 17.0                                              #!
    Dmax = 0.02                                           #!
    Delta = [0.01 0.01]                                   #!
    jq = [0.5 0.5]                                        #!
    jc = [6.0 8.0]                                        #!
    ns = [6.0 8.0]                                        #!
end
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


sigma = r./ns
m0 = floor(Int,r./ns)
rho = r./ns - m0

# initialize the "dummy function"
F0 = NFV.VNF(1, N, dt, r, 0.0, 1.0, 0.0, 0.0, 0.0 )
F0.m0 = 1

# create the forwarding graph
FG = []
for i = 1:m
    append!(FG, [NFV.VNF(i, N, dt, ns[i], rho[i], sigma[i], Delta[i], jq[i], jc[i])])
end

# compute qmax and delta 
for i in 1:length(FG)
    FG[i].m0 = m0[i]
    if i == 1
        NFV.qmax!(F0, FG[i])
        NFV.delta!(r, F0, FG[i])
    else
        NFV.qmax!(FG[i-1], FG[i])
        NFV.delta!(r, FG[i-1], FG[i])

    end

    FG[i].acost = FG[i].jq*FG[i].qmax
    
end

# compute a and c
a = sum(map(F-> F.acost, FG))
c = Dmax/sum(map(F-> F.delta, FG))

# ------------------------------------------------------------
#             SOLVE THE OPTIMIZATION PROBLEM
# ------------------------------------------------------------

println("solving optimization problem")

# Create a new vector and reorder it based on the values of Tbar
FGs = sort(FG, by=x->x.Tbar)

# Compute the lower bound of the cost
Jlb = sum(map(x-> x.jc*r/x.ns, FGs))

# Pick out the functions having less critical period than the
# largest possible period
FGf = filter(x-> x.Tbar < c, FGs)

# Define all points in [0,c] in which the cost function is not differentiable
C = [0 ; map(x-> x.Tbar, FGf) ; c]

# Look at the interior point for all intervals where the derivative is zero
Cstar = []
for i in 1:length(FGf)-1
    cstar = sqrt(sum(map(x-> x.jc*x.Delta, FGf[1:i]))/a)
    # check if cstar is withing the interval!
    if FGf[i].Tbar < cstar && cstar < FGf[i+1].Tbar
        append!(Cstar, cstar)
    end
end

# Find the smallest point => yields the period.

# Investigate all the points in C and Cstar and compute the cost
# choose the period yielding the smallest cost!

append!(C, Cstar)
J = zeros(length(C))
for i in 1:length(C)
    J[i] = NFV.cost(C[i], a, Jlb, FGs)
end

(Jmin,Tstar_i) = findmin(J)
Tstar = C[Tstar_i] # the optimal period

# Update the values (i.e., qmax, Ton, Toff, ...) of the functions
# in the chain to also include the period
for F in FG
    F.T = Tstar
    F.qmax = F.qmax*Tstar
    F.Ton = Tstar*F.rho
    F.Toff = Tstar*(1-F.rho)    
end


# ------------------------------------------------------------
#                GENERATE THE DESIRED SCHEDULE
# ------------------------------------------------------------


# Compute the start-times ton and queue-size at on-time
if FG[1].rho < 1e-5
    FG[1].ton = tend
elseif FG[1].rho > (1-1e-5)
    FG[1].ton = dt
else
    FG[1].ton = FG[1].qmax/(r - FG[1].m0*FG[1].ns)
end

NFV.qon!(F0,FG[1])

if length(FG) > 1
    for i in 2:length(FG)
        NFV.ton!(FG[i-1], FG[i])
        NFV.qon!(FG[i-1], FG[i])
    end
end


Tk = max(1,floor(Int,Tstar/dt)) # discretize the period 

for j in 1:length(FG)
    try
        FG[j].tonk = floor(Int, FG[j].ton/dt)
    catch e
        @show e
    end

    if FG[j].ton < 0.0
        FG[j].ton = FG[j].ton + FG[j].T
        FG[j].tonk = floor(Int, FG[j].ton/dt)
        FG[j].ton_changed = true

    end
    
    if FG[j].ton > FG[j].T
        FG[j].ton = FG[j].ton - FG[j].T
        FG[j].tonk = floor(Int, FG[j].ton/dt)
        FG[j].ton_changed = true
    end
end


# ------------------------------------------------------------
#                     SIMULATE THE SYSTEM
# ------------------------------------------------------------

Dmaxk = floor(Int, Dmax/dt) # discretize the deadline

# add input to chain if its the first function in the chain
if stochastic
    # stochastic input rate
    FG[1].a = r*dt*(1-randomness/2) + randomness*r*dt*rand(N)
else
    # deterministic input rate
    FG[1].a = r*dt*ones(N)
end

sim = 1
M = 1
@time NFV.sim_fg!(FG, r, Tk, Dmaxk, dt, N, feedback, stochastic, M, sim, deadline_violation)


# ------------------------------------------------------------
#                     PLOT VARIOUS METRICS
# ------------------------------------------------------------
close("all")


# PLOT THE STATESPACE TRAJECTORY OF QUEUE-1 AND QUEUE-2
figure("Statespace of queue trajectory", figsize=(10,5))
PyPlot.plot(FG[1].q, FG[2].q,  label="T(t)" )
xlabel("queue 1")
ylabel("queue 2")
title("statespace of queue-trajectory")


# PLOT QUEUE-SIZED FOR EVERY FUNCTION IN THE CHAIN
for F in FG
    figure("Queue size #$(F.id)", figsize=(10,5))
    PyPlot.plot(t, F.q,  label="q_$(F.id)(t)" )
    xlabel("time (s)")
    ylabel("q(t) ")
    title("queue size #$(F.id)")
end
