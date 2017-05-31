#=

Simualtion of the forwarding graphs in

"Feedback for increased robustness of forwarding graphs in the Cloud"
by Victor Millnert, Johan Eker, and Enrico Bini.


The simulation allows for random generation of a forwarding graph of
length n. It the solves for the optimal schedule for that forwarding
graph and ends by simulating if with and without stochastic input
traffic.

=#

include("./NFV.jl")
using PyPlot

# ------------------------------------------------------------
#          SPECIFY VARIOUS HELP FUNCTIONS
# ------------------------------------------------------------


#=
Computes the maximum queue-size for function F1 given function F0
=#
function qmax!(F0, F1)

    F1.qmax = max(max(F1.rho*(F1.ns*(1-F1.rho) - F0.ns*(1-F0.rho)),
            (1-F0.rho)*(F0.ns*F0.rho - F1.ns*F1.rho)),
            max(F0.rho*(F0.ns*(1-F0.rho) - F1.ns*(1-F1.rho)),
            (1-F1.rho)*(F1.ns*F1.rho - F0.ns*F0.rho)));
end


#=
Computes the extra latency that function F1 add to the chain given F0
=#
function delta!(r, F0, F1)

    if (F1.m0+1)*F1.ns >= (F0.m0+1)*F0.ns && F1.m0*F1.ns >= F0.m0*F0.ns
        # Case 1a

        F1.delta = 1/r*F1.ns*F1.rho^2*(F1.ns*(1-F1.rho) -
        F0.ns*(1-F0.rho)) / (F0.ns*(1-F0.rho) + F1.ns*F1.rho)

    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns >= F0.m0*F0.ns
        # Case 1b
        F1.delta = 0.0

    elseif (F1.m0+1)*F1.ns >= (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2a

        if F1.Ton >= F0.Ton

            F1.delta = 1/r*(F0.ns*F0.rho*(F0.rho-F1.rho) + (1-F1.rho)*(F1.ns*F1.rho - F0.ns*F0.rho))

        else

            F1.delta = 1/r*(F1.rho*(F1.ns*(1-F1.rho) - F0.ns*(1-F0.rho)) + F0.ns*(F0.rho-1)*(F0.rho-F1.rho))
            
        end
    
    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2b

        F1.delta = 1/r*F1.ns*(1-F1.rho)^2*(F1.ns*F1.rho + F0.ns*F0.rho) / (F1.ns*(1-F1.rho) + F0.ns*F0.rho)
        
    end
end

#=
Computes t_on for function F1 given function F0
=#
function ton!(F0, F1)

    if (F1.m0+1)*F1.ns >= (F0.m0+1)*F0.ns && F1.m0*F1.ns >= F0.m0*F0.ns
        # Case 1a

        F1.ton = F0.ton + F1.T*F1.rho*(F1.ns*(1-F1.rho) -
        F0.ns*(1-F0.rho)) / (F0.ns*(1-F0.rho) + F1.ns*F1.rho)

    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns >= F0.m0*F0.ns
        # Case 1b
        F1.ton = F0.ton

    elseif (F1.m0+1)*F1.ns >= (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2a

        F1.ton = F0.ton + F1.T*(F0.rho - F1.rho)

    
    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2b

        # F1.ton = F0.ton + F1.T*(1-F1.rho)

        F1.ton = F0.ton - F1.T*(1-F1.rho)*(F1.ns*F1.rho -
        F0.ns*F0.rho)/(F1.ns*(1-F1.ns) + F0.ns*F0.rho)

    end
end

#=
Computes the extra latency that function F1 add to the chain given F0
=#
function qon!(F0, F1)
    
    if (F1.m0+1)*F1.ns >= (F0.m0+1)*F0.ns && F1.m0*F1.ns >= F0.m0*F0.ns
        # Case 1a
        F1.case = "case 1a"
        F1.qon = F1.T*F1.rho*(F1.ns*(1-F1.rho) - F0.ns*(1-F0.rho))

    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns >= F0.m0*F0.ns
        # Case 1b
        F1.qon = 0.0

        F1.case = "case 1b"

    elseif (F1.m0+1)*F1.ns >= (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2a

        if F1.Ton >= F0.Ton

            F1.qon = F1.T*(1 - F1.rho)*(F1.ns*F1.rho - F0.ns*F0.rho)
        F1.case = "case 2a1"

        else

            F1.qon = F1.T*F1.rho*(F1.ns*(1-F1.rho) - F0.ns*(1-F0.rho))
            F1.case = "case 2a2"

        end
    
    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2b
        F1.case = "case 2b"

        F1.qon = F1.T*(1-F1.rho)*(F1.ns*F1.rho - F0.ns*F0.rho)
        
    end
end


"""

Computes the cost for a period T given a, lower bound on the cost Jlb
and a specific forwarding graph.

"""
function cost(T, a, Jlb, FG)
    J = a*T + Jlb

    for F in FG
        if T < F.Tbar
            J = J + F.jc*(1-F.rho)
        elseif T >= F.Tbar
            J = J + F.jc*F.Delta/T
        end
        
    end
    
    return J
    
end


# ------------------------------------------------------------
#          GENERATE THE FORWARDING GRAPHS
# ------------------------------------------------------------

# if the input should be stochastic
stochastic = false

# if the functions should use feedback to compute the necessary
# on-time
feedback = true

dt = 1e-5
t = 0:dt:5
N = length(t)

m = 2 # number of functions in the chain

# Deadline
Dmax = 0.05;

Delta_min = 0.01 # minumum time-overhead
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
r_min = 10.0 
r_max = 50.0
r_d   = 0.1
r = rand(r_min:r_d:r_max)

# define the service capacity for one VM
ns_min = 5.0
ns_max = 10.0
ns_d   = 0.01
ns = rand(ns_min:ns_d:ns_max, m) 


# Parameters used for the example (uncomment to run)
# r = 17
# Delta = [0.01 0.01]
# jq = [0.5 0.5]
# jc = [6.0 8.0]
# ns = [6.0 8.0]


sigma = r./ns
m0 = floor(Int,r./ns)
rho = r./ns - m0

# initialize the "dummy function"
F0 = NFV.VNF(1, N, dt, r, 0.0, 1.0, 0.0, 0.0, 0.0 )

FG = []

for i = 1:m
    append!(FG, [NFV.VNF(i, N, dt, ns[i], rho[i], sigma[i], Delta[i], jq[i], jc[i])])
end

F0.m0 = 1

# compute qmax and delta
for i in 1:length(FG)
    FG[i].m0 = m0[i]
    if i == 1
        qmax!(F0, FG[i])
        delta!(r, F0, FG[i])
    else
        qmax!(FG[i-1], FG[i])
        delta!(r, FG[i-1], FG[i])

    end

    FG[i].acost = FG[i].jq*FG[i].qmax
    
end


a = sum(map(F-> F.acost, FG))

c = Dmax/sum(map(F-> F.delta, FG))

# @show a
# @show F1.Tbar
# @show F2.Tbar
# @show F1.delta
# @show F2.delta
# @show c



# ------------------------------------------------------------
#             SOLVE THE OPTIMIZATION PROBLEM
# ------------------------------------------------------------

# Create a new vector and reorder it based on the values of Tbar
FGs = sort(FG, by=x->x.Tbar)

# Compute the lower bound of the cost
Jlb = sum(map(x-> x.jc*r/x.ns, FGs))

# pick out the functions having less critical period than the largest possible period
FGf = filter(x-> x.Tbar < c, FGs)

# Step 1:
# define all points in [0,c] in which the cost function is not differentiable
C = [0 ; map(x-> x.Tbar, FGf) ; c]

# Step 2:
# look at the interior point for all intervals where the derivative is zero
Cstar = []

for i in 1:length(FGf)-1
    cstar = sqrt(sum(map(x-> x.jc*x.Delta, FGf[1:i]))/a)
    # check if cstar is withing the interval!
    if FGf[i].Tbar < cstar && cstar < FGf[i+1].Tbar
        append!(Cstar, cstar)
    end
end

# Step 3:
# find the smallest point => yields the period.

# investigate all the points in C and Cstar and compute the cost
# choose the period yielding the smallest cost!

append!(C, Cstar)
J = zeros(length(C))
for i in 1:length(C)
    J[i] = cost(C[i], a, Jlb, FGs)
end

(Jmin,Tstar_i) = findmin(J)
Tstar = C[Tstar_i]
# @show Jmin, Tstar


# Now we have to update the values (i.e., qmax, Ton, Toff, ...) of the
# functions in the chain to also include the period
for F in FG
    F.T = Tstar
    F.qmax = F.qmax*Tstar
    F.Ton = Tstar*F.rho
    F.Toff = Tstar*(1-F.rho)    
end


###################################

# Compute the start-times ton and queue-size at on-time
FG[1].ton = FG[1].qmax/(r - FG[1].m0*FG[1].ns)
FG[1].qon = FG[1].qmax

if length(FG) > 1
    for i in 2:length(FG)
        ton!(FG[i-1], FG[i])
        qon!(FG[i-1], FG[i])
    end
end

# @show FG[1].ton, FG[1].qon
# @show FG[2].ton, FG[2].qon

###################################



# ------------------------------------------------------------
#                     SIMULATE THE SYSTEM
# ------------------------------------------------------------

# start by updating the parameters to include the period
Tk = max(1,floor(Int,Tstar/dt))
    
for F in FG
    F.tonk = floor(Int, F.ton/dt)
end

Dmaxk = floor(Int, Dmax/dt)
Delay_violation = zeros(Bool, N)


for i=1:N
    
    k = floor(Int, i/Tk) # find what period we're in
    
    for F in FG

        # add input to chain if its the first function in the chain
        if F.id == 1
            if stochastic
                # stochastic input rate
                F.a[i] = r*dt*0.9 + 0.2*r*dt*rand()
            else
                # deterministic input rate
                F.a[i] = r*dt
            end
        end

        # compute cumulative arrivals
        if i == 1
            F.A[i] = F.a[i]
        else
            F.A[i] = F.A[i-1] + F.a[i]
        end

        # update the queue
        if i == 1
            F.q[i] = F.a[i]
        else
            F.q[i] = F.q[i-1] + F.a[i]
        end
        
        # compute the on-time (in samples) based on the queue-size from i-1
        if i == k*Tk + F.tonk 
            extra_ontime = 0.0
            if i > 1 && feedback
                # we should use feedback to compute necessary on-time
                extra_ontime = (F.q[i] - F.qon)/F.ns
            # else
            #     # Don't use feedback...
            #     extra_ontime = 0.0
            end
            F.Ton_tilde = F.Ton + extra_ontime
            F.steps_left = floor(Int, F.Ton_tilde/dt)
        end


        # need special case for the initialization of case 2b
        if i == 1 && F.case == "case 2b"
            extra_ontime = 0.0
            F.Ton_tilde = F.Ton + extra_ontime
            F.steps_left = floor(Int, F.Ton_tilde/dt)
        end
        # compute the service to the queue
        if F.steps_left > 1
            # the additional machine is on
            F.s[i] = (F.m0+1)*F.ns*dt
            F.steps_left = F.steps_left-1
        else
            # the additional machine is off
            F.s[i] = F.m0*F.ns*dt
        end

        if i == 1
            F.S[i] = F.s[i]
        else
            F.S[i] = F.S[i-1] + F.s[i]
        end

      

        # compute the departures
        F.d[i] = min(F.q[i], F.s[i])

        if i == 1
            F.D[i] = F.d[i]
        else
            F.D[i] = F.D[i-1] + F.d[i]
        end
                
        # update the queue-size
        F.q[i] = F.A[i] - F.D[i]

        # add the departures to the next function in the chain
        if F.id < length(FG)
            FG[F.id+1].a[i] = F.d[i]
        end

    end

    # Compute if the deadline is violated or not
    if i > Dmaxk && FG[end].D[i] < FG[1].A[i-Dmaxk]
        
       Delay_violation[i] = 1
    end

end


# ------------------------------------------------------------
#                     PLOT VARIOUS METRICS
# ------------------------------------------------------------
close("all")
figure("Statespace of queue trajectory", figsize=(10,5))
PyPlot.plot(FG[1].q, FG[2].q,  label="T(t)" )
xlabel("queue 1")
ylabel("queue 2")
title("statespace of queue-trajectory")

figure("Queue size 1", figsize=(10,5))
PyPlot.plot(t, FG[1].q,  label="q(t)" )
xlabel("time (s)")
ylabel("q(t) ")
title("queue size for 1st function")

figure("Queue size 2", figsize=(10,5))
PyPlot.plot(t, FG[2].q,  label="q(t)" )
xlabel("time (s)")
ylabel("q(t) ")
title("queue size for 2nd function")



# figure("End-to-end delay for chain", figsize=(10,5))
# # PyPlot(t, Delay, label="Delay")
# PyPlot.plot(t, Delay_violation, label="violation")
# xlabel("time (s)")
# ylabel("delay (s)")
# title("end-to-end delay")
