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

        F1.ton = F0.ton + F1.T*(1-F1.rho)
        
    end
end

#=
Computes the extra latency that function F1 add to the chain given F0
=#
function qon!(F0, F1)

    if (F1.m0+1)*F1.ns >= (F0.m0+1)*F0.ns && F1.m0*F1.ns >= F0.m0*F0.ns
        # Case 1a

        F1.qon = F1.T*F1.rho*(F1.ns*(1-F1.rho) - F0.ns*(1-F0.rho))

    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns >= F0.m0*F0.ns
        # Case 1b
        F1.qon = 0.0

    elseif (F1.m0+1)*F1.ns >= (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2a

        if F1.Ton >= F0.Ton

            F1.qon = F1.T*(1 - F1.rho)*(F1.ns*F1.rho - F0.ns*F0.rho)

        else

            F1.qon = F1.T*F1.rho*(F1.ns*(1-F1.rho) - F0.ns*(1-F0.rho))
            
        end
    
    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2b

        F1.qon = F1.T*(1-F1.rho)*(F1.ns*F1.rho + F0.ns*F0.rho)
        
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

dt = 1e-5
t = 0:dt:1
N = length(t)

m = 5 # number of functions in the chain

Delta_min = 0.01 # minumum time-overhead
Delta_max = 0.10 # maximum time-overhead
Delta_d   = 0.0001 # steps for the difference
Delta = rand(Delta_min:Delta_d:Delta_max, m)

# Deadline
Dmax = 0.02;


# Define costs
# buffer cost
jq_min = 0.1
jq_max = 1.0
jq_d   = 0.001
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
