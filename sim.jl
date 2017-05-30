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
function delta!(F0, F1)


end


# ------------------------------
# Specify the service parameters
# ------------------------------

dt = 1e-5
t = 0:dt:100
N = length(t)

Delta = 0.01;
Delta1 = Delta;
Delta2 = Delta;

# Deadline
Dmax = 0.02;


# Define costs
jq1 = 0.5;
jq2 = 0.5;
jc1 = 6.0;
jc2 = 8.0;

r = 17.0 # the input rate to the service-chain
s1 = 6.0
s2 = 8.0
sigma1 = r/s1
sigma2 = r/s2
m1 = floor(r/s1)
m2 = floor(r/s2)
rho1 = r/s1 - m1
rho2 = r/s2 - m2


F0 = NFV.VNF(1, N, dt, r, 0.0, 1.0, 0.0, 0.0, 0.0 )
F1 = NFV.VNF(1, N, dt, s1, rho1, sigma1, Delta, jq1, jc1)
F2 = NFV.VNF(1, N, dt, s2, rho2, sigma2, Delta, jq2, jc2)


qmax!(F0, F1)
qmax!(F1, F2)

F1.acost = F1.jq*F1.qmax
F2.acost = F2.jq*F2.qmax
a = F1.acost+F2.acost
@show a
@show F1.Tbar
@show F2.Tbar

F1.delta =  F1.rho*(F1.ns*(1-F1.rho))/r
F2.delta = F2.ns/r*F2.rho^2*((1-F2.rho)*F2.ns-(1-F1.rho)*F1.ns)/((1-F1.rho)*F1.ns + F2.rho*F2.ns)

c = Dmax/(F1.delta + F2.delta)

@show F1.delta
@show F2.delta
@show c


# -------------------------------
# Solve the optimization problem
# -------------------------------

# This has to be done automatically!!!

# Step 1:
# define all points in [0,c] in which the cost function is not differentiable


# Step 2:
# look at the interior point for all intervals where the derivative is zero

# Step 3:
# find the smallest point => yields the period.




# -------------------------------
# Simulate the system
# -------------------------------

# start by updating the parameters to include the period