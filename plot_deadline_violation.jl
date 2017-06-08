#=

Simple script to plot the probability of the deadline_violation.

/Victor Millnert


Assumptions:
- the data that will be saved is "deadline_violation"
- format of dedaline_violation is MxN
- M::Int64 is the number of simulations
- N::Int64 is the length of each simulation
- feedback::Bool indicating whether feedback was used or not
- t::Vector{Float64, N} time-vector for the simulation

=#
using PyPlot

# 1. compute missed dedaline-probability
missed_deadline_prob = sum(deadline_violation, 1)[:]./M

title_string = feedback ? "prob. with feedback" : "prob. without feedback"

# PLOT THE STATESPACE TRAJECTORY OF QUEUE-1 AND QUEUE-2
figure("$(title_string)", figsize=(10,5))
PyPlot.plot(t, missed_deadline_prob)
xlabel("time (s)")
ylabel("probability")
title("$(title_string)")

