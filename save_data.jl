#=

Simple script to write the current simulation data to a text-file that
will later be used to generate the TiKZ-plot.

/Victor Millnert


Assumptions:
- the data that will be saved is "deadline_violation"
- format of dedaline_violation is MxN
- M::Int64 is the number of simulations
- N::Int64 is the length of each simulation
- feedback::Bool indicating whether feedback was used or not
- t::Vector{Float64, N} time-vector for the simulation

=#

# 1. compute missed dedaline-probability
missed_deadline_prob = sum(deadline_violation, 1)[:]./M

# 2. write it to a text-file
filename = feedback ? "feedback.txt" : "no_feedback.txt"
f = open("$(filename)", "w")
for i in 1:N
    write(f, "$(t[i]) \t $(missed_deadline_prob[i]) \n")
end
close(f)
