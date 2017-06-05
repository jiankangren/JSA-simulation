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

    # F1.qmax = max(max(F1.rho*(F1.ns*(1-F1.rho) - F0.ns*(1-F0.rho)),
    #                   (1-F0.rho)*(F0.ns*F0.rho - F1.ns*F1.rho)),
    #               max(F0.rho*(F0.ns*(1-F0.rho) - F1.ns*(1-F1.rho)),
    #                   (1-F1.rho)*(F1.ns*F1.rho - F0.ns*F0.rho)));

     if (F1.m0+1)*F1.ns >= (F0.m0+1)*F0.ns && F1.m0*F1.ns >= F0.m0*F0.ns
        # Case 1a

         F1.qmax = max(F1.rho*(F1.ns*(1-F1.rho) - F0.ns*(1-F0.rho)),
                       (1-F0.rho)*(F0.ns*F0.rho - F1.ns*F1.rho))

                       
    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns >= F0.m0*F0.ns
        # Case 1b

         F1.qmax = max(F0.rho*(F0.ns*(1-F0.rho) - F1.ns*(1-F1.rho)),
                       (1-F0.rho)*(F0.ns*F0.rho - F1.ns*F1.rho))
         
    elseif (F1.m0+1)*F1.ns >= (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2a

        if F1.Ton >= F0.Ton

            F1.qmax = (1-F1.rho)*(F1.ns*F1.rho - F0.ns*F0.rho)
        else

            F1.qmax = F1.rho*(F1.ns*(1-F1.rho) - F0.ns*(1-F0.rho))
            
        end
    
    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2b

         F1.qmax = max(F0.rho*(F0.ns*(1-F0.rho) - F1.ns*(1-F1.rho)),
                       (1-F1.rho)*(F1.ns*F1.rho - F0.ns*F0.rho))

     end
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

        F1.delta = 1/r*F1.ns*(1-F1.rho)^2*(F1.ns*F1.rho - F0.ns*F0.rho)/(F1.ns*(1-F1.rho) + F0.ns*F0.rho)
        
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

        F1.ton = F0.ton - F1.T*(1-F1.rho)*(F1.ns*F1.rho - F0.ns*F0.rho)/(F1.ns*(1-F1.rho) + F0.ns*F0.rho)
        # println("-----------")
        # @show F1.id
        # @show F0.rho
        # @show F0.ton
        # @show F0.ns
        # @show F1.rho
        # @show F1.ns
        # @show F1.T
        # @show F1.ton
        
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



function sim_fg!(FG, r::Float64, Tk::Int64, Dmaxk::Int64, dt::Float64, N::Int64, feedback::Bool, stochastic::Bool, M::Int64, sim::Int64, deadline_violation, debug)

    if debug
        deadline_violation = zeros(N)
    end
    for i in 1:N
        
        k = floor(Int, i/Tk) # find what period we're in
        
        for j in 1:length(FG)


            # compute cumulative arrivals
            if i == 1
                FG[j].A[i] = FG[j].a[i]
            else
                FG[j].A[i] = FG[j].A[i-1] + FG[j].a[i]
            end

            # update the queue
            if i == 1
                FG[j].q[i] = FG[j].a[i]
            else
                FG[j].q[i] = FG[j].q[i-1] + FG[j].a[i]
            end
            
            # compute the on-time (in samples) based on the queue-size from i-1
            if i == k*Tk + FG[j].tonk 
                extra_ontime = 0.0
                if i > 1 && feedback
                    # we should use feedback to compute necessary on-time
                    extra_ontime = (FG[j].q[i] - FG[j].qon)/FG[j].ns
                end
                
                FG[j].Ton_tilde = FG[j].Ton + max(extra_ontime,0.0)
                FG[j].steps_left[i] = floor(Int, FG[j].Ton_tilde/dt)
            end


            # if i == 1 && FG[j].case == "case 2b"
            #     # FG[j].Ton_tilde = FG[j].Ton
            #     FG[j].Ton_tilde = FG[FG[j].id-1].ton - FG[j].Toff+2*dt
            #     FG[j].steps_left[i] = floor(Int, FG[j].Ton_tilde/dt)
            # end

            # We have to initialize the system in a good state some
            # functions should be initialized with their additional
            # machine on
            if i == 1
                if FG[j].ton + FG[j].Ton > FG[j].T
                    FG[j].Ton_tilde = FG[j].ton + FG[j].Ton - FG[j].T
                    FG[j].steps_left[i] = floor(Int, FG[j].Ton_tilde/dt)
                end
                
                # FG[j].Ton_tilde = FG[j].Ton
                # FG[j].steps_left[i] = floor(Int, FG[j].Ton_tilde/dt)
            end

            # compute the service to the queue
            if FG[j].steps_left[i] > 1
                # the additional machine is on
                FG[j].s[i] = (FG[j].m0+1)*FG[j].ns*dt
                if i < N
                    FG[j].steps_left[i+1] = FG[j].steps_left[i]-1
                end
                FG[j].m[i] = 1
            else
                # the additional machine is off
                FG[j].s[i] = FG[j].m0*FG[j].ns*dt
            end

            if i == 1
                FG[j].S[i] = FG[j].s[i]
            else
                FG[j].S[i] = FG[j].S[i-1] + FG[j].s[i]
            end

            # compute the departures
            FG[j].d[i] = min(FG[j].q[i], FG[j].s[i])

            if i == 1
                FG[j].D[i] = FG[j].d[i]
            else
                FG[j].D[i] = FG[j].D[i-1] + FG[j].d[i]
            end
            
            # update the queue-size
            FG[j].q[i] = FG[j].A[i] - FG[j].D[i]

            # add the departures to the next function in the chain
            if FG[j].id < length(FG)
                FG[FG[j].id+1].a[i] = FG[j].d[i]
            end

        end # end j in 1:length(FG)

        # Compute if the deadline is violated or not
        if i > Dmaxk+2 && FG[end].D[i] < FG[1].A[i-Dmaxk-2]
            deadline_violation[i] =  deadline_violation[i] + 1.0/M
        end

    end # end simulating FG


    # --------------------------------------------
    #    FINAL CHECK FOR CATCHING THAT GODDAMN BUG
    # --------------------------------------------

    if debug
    @show sum(deadline_violation)/N

    
    if sum(deadline_violation)/N > 0.5
        close("all")
        for F in FG

            if F.id==1
                figure("Queue size #$(F.id)", figsize=(10,5))
                PyPlot.plot(t, F.q,  label="q_$(F.id)(t)" )
                PyPlot.plot(t, F.m,  label="m_$(F.id) is on?")
                xlabel("time (s)")
                ylabel("q(t) ")
                title("queue size #$(F.id)")

            else
                figure("Queue size #$(F.id)", figsize=(10,5))
                PyPlot.plot(t, F.q,  label="q_$(F.id)(t)" )
                PyPlot.plot(t, F.m,  label="m_$(F.id) is on?")
                PyPlot.plot(t, FG[F.id-1].m, label="m_$(F.id-1) is on?")
                xlabel("time (s)")
                ylabel("q(t) ")
                title("queue size #$(F.id)")
            end
            # figure("service #$(F.id)", figsize=(10,5))
            # PyPlot.plot(t, F.A,  label="A_$(F.id)(t)" )
            # PyPlot.plot(t, F.D,  label="A_$(F.id)(t)" )
            # xlabel("time (s)")
            # ylabel("s(t) ")
            # title("service #$(F.id)")

            println("---------------------")
            @show F.id
            @show F.qmax
            @show maximum(F.q)
            @show r
            @show F.case
            @show F.ns
            @show F.m0
            @show F.qon
            @show F.ton
            @show F.ton_changed
            @show F.T
            @show F.Ton
            @show F.Toff
            @show F.tonk
            @show F.rho
            @show r - F.m0*F.ns
            println("---------------------")
            
        end

        figure("Deadline plots", figsize=(10,5))
        PyPlot.plot(t, FG[1].A,  label="A_1(t)" )
        PyPlot.plot(t, FG[end].D,  label="D_n(t)" )
        PyPlot.plot(t[Dmaxk:end-1], FG[1].A[1:end-Dmaxk], label="deadline")
        xlabel("time (s)")
        ylabel("s(t) ")
        title("Deadline plots")

        figure("Deadline violation probability", figsize=(10,5))
        PyPlot.plot(t, deadline_violation, label="feedback")
        xlabel("time (s)")
        ylabel("violation probability")
        title("Deadline violation probability")
        error("break!!!")

    end

    end
    # --------------------------------------------
    #                END OF THE CHECK
    # --------------------------------------------

    
end #sim_fg!()



    
# ------------------------------------------------------------
#          GENERATE THE FORWARDING GRAPHS
# ------------------------------------------------------------

debug = false

M = 100 # number of simulations to be run

randomness = 0.4
# if the input should be stochastic
stochastic = true

# if the functions should use feedback to compute the necessary
# on-time
feedback = false

dt = 1e-3
tend = 300
t = 0:dt:tend
N = length(t)

deadline_violation = zeros(N)
input = zeros(N)

for sim = 1:M

    @show sim
    # ------------------------------------------------------------
    #          GENERATE THE FORWARDING GRAPHS
    # ------------------------------------------------------------
    # println("initializing forwarding graph")
    m = 6 # number of functions in the chain

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

    # println("solving optimization problem")
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
    if FG[1].rho < 1e-5
        FG[1].ton = tend
    elseif FG[1].rho > (1-1e-5)
        FG[1].ton = dt
    else
        FG[1].ton = FG[1].qmax/(r - FG[1].m0*FG[1].ns)
    end

    # if r - FG[1].m0*FG[1].ns == 0.0
    #     FG[1].ton = tend
    # else
    #     FG[1].ton = FG[1].qmax/(r - FG[1].m0*FG[1].ns)
    # end
    # # FG[1].ton = Inf
    # FG[1].qon = FG[1].qmax

    
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

fuck = false
    # for F in FG
    for j in 1:length(FG)
        try
            FG[j].tonk = floor(Int, FG[j].ton/dt)
        catch e
            println("---------------------")
            @show e
            @show sim
            @show FG[j].id
            @show r
            @show FG[j].case
            @show FG[j].ns, FG[j].m0
            @show FG[j].qon
            @show FG[j].ton
            @show FG[j].T, FG[j].Ton, FG[j].Toff
            @show FG[j].rho
            @show r - FG[j].m0*FG[j].ns
            # FG[j].ton = tend
            # FG[j].tonk = floor(Int, FG[j].ton/dt)
            println("---------------------")
        end

        if FG[j].ton < 0.0
            println("---------------------")
            println("t_on < 0")
            @show FG[j].id
            # @show r
            # @show FG[j].case
            # @show FG[j].ns, FG[j].m0
            # @show FG[j].qon
            # @show FG[j].ton
            # @show FG[j].T, FG[j].Ton, FG[j].Toff
            # @show FG[j].rho
            # @show r - FG[j].m0*FG[j].ns
            # # FG[j].ton = tend
            # println("---------------------")
            FG[j].ton = FG[j].ton + FG[j].T
            FG[j].tonk = floor(Int, FG[j].ton/dt)
            FG[j].ton_changed = true
            @show FG[j].ton_changed

        end
        
        if FG[j].ton > FG[j].T
            # fuck = true
            println("---------------------")
            println("t_on is greater than T")
            @show FG[j].id
            # @show r
            # @show FG[j].case
            # @show FG[j].ns, FG[j].m0
            # @show FG[j].qon
            # @show FG[j].ton
            # @show FG[j].T, FG[j].Ton, FG[j].Toff
            # @show FG[j].rho
            # @show r - FG[j].m0*FG[j].ns
            # # FG[j].ton = tend
            # println("---------------------")
            FG[j].ton = FG[j].ton - FG[j].T
            FG[j].tonk = floor(Int, FG[j].ton/dt)
            FG[j].ton_changed = true
            @show FG[j].ton_changed
        end
    end

# if !fuck
#     error("we want errorsssss...")
# end

Dmaxk = floor(Int, Dmax/dt)
# deadline_violation = zeros(Bool, N)


# add input to chain if its the first function in the chain
if stochastic
    # stochastic input rate
    FG[1].a = r*dt*(1-randomness/2) + randomness*r*dt*rand(N)
else
    # deterministic input rate
    FG[1].a = r*dt*ones(N)
end
input = input + FG[1].a./M

# println("starting simulation")

sim = 1

@time sim_fg!(FG, r, Tk, Dmaxk, dt, N, feedback, stochastic, M, sim, deadline_violation, debug)


end # end multiple sim


# ------------------------------------------------------------
#                     PLOT VARIOUS METRICS
# ------------------------------------------------------------
# close("all")

# figure("Statespace of queue trajectory", figsize=(10,5))
# PyPlot.plot(FG[1].q, FG[2].q,  label="T(t)" )
# xlabel("queue 1")
# ylabel("queue 2")
# title("statespace of queue-trajectory")

# if maximum(deadline_violation)==1
#     for F in FG

#         figure("Queue size #$(F.id)", figsize=(10,5))
#         PyPlot.plot(t, F.q,  label="q_$(F.id)(t)" )
#         xlabel("time (s)")
#         ylabel("q(t) ")
#         title("queue size #$(F.id)")

#         # figure("service #$(F.id)", figsize=(10,5))
#         # PyPlot.plot(t, F.s,  label="s_$(F.id)(t)" )
#         # xlabel("time (s)")
#         # ylabel("s(t) ")
#         # title("service #$(F.id)")

#         println("---------------")
#         @show F.id
#         @show F.case
#         @show F.qmax
#         @show maximum(F.q)
#     end
# end

# figure("on-time", figsize=(10,5))
# PyPlot.plot(t, FG[1].s,  label="s1(t)" )
# PyPlot.plot(t, FG[2].s,  label="s2(t)" )
# xlabel("time (s)")
# ylabel("s(t) ")
# title("Service rate for the functions")

# figure("steps_left", figsize=(10,5))
# PyPlot.plot(t, FG[2].steps_left,  label="steps_left" )
# xlabel("time (s)")
# ylabel("steps ")
# title("steps left")

# figure("Deadline violation", figsize=(10,5))
# PyPlot.plot(t, deadline_violation, color="red", marker="*", linestyle="nothing")
# xlabel("time (s)")
# ylabel("violation")
# title("Deadline violation")



# @show FG[2].case



figure("Deadline violation probability", figsize=(10,5))
PyPlot.plot(t, deadline_violation, label="feedback")
xlabel("time (s)")
ylabel("violation probability")
title("Deadline violation probability")
