module NFV

"""

    VNF(id::Int, N::Int, dt::Float, ns::Float64, rho::Float64,
                  sigma::Float64, Delta::Float64, jq::Float64, jc::Float64)

Special type for representing a virtual network function

    id    - the functions place in the chain (i.e., id=i is the i'th function)
    N     - number of simulation points
    dt    - time-step for the simulation
    ns    - service capacity for one VM
    Delta - time-overhead for realizeing the service-rate (i.e. time
            for starting/stopping a new instance)

# Example
```julia

julia> VNF(1, 10000, 0.001, 0.31, 0.112, 0.13)

```
    
"""
type VNF

    id::Int64 # id of the function, i.e. what place it is in the
              # service-chain
    
    N::Int64 # size of the vectors
    dt::Float64 # time-step used for the simulation

    # Cost parameters
    jq::Float64 # cost for queue allocation
    jc::Float64 # cost for computation
    
    # Input 
    a::Vector{Float64} # arrival-rate
    A::Vector{Float64} # cumulative arrivals

    # Buffer and queue related
    q::Vector{Float64} # queue-size
    
    # Service
    s::Vector{Float64} # service rate
    S::Vector{Float64} # cumulative amount of served packets

    # Departure
    d::Vector{Float64} # departure rate
    D::Vector{Float64} # cumulative amount of departured packets
    
    # Instances
    Delta::Float64 # time-overhead for changing the service-rate

    # Service parameters
    ns::Float64 # nominal service rate for one instance
    rho::Float64 # normalized residual rate
    sigma::Float64 # exact number of machines needed
    m::Vector{Int64} # 1 if additional machine is running
    m0::Int64 # number of instances always kept on
    ton_changed::Bool # boolean whether the on-time was changed or not..
    
    # # Delay
    # delay::Vector{Float64} # delay for passing through the function

    # Constants used for solving the optimization problem
    Tbar::Float64 # limit for when there is no more time to switch off the machine
    qmax::Float64 # maximum queue-size
    acost::Float64
    delta::Float64 # extra latency that this function adds to the chain
    
    
    # Variables used for the feedback-law
    T::Float64 # Optimal period for the schedule
    Tk::Int64  # Optimal period (expressed in floor(T/dt)
    Ton::Float64 # time the extra machine should be on
    Ton_tilde::Float64 # the altered on-time for the additional machine
    Toff::Float64 # time the extra machine should be off
    ton::Float64  # absolute time when the machine should be switched on
    tonk::Int64   # absolute time (in number of samples) that the machine should be switched on
    toff::Float64 # absolute time when the machine should be switched off
    
    qon::Float64 # desired queue-size when the machine is switched on
    steps_left::Vector{Int64} # simple counter of how many simulation steps
                      # more the additional machine should be switched
    # on for

    case::String # what case this function belongs to
    
    VNF(id,
        N,
        dt,
        ns,
        rho,
        sigma,
        Delta,
        jq,
        jc) =
            new(id,
                N,
                dt,
                jq, # jq
                jc, # jc
                zeros(N), # a
                zeros(N), # A
                zeros(N), # q
                zeros(N), # s
                zeros(N), # S
                zeros(N), # d
                zeros(N), # D
                Delta,    # Delta
                ns,       # ns
                rho,      # rho
                sigma,    # sigma
    zeros(Int, N), # m
    0, # m0
    false, #ton_changed
    # zeros(N), # delay
    Delta/(1-rho), # Tbar
    Inf, # qmax
    Inf, # acost
    Inf, # delta
    0.0, # T
    0.0, # Tk
    0.0, # Ton
    0.0, # Ton_tilde
    0.0, # Toff
    0.0, # ton
    0, # tonk
    0.0, # toff
    0.0, # qon
    zeros(Int, N),   # steps_left
    "")  # case 
    
    

end




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
        
        F1.delta = 1/r*F1.ns*F1.rho^2*(F1.ns*(1-F1.rho) - F0.ns*(1-F0.rho)) / (F0.ns*(1-F0.rho) + F1.ns*F1.rho)

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
        
        F1.ton = F0.ton + F1.T*F1.rho*(F1.ns*(1-F1.rho) - F0.ns*(1-F0.rho)) / (F0.ns*(1-F0.rho) + F1.ns*F1.rho)

    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns >= F0.m0*F0.ns
        # Case 1b
        F1.ton = F0.ton

    elseif (F1.m0+1)*F1.ns >= (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2a
        F1.ton = F0.ton + F1.T*(F0.rho - F1.rho)

    elseif (F1.m0+1)*F1.ns < (F0.m0+1)*F0.ns && F1.m0*F1.ns < F0.m0*F0.ns
        # Case 2b
        F1.ton = F0.ton - F1.T*(1-F1.rho)*(F1.ns*F1.rho - F0.ns*F0.rho)/(F1.ns*(1-F1.rho) + F0.ns*F0.rho)
        
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


#=
Computes the cost for a period T given a, lower bound on the cost Jlb
and a specific forwarding graph.
=#
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



function sim_fg!(FG, r::Float64, Tk::Int64, Dmaxk::Int64, dt::Float64, N::Int64, feedback::Bool, stochastic::Bool, M::Int64, sim::Int64, deadline_violation)

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
                    println("-------------")
                    @show i
                    @show Tk
                    @show FG[j].id
                    @show FG[j].q[i]
                    @show FG[j].qon
                    @show extra_ontime
                end
                
                # FG[j].Ton_tilde = FG[j].Ton + extra_ontime
                FG[j].Ton_tilde = FG[j].Ton + max(extra_ontime,0.0)
                FG[j].steps_left[i] = floor(Int, FG[j].Ton_tilde/dt)
            end

            # We have to initialize the system in a good state some
            # functions should be initialized with their additional
            # machine on
            if i == 1
                if FG[j].ton + FG[j].Ton > FG[j].T
                    FG[j].Ton_tilde = FG[j].ton + FG[j].Ton - FG[j].T
                    FG[j].steps_left[i] = floor(Int, FG[j].Ton_tilde/dt)
                end                
            end

            # compute the service to the queue
            if FG[j].steps_left[i] > 0
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
            deadline_violation[sim,i] =  1.0
        end

    end # end simulating FG
    
end #sim_fg!()





end


