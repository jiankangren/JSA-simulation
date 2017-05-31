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
    m::Vector{Int64} # number of running instances
    m0::Int64 # number of instances always kept on
    
    # Delay
    delay::Vector{Float64} # delay for passing through the function

    # Constants used for solving the optimization problem
    Tbar::Float64 # limit for when there is no more time to switch off the machine
    qmax::Float64 # maximum queue-size
    acost::Float64
    delta::Float64 # extra latency that this function adds to the chain
    
    
    # Variables used for the feedback-law
    T::Float64 # Optimal period for the schedule
    Ton::Float64 # time the extra machine should be on
    Ton_tilde::Float64 # the altered on-time for the additional machine
    Toff::Float64 # time the extra machine should be off
    ton::Float64  # absolute time when the machine should be switched on
    toff::Float64 # absolute time when the machine should be switched off
    
    qon::Float64 # desired queue-size when the machine is switched on
    steps_left::Int64 # simple counter of how many simulation steps
                      # more the additional machine should be switched
                      # on for
    
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
    zeros(N), # delay
    Delta/(1-rho), # Tbar
    Inf, # qmax
    Inf, # acost
    Inf, # delta
    0, # T
    0, # Ton
    0, # Ton_tilde
    0, # Toff
    0, # ton
    0, # toff
    0, # qon
    0) # steps_left
    
    

end

end
