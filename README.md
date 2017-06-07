# JSA-simulation
Simulation for the Journal of Systems Architecture


## Recreating the state-space trajectory

To recreate the example chain and the state-space trajectory, use
`sim.jl`. Make sure to have `recreate_example = true`. You can then
choose whether to use feedback or not, and whether to have a
stochastic input or no.

You can also try out the randomly generated service chains.

To run a file you need to start `Julia`[https://julialang.org]. Then
run a script as below (which will run `sim.jl`)

```
julia> include("sim.jl")
```

## Recreating missed deadline probabilities

To recreate the missed deadline probabilities you need to run
`sim_feedback.jl` (for functions using feedback) or
`sim_nofeedback.jl` (to have the functions not use feedback).  Finally
one can use `plot_dedaline_violation.jl` after all the simulations has
stopped in order to plot the result.


