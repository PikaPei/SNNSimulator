Usage
=====

Main Features
--------------

#### Modules
``` julia
include("elements.jl")
include("simulator.jl")
```

#### Create Network
``` julia
Net = Network()
```

#### Add Neuron and Set Parameters
``` julia
add_neuron(Net, "neux")
add_neuron(Net, "neuy")

# (N, Cm, taum, rest, reset, threshold)
set_neuron_param_all(Net, 1, 0.4, 14.0, -70.0, -55.0, -50.0)
```

#### Add Receptors and Attach To Neuron
``` julia
# (name, tau, reversal, external_freq, mean_external_efficacy, mean_external_connection)
# External frequency is for background input. (Not yet implemented)
add_receptor(Net, "AMPA", 2.0, 0.0, 0.0, 0.0, 0.0)
add_receptor(Net, "GABA", 5.0, -90.0, 0.0, 0.0, 0.0)
set_neuron_receptor_all(Net, "AMPA", "GABA")
```

#### Add Synapse
``` julia
# (Pre_syn, Post_syn, target_receptor, mean_efficacy, weight)
add_synapse(Net, "neux", "neuy", "AMPA", 1.0, 1.0)
```

#### Add Event
+ MembraneNoise
+ ExtFreq
+ Spike
+ EmitSpike
+ EndTrial _(sole & necessary)_

``` julia
# (time, "MembraneNoise", population, Gauss_mean, Gauss_std)  # std hasn't been implemented well
# (time, "ExtFreq", population, receptor, external_freq)
# (time, "Spike", population, receptor)
# (time, "EmitSpike", population)
# (time, "EndTrial")
# "ALL" is a build-in group. (Not yet implemented)
add_event(Net, 0.1, "MembraneNoise", "ALL", 1.0, 0.1)  # Time cannot start from 0.0
add_event(Net, 20.0, "ExtFreq", "neux", "AMPA", 50.0)
add_event(Net, 100.0, "EndTrial")
```

#### Add Output Format (Not yet implemented)
+ Membrane Potential: "MemPot" or 'M'
+ Spike: "Spike" or 'S'
+ Firing Rate: "FiringRate" or 'F'

``` julia
# (filename, "MemPot", population)
# (filename, "Spike", population)
# (filename, "FiringRate", population, window, step)
# "ALL" is a build-in group. (Not yet implemented)
add_output(Net, "MemPot.dat", "MemPot", "ALL")
add_output(Net, "Spike.dat", "Spike", "ALL")
add_output(Net, "FiringRate.dat", "FiringRate", "neux", 100, 10)
```

#### Add Group (Not yet implemented)
``` julia
# (group_name, *members)
add_group(Net, "neu_group", "neux", "neuy")
```


Run Simulation
--------------

#### Solve ODEProblem
``` julia
prob, tstops, save_idxs = gen_problem(Net)
sol = solve(prob, tstops=tstops, saveat=0.1, save_idxs=save_idxs)
# Currently, timestep can only be 0.1
```

#### Plot Membrane Potential
```julia
using Plots
p = plot(sol, plotdensity=10000, leg=false, dpi=200)
savefig(p, "MemPot.png")
```
