using Plots

include("elements.jl")
include("simulator.jl")


Net = Network()


# (name, N, Cm, taum, rest, reset, threshold)
add_neuron(Net, "neux")
add_neuron(Net, "neuy")
set_neuron_param_all(Net, 1, 0.4, 14.0, -70.0, -55.0, -50.0)


# (name, tau, reversal, external_freq, mean_external_efficacy, mean_external_connection)
add_receptor(Net, "Ach", 20.0, 0.0, 0.0, 2.1, 1.0)
add_receptor(Net, "AMPA", 2.0, 0.0, 0.0, 2.1, 1.0)
add_receptor(Net, "NMDA", 100.0, 0.0, 0.0, 2.1, 1.0)
add_receptor(Net, "GABA", 5.0, -90.0, 0.0, 0.0, 0.0)
set_neuron_receptor_all(Net, "AMPA", "GABA")


# (Pre_syn, Post_syn, target_receptor, mean_efficacy, weight)
add_synapse(Net, "neux", "neuy", "AMPA", 2.0, 1.0)


# (time, "MembraneNoise", population, Gauss_mean, Gauss_std)
# (time, "ExtFreq", population, receptor, external_freq)
# (time, "Spike", population, receptor)
# (time, "EmitSpike", population)
# (time, "EndTrial")
add_event(Net, 0.1, "MembraneNoise", "neux", 1.0, 0.0)
add_event(Net, 100.0, "EndTrial")


prob, tstops, save_idxs = gen_problem(Net)
sol = solve(prob, tstops=tstops, saveat=0.1, save_idxs=save_idxs, abstol=1e-9, reltol=1e-9)

p = plot(sol, plotdensity=10000, leg=false, dpi=200)
savefig(p, "MemPot.png")
