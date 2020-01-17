using DifferentialEquations
using ModelingToolkit
using DataStructures
using Random

include("elements.jl")


@parameters t
@derivatives D'~t


mutable struct NetworkIndex
	vars::OrderedDict{String,Operation}
	ode_idx::OrderedDict{String,Int}
	param_idx::OrderedDict{String,Int}
	refractory_recovery_time::OrderedDict{String,Float64}

	NetworkIndex() = new(OrderedDict{String,Operation}(),
						 OrderedDict{String,Int}(),
						 OrderedDict{String,Int}(),
						 OrderedDict{String,Float64}())
end

network_idx = NetworkIndex()


function gen_variables(net::Network)
	vars = OrderedDict{String,Operation}()

	param_index = OrderedDict{String,Int}()
	counter = 1

	for neu in values(net.neu)
		vars[neu.name] = Variable(Symbol(neu.name))(t)

		current = neu.name * "_I"
		vars[current] = Variable(Symbol(current))()
		param_index[current] = counter
		counter += 1

		refractory_status = neu.name * "_refractory_status"
		vars[refractory_status] = Variable(Symbol(refractory_status))()
		param_index[refractory_status] = counter
		counter += 1

		for rec in neu.receptors
			neu_rec = neu.name * '_' * rec.name
			vars[neu_rec] = Variable(Symbol(neu_rec))(t)
		end
	end

	network_idx.param_idx = param_index
	network_idx.vars = vars
end


function LIF(neu::NeuralPopulation)
	-1.0 / neu.taum * (network_idx.vars[neu.name] - neu.rest)
end


function receptor_model(neu::NeuralPopulation, rec::Receptor)
	neu_rec = neu.name * '_' * rec.name
	-network_idx.vars[neu_rec] * (network_idx.vars[neu.name] - rec.reversal) / 1000
end


function gen_ode(net::Network)

	ode_idx = OrderedDict{String,Int}()
	counter = 1

	eqs = []

	for neu in values(net.neu)

		# receptor
		receptor_current = Operation[]
		for rec in neu.receptors
			neu_rec = neu.name * '_' * rec.name
			eq = D(network_idx.vars[neu_rec]) ~ -network_idx.vars[neu_rec] / rec.tau

			push!(eqs, eq)
			ode_idx[neu_rec] = counter
			counter += 1

			push!(receptor_current, receptor_model(neu, rec))
		end

		# neuron
		current = neu.name * "_I"
		refractory_status = neu.name * "_refractory_status"
		eq = D(network_idx.vars[neu.name]) ~ network_idx.vars[refractory_status] *
											 (LIF(neu) + ( +(receptor_current..., network_idx.vars[current]) / neu.Cm ))

		push!(eqs, eq)
		ode_idx[neu.name] = counter
		counter += 1

	end

	network_idx.ode_idx = ode_idx
	return eqs
end


function callback_fire(neu::NeuralPopulation)
	condition(u, t, integrator) = (u[ network_idx.ode_idx[neu.name] ] >= neu.threshold)

	function affect!(integrator)
		integrator.u[ network_idx.ode_idx[neu.name] ] = neu.reset

		if integrator.p[ network_idx.param_idx[neu.name*"_refractory_status"] ] == 1.0
			integrator.p[ network_idx.param_idx[neu.name*"_refractory_status"] ] = 0.0
			network_idx.refractory_recovery_time[neu.name] = integrator.t + neu.refractory
		end

		if neu.spike_delay == 0.0
			for tar in values(neu.targets)
				integrator.u[ network_idx.ode_idx[tar.name*'_'*tar.target_receptor] ] += (tar.mean_efficacy * tar.weight)  # No Upper Bond
			end
		else
			enqueue!(neu.spike_delay_storage, integrator.t + neu.spike_delay)
		end
	end

	DiscreteCallback(condition, affect!, save_positions=(false,false))
end


function callback_refractory_detection(neu::NeuralPopulation)
	function condition(u, t, integrator)
		integrator.p[ network_idx.param_idx[neu.name*"_refractory_status"] ] == 0.0 &&
		integrator.t >= network_idx.refractory_recovery_time[neu.name]
	end

	affect!(integrator) = integrator.p[ network_idx.param_idx[neu.name*"_refractory_status"] ] = 1.0

	DiscreteCallback(condition, affect!, save_positions=(false,false))
end


function callback_spike_delay_detection(neu::NeuralPopulation)
	function condition(u, t, integrator)
		!isempty(neu.spike_delay_storage) &&
		integrator.t >= first(neu.spike_delay_storage)
	end

	function affect!(integrator)
		for tar in values(neu.targets)
			integrator.u[ network_idx.ode_idx[tar.name*'_'*tar.target_receptor] ] += (tar.mean_efficacy * tar.weight)  # No Upper Bond
		end
		dequeue!(neu.spike_delay_storage)
	end

	DiscreteCallback(condition, affect!, save_positions=(false,false))
end


function callback_event(event::EventSpike)
	condition(u, t, integrator) = (t == event.time)

	affect!(integrator) = (integrator.u[ network_idx.ode_idx[event.population*'_'*event.receptor] ] += 1.0)
	# affect!(integrator) =  # Set Upper Bond
	#     integrator.u[ network_idx.ode_idx[event.population*'_'*event.receptor] ] +=
	#         (1.0 - integrator.u[ network_idx.ode_idx[event.population*'_'*event.receptor] ])

	DiscreteCallback(condition, affect!, save_positions=(false,false))
end


function callback_event(event::EventEmitSpike)
	condition(u, t, integrator) =
		t == event.time && integrator.p[ network_idx.param_idx[event.population*"_refractory_status"] ] == 1.0

	affect!(integrator) = integrator.u[ network_idx.ode_idx[event.population] ] = 0.0

	DiscreteCallback(condition, affect!, save_positions=(false,false))
end


function callback_event(event::EventMembraneNoise)
	condition(u, t, integrator) = (t == event.time)

	affect!(integrator) =
		integrator.p[ network_idx.param_idx[event.population*"_I"] ] = (event.Gauss_mean + event.Gauss_std * randn())

	DiscreteCallback(condition, affect!, save_positions=(false,false))
end


function poisson_point_process(firing_rate; duration=1.0, dt=0.001)
	# unit >> (ms)
	firing_rate *= dt
	duration /= dt

	t_spikes = Float64[]
	time = randexp()

	while time < duration
		push!(t_spikes, time)
		time += randexp() / firing_rate
	end

	# unit >> (s)
	t_spikes * dt
end


function frequency_to_spike(net::Network)
	freq_temp = OrderedDict{Tuple{String,String},Vector{Tuple{Float64,Float64}}}()
	for e in net.event_temp
		if !haskey(freq_temp, (e.population, e.receptor))
			freq_temp[(e.population, e.receptor)] = []
		end
		push!(freq_temp[(e.population, e.receptor)], (e.time, e.external_freq))
	end

	for (neu_rec, t_freq) in freq_temp
		for i in 1:(length(t_freq)-1)
			if t_freq[i][2] != 0.0
				for t_spike in poisson_point_process(t_freq[i][2], duration=(t_freq[i+1][1]-t_freq[i][1]))
					add_event(net, t_freq[i][1]+t_spike, "Spike", neu_rec[1], neu_rec[2])
				end
			end
		end

		# last component
		if t_freq[end][2] != 0.0
			for t_spike in poisson_point_process(t_freq[end][2], duration=((net.event_end.time)-t_freq[end][1]))
				add_event(net, t_freq[end][1]+t_spike, "Spike", neu_rec[1], neu_rec[2])
			end
		end
	end

	sort!(net.event, by=x->x.time)
end


function initialize(var::String, net::Network)
	var in keys(net.neu) ? net.neu[var].rest : 0.0
end


function gen_tstops(net::Network)
	unique([ 0.0:0.1:net.event_end.time ; [event.time for event in net.event] ]) |> sort  # 0.1 is time step
end


function gen_problem(net::Network)
	# ODE equations
	gen_variables(net)
	eqs = gen_ode(net)
	de = ODESystem(eqs)
	f = ODEFunction(de,
					[network_idx.vars[ode] for ode in keys(network_idx.ode_idx)],
					[network_idx.vars[p] for p in keys(network_idx.param_idx)])

	net.event_end.time != 0.0 || begin @error "You haven't set the end time yet!"; exit(1) end
	frequency_to_spike(net)

	# Callback Functions
	refractory_callbacks = [callback_refractory_detection(neu) for neu in values(net.neu)]
	spike_delay_callbacks = [callback_spike_delay_detection(neu) for neu in values(net.neu)]
	event_callbacks = [callback_event(event) for event in net.event]
	fire_callbacks = [callback_fire(neu) for neu in values(net.neu)]
	cb = CallbackSet(refractory_callbacks..., spike_delay_callbacks..., event_callbacks..., fire_callbacks...)

	u0 = [initialize(ode, net) for ode in keys(network_idx.ode_idx)]
	tspan = (0.0, net.event_end.time)
	p = repeat([0.0, 1.0], length(net.neu))

	prob = ODEProblem(f, u0, tspan, p, callback=cb)
	save_idxs = [network_idx.ode_idx[neu.name] for neu in values(net.neu)]
	(prob, gen_tstops(net), save_idxs)

	# # MethodError
	# sol = solve(prob, tstops=tstops)
	# sol
end
