using DataStructures


"""
NeuralPopulation
"""
mutable struct NeuralPopulation{T1,T2}
	name::String
	N::Int
	Cm::Float64
	taum::Float64
	rest::Float64
	reset::Float64
	threshold::Float64
	refractory::Float64
	spike_delay::Float64
	spike_delay_storage::Queue{Float64}
	receptors::Vector{T1}
	targets::OrderedDict{String,T2}
	LTP::Bool
	self_connection::Bool
end

NeuralPopulation(name, N=1,
				 Cm=0.0,
				 taum=0.0,
				 rest=0.0,
				 reset=0.0,
				 threshold=0.0;
				 refractory=0.0,
				 spike_delay=0.0,
				 spike_delay_storage=Queue{Float64}(),
				 receptors=Vector{Receptor}(),
				 targets=OrderedDict{String,Target}(),
				 LTP=false, self_connection=false) = NeuralPopulation(
								name, N, Cm, taum, rest, reset,
								threshold, refractory,
								spike_delay, spike_delay_storage,
								receptors, targets, LTP, self_connection)


"""
Receptor
"""
mutable struct Receptor
	name::String
	tau::Float64
	reversal::Float64
	external_freq::Float64
	mean_external_efficacy::Float64
	mean_external_connection::Float64
end


"""
Target
"""
mutable struct Target
	name::String
	target_receptor::String
	mean_efficacy::Float64
	weight::Float64
end


"""
Group
"""
mutable struct Group
	name::String
	members::Vector{String}
end


"""
Event
"""
abstract type AbstractEvent end

mutable struct EventMembraneNoise <: AbstractEvent
	time::Float64
	event_type::String
	population::String
	Gauss_mean::Float64
	Gauss_std::Float64
end

mutable struct EventExtFreq <: AbstractEvent
	time::Float64
	event_type::String
	population::String
	receptor::String
	external_freq::Float64
end

mutable struct EventSpike <: AbstractEvent
	time::Float64
	event_type::String
	population::String
	receptor::String
end

mutable struct EventEmitSpike <: AbstractEvent
	time::Float64
	event_type::String
	population::String
end

mutable struct EventEndTrial <: AbstractEvent
	time::Float64
	event_type::String
end


"""
Output
"""
abstract type AbstractOutput end

mutable struct OutputMemPot <: AbstractOutput
	filename::String
	output_type::String
	population::String
end

mutable struct OutputSpike <: AbstractOutput
	filename::String
	output_type::String
	population::String
end

mutable struct OutputFiringRate <: AbstractOutput
	filename::String
	output_type::String
	population::String
	firing_rate_indow::Int
	print_step::Int
end


"""
Network
"""
mutable struct Network
	neu::OrderedDict{String,NeuralPopulation}
	receptors::OrderedDict{String,Receptor}
	group::Vector{Group}
	event::Vector{AbstractEvent}
	output::Vector{AbstractOutput}
	event_temp::Vector{EventExtFreq}
	event_end::EventEndTrial

	#  Inner Constructor
	Network() = new(OrderedDict{String,NeuralPopulation}(),
					OrderedDict{String,Receptor}(),
					Vector{Group}(),
					Vector{AbstractEvent}(),
					Vector{AbstractOutput}(),
					Vector{EventExtFreq}(),
				    EventEndTrial(0.0, "EndTrial"))
end


#  --------------------------------------------------------


function add_neuron(Net::Network, name::String)
	neuron = NeuralPopulation(name)
	Net.neu[name] = neuron
end


function add_receptor(Net::Network,
					  name::String,
					  tau::Real,
					  reversal::Real,
					  external_freq::Real,
					  mean_external_efficacy::Real,
					  mean_external_connection::Real)

	rec = Receptor(name, tau, reversal, external_freq, mean_external_efficacy, mean_external_connection)
	Net.receptors[name] = rec
end


function set_neuron_param_all(Net::Network,
							  N::Int,
							  Cm::Real,
							  taum::Real,
							  rest::Real,
							  reset::Real,
							  threshold::Real)

	for neu in values(Net.neu)
		neu.N = N
		neu.Cm = Cm
		neu.taum = taum
		neu.rest = rest
		neu.reset = reset
		neu.threshold = threshold
	end
end


function set_neuron_receptor_all(Net::Network, args...)
	for neu in Net.neu, rec in args
		push!(neu.second.receptors, Net.receptors[rec])
	end
end


function add_synapse(Net::Network,
					pre_syn::String,
					post_syn::String,
					target_receptor::String,
					mean_efficacy::Float64,
				    weight::Float64)

	post_syn in keys(Net.neu) || begin @error "No target neural population can be found."; exit(1) end
	Net.neu[pre_syn].targets[post_syn] = Target(post_syn, target_receptor, mean_efficacy, weight)
end


function add_event(Net::Network, time::Real, event_type::String, args...)
	if event_type == "MembraneNoise"
		event = EventMembraneNoise(time, event_type, args...)
		push!(Net.event, event)

	elseif event_type == "Spike"
		event = EventSpike(time, event_type, args...)
		push!(Net.event, event)

	elseif event_type == "EmitSpike"
		event = EventEmitSpike(time, event_type, args...)
		push!(Net.event, event)

	elseif event_type == "ExtFreq"
		event = EventExtFreq(time, event_type, args...)
		push!(Net.event_temp, event)

	elseif event_type == "EndTrial"
		event = EventEndTrial(time, event_type, args...)
		Net.event_end = event
	end
end
