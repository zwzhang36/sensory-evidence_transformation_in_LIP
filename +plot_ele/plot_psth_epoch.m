function F = plot_psth_epoch(neuron, bin)

neuron.spike_group(bin,'align','fixation');
timeline = neuron.get_timeline(bin).mean;
condition= neuron.race.condition;
spCount  = neuron.spCount.race_fix.mean;
subWeight= neuron.subweight;


F = plot_ele.plot_psth_epoch_weight(spCount, timeline, condition, subWeight);
end