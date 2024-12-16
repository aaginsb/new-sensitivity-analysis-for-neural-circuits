function [t_vec, smooth_m, t_vec_base1, t_vec_stim, t_vec_base2, spike_mat] = make_smooth_avg_fr(stim_info, t_fin)
% early drafts of this code was written by Kevin Hannay, and later drafts were written by Alex Ginsberg

% parse out the information from stimInfo:
    baseline_freq = stim_info(1); % frequency of incoming inputs in absence of stimulus
    stim_freq = stim_info(2);  % frequency of incoming inputs in the presence of stimulus from the specified type of input
    number_fibers = stim_info(3); % number of fibers over which we 
    duration_of_baseline = stim_info(4); %time over which we simulate incoming input before stimulus arrives
    duration_of_stimulus = stim_info(5);

%change to per ms (because duration is in ms):
    baseline_freq = baseline_freq/1000;
    stim_freq = stim_freq/1000;

% simulate the baseline period using poisson_spike_gen(fr, tSim, nTrials)
number_realizations = 1;
summing_firing_rates = 0;
for realization = 1 : number_realizations
    % simulate the incoming spike train before the stimulus starts
    [spike_mat_base1, t_vec_base1] = poisson_spike_gen(baseline_freq, duration_of_baseline, number_fibers);

    % simulate the stimulus period--t_vec_stim gives stimulus
    [spike_mat_stim, t_vec_stim] = poisson_spike_gen(stim_freq, duration_of_stimulus, number_fibers);
    t_vec_stim = (t_vec_stim + duration_of_baseline);

    % simulate the remaining(t_fin*1000 - duration_of_stimulus - duration_of_baseline) ms after the stimulus arrives
    remainder_baseline_duration = round(t_fin*1000 - (duration_of_baseline + duration_of_stimulus),1);
    [spike_mat_base2, t_vec_base2] = poisson_spike_gen(baseline_freq, remainder_baseline_duration, number_fibers);
    t_vec_base2 = (t_vec_base2 + duration_of_baseline + duration_of_stimulus);
    
    % put the baseline and stimulus periods together
    spike_mat = [spike_mat_base1 spike_mat_stim spike_mat_base2];
    
    % this has dt = 1ms (to change, see PoissonSpikeGen)
    t_vec = [t_vec_base1 t_vec_stim t_vec_base2];

    % convert to seconds:
    t_vec = t_vec./1000.0;

    %% Finish smoothing stuff
    % make an average firing rate from the matrix of 0 and 1 but simply adding
    % the rows (adding the activity of the fibers for each time bin)
    numb_spikes_vec = sum(spike_mat,1); %Number of 1's in a column of spikeMat = Number of spikes along all fibers in 1 time bin
        %Total Number of 1's = total number of spikes across all fibers
    
    % convert this sum to a firing rate:
    total_fr = numb_spikes_vec/(t_vec(2) - t_vec(1))/number_fibers; %number spikes per unit time per fiber--firing rate   
        smooth_width = floor(10/(t_vec(2) - t_vec(1))/1000); %smoothing window in ms (100 ms is the value used in most simulations)
        smoothM = smooth(total_fr, smooth_width); % built-in matlab function
    summing_firing_rates = summing_firing_rates + smoothM; %this is simply smoothM if averaging over multiple time courses of simulations of firing-rates 
        % along all neurons in the bundle of ABeta fibers
end

smooth_m = summing_firing_rates / number_realizations;


end

