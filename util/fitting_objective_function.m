function dif = fitting_objective_function(params, frame, timestamps, baseline, signal_ID_per_chan, signal_code_table, event_signal_IDs, rise_time)
    % this function is meant to parse params into separate input arguments to be passed to forward_model
    % then it returns the difference between y_hat produced by forward_model and the measure y a.k.a. Frame
    % params is a single vector because that's what lsqnonlin wants
    % it is a column vector structured as follows:
    % [event_t_centers; event_t_transits; event_amps; event_node_amp_factors; event_node_dur_factors; event_accel_factors]
    
    if mod(length(params), 6) ~= 0
        error('Length of params must be divisible by 6 because there are 6 params per event in fancy-fit model!');
    elseif length(params) ~= length(event_signal_IDs) * 6
        error('Length of params inconsistent with number of events according to event_signal_IDs');
    end
    
    event_t_centers = params(1:6:end);
    event_t_transits = params(2:6:end);
    event_amps = params(3:6:end);
    event_node_amp_factors = params(4:6:end);
    event_node_dur_factors = params(5:6:end);
    event_accel_factors = params(6:6:end);
    
    y_hat = forward_model(timestamps, baseline, signal_ID_per_chan, signal_code_table, ...
        event_signal_IDs, event_t_centers, event_t_transits, event_amps, ...
        event_node_amp_factors, event_node_dur_factors, event_accel_factors, rise_time);
    
    dif = frame - y_hat;
%     y_hat_d = detrend(y_hat, 'linear');
%     dif = frame - (lowpass(y_hat_d, 200, 10e3) + (y_hat - y_hat_d));

end