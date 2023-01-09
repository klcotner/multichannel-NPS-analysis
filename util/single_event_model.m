% helper function to generate the model of a single event given shape parameters
% 
% INPUTS
%   * signal_code = vector of 1's and 0's (nominal signal shape)
%   * timestamps = vector of the frame's timestamps [sec]
%   * t_center = timestamp of the center of the event [sec]
%   * t_transit = duration of the event [sec]
%   * amp = event amplitude (R_signal_pore - R_baseline) [data units]
%   * node_amp_factor = ratio (R_signal_node / R_signal_pore)
%   * node_dur_factor = ratio (dT_node / dT_pore)
%   * accel_factor = ratio (V_end / V_start) (how much the cell is speeding up or slowing down during the event)
%   * rise_time = time to rise from 10% to 90% of the difference between R_signal_node and R_signal_pore
% 
% OUTPUT
%   * event_model = vector of frame delta_R for the single event modeled in a single channel

function event_model = single_event_model(signal_code, timestamps, t_center, t_transit, amp, node_amp_factor, node_dur_factor, accel_factor, rise_time)

    % construct normalized edge_times [0, 1]
    seg_times = run_length_encode(signal_code); % duration of each segment
    seg_times = seg_times / sum(seg_times); % normalized to sum to 1
    node_inds = (1+signal_code(1)):2:length(seg_times); % which indices in seg_times correspond to nodes?
    seg_times(node_inds) = seg_times(node_inds) * node_dur_factor; % apply time scaling according to node duration factor
    seg_times = seg_times / sum(seg_times); % normalized to sum to 1
    edge_times = [0, cumsum(seg_times)]; % normalized edge times [0, 1]
    if accel_factor ~= 1
        accel_warp = @(t) (sqrt((accel_factor-1)*(1+accel_factor)*t + 1) - 1) / (accel_factor-1);
        edge_times = accel_warp(edge_times);
    end
    % convert to real edge_times
    edge_times = edge_times * t_transit + (t_center - t_transit/2);
    
    % solve for logistic piece parameters
    dur_edge = rise_time * 3; % define duration of edge as 3X rise_time
    den = 5 * sqrt(7/11) / 4; % modified logistic denominator
    off = (44 - 5*sqrt(77)) / 88; % modified logistic vertical offset
    a = (2*log((9 + sqrt(77))/2)) / rise_time; % logistic parameter
    
    edge_func = @(left, right, t_mid, t) ((1./(1+exp(-a.*(t - t_mid))) - off) ./ den) * (right-left) + left; % logistic function that rises from left to right

    % construct signal (event_model)
    event_model = zeros(size(timestamps));
    
    prev_amp = 0;
    for i = 1:length(edge_times)
        if i==length(edge_times) % this is the last edge
            curr_amp = 0;
        elseif any(node_inds==i) % is node
            curr_amp = node_amp_factor * amp;
        else % is pore
            curr_amp = amp;
        end
        et = edge_times(i);
        edge_mask = timestamps >= (et - dur_edge/2) & timestamps <= (et + dur_edge/2);
        event_model(edge_mask) = edge_func(prev_amp, curr_amp, et, timestamps(edge_mask));
        rest_mask = timestamps > (et + dur_edge/2);
        event_model(rest_mask) = curr_amp;
        prev_amp = curr_amp;
    end
end

function rcode = run_length_encode(code)
    rcode = 1;
    for i = 2:length(code)
        if code(i) == code(i-1)
            rcode(end) = rcode(end) + 1;
        else
            rcode(end+1) = 1;
        end
    end
end