% function to calculate the list of diameters and the measures of central tendancy
% > NB: num_ch can be a vector of the same length as eventlist, if needed
% * depends on d_sphere function from coded-nps-software/Util

function [diameter_list, d_mean, d_median, d_peak] = calc_diameters(D_eff, event_list, num_ch, L_ch)

    % calculate diameter list
    scamp_list = [event_list.amp_singlechannel]';
    blamp_list = [event_list.baseline_amp_center]';
    diameter_list = d_sphere(D_eff, scamp_list, blamp_list.*num_ch, L_ch);

    if nargout>1
        % calculate mean diameter
        d_mean = mean(diameter_list);
    
        if nargout>2
            % calculate median diameter
            d_median = median(diameter_list);
        
            if nargout>3
                % calculate peak of diameter kernal density estimate
                [ksfreq, ksx] = ksdensity(diameter_list,'numpoints',1000);
                [~,xi] = max(ksfreq);
                d_peak = ksx(xi); % um
            end
        end
    end
    
end