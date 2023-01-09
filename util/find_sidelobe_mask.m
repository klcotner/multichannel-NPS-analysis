function [sidelobe_mask, mainlobe_mask] = find_sidelobe_mask(r, c, sz_corr_map, index_dict_lr, signal_IDs, code_len)
    % helper function to find the indices within correlation map that
    % contain energy originating from the detected signal template at (r,c)

    sig_len = index_dict_lr(r,1); % length of detected template [samp]
%     halfwidth_per_row = floor(sig_len/2 + index_dict_lr(:,1)/2) - 1; % length in samples to look to the left and right in each row
    quarterwidth_per_row = floor(sig_len/4 + index_dict_lr(:,1)/4) - 1; % length in samples to look to the left and right in each row
    pulsewidth_per_row = ceil((sig_len/2 + index_dict_lr(:,1)/2) / code_len) - 1; % length of single pulse {0 or 1} in samples
    sidelobe_mask = false(sz_corr_map); % initialize mask
    mainlobe_mask = false(sz_corr_map);
    signal_ID = signal_IDs(r); % get signal ID of detected template
    r_inds = find(signal_IDs == signal_ID); % rows indices of correlation map to consider
    for i = 1:length(r_inds) % iterate over rows
        r_ind = r_inds(i);
        c_inds = (c-quarterwidth_per_row(r_ind)):(c+quarterwidth_per_row(r_ind)); % column indices to consider looking for sidelobes
        c_inds(c_inds<1) = 1;
        c_inds(c_inds>sz_corr_map(2)) = sz_corr_map(2);
        sidelobe_mask(r_ind, c_inds) = true;
        c_inds = (c-pulsewidth_per_row(r_ind)):(c+pulsewidth_per_row(r_ind)); % column indices of main lobe to exclude
        c_inds(c_inds<1) = 1;
        c_inds(c_inds>sz_corr_map(2)) = sz_corr_map(2);
        sidelobe_mask(r_ind, c_inds) = false;
        mainlobe_mask(r_ind, c_inds) = true;
    end
end