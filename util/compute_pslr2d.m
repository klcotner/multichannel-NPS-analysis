function pslr2d = compute_pslr2d(r, c, corr_map, index_dict_lr, signal_IDs, code_len, robust, positive)
    % computes 2D peak-to-sidelobe ratio at a given point (r,c) of a correlation map (corr_map)
    % when robust is true (default), sidelobe is computed using median instead of max
    % when positive is true (default), negative sidelobes are excluded
    if nargin < 8
        positive = true;
    end
    if nargin < 7
        robust = true;
    end
    mask = find_sidelobe_mask(r, c, size(corr_map), index_dict_lr, signal_IDs, code_len);
    vals = corr_map(mask);
    if positive % exclude negative sidelobes
        vals = vals(vals > 0);
    end
    if robust
        sidelobe = median(abs(vals));
    else
        sidelobe = max(abs(vals));
    end
    peak = corr_map(r, c);
    pslr2d = peak / sidelobe;
end