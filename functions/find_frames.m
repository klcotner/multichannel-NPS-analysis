% computes automatic frame start and end points based on envelope heuristic

% TODO:
%   * ensure that noise-only chunks also adhere to min frame length
%   * allow multiple noise-only chunks between events (so they can be shorter)

function [bp_inds, Nstarts, Nstops] = find_frames(noise_mask, fs, bp_params)
% Inputs:
%   * noise mask = computed by preprocess_data
%   * fs = sampling rate of noise_mask [Hz]
%   * bp_params = structure containing options
%       * min_framelength_sec = minimum frame length [sec]
%       * min_framemargin_sec = minimum frame margin [sec]
% Outputs:
%   * bp_inds = indices of frame boundaries
%   * Nstarts = starting indices of frames
%   * Nstops = stopping indices of frames

min_framelength_samples = ceil(bp_params.min_framelength_sec * fs);
min_framemargin_samples = ceil(bp_params.min_framemargin_sec * fs);

noise_chunk_dif = diff([noise_mask(1); noise_mask]);
noise_chunk_starts = find(noise_chunk_dif == 1); % extract noise chunk start indices from noise_mask
noise_chunk_stops = find(noise_chunk_dif == -1); % extract noise chunk stop indices from noise_mask
if noise_chunk_stops(1) < noise_chunk_starts(1) % first stop index is before first start index..
    fprintf('find_frames: removing noise chunk stop indices that are before the first start index...\n');
    noise_chunk_stops(noise_chunk_stops < noise_chunk_starts(1)) = [];
end
if noise_chunk_starts(end) > noise_chunk_stops(end) % last start index is after last stop index..
    fprintf('find_frames: removing noise chunk start indices that are after the last stop index...\n');
    noise_chunk_starts(noise_chunk_starts > noise_chunk_stops(end)) = [];
end
if length(noise_chunk_starts) ~= length(noise_chunk_stops)
    error('noise chunk start and stop indices not matched up!');
end

noise_chunk_lens = noise_chunk_stops - noise_chunk_starts + 1; % length of noise chunks [samp]
margin_mask = noise_chunk_lens >= 2 * min_framemargin_samples;
noise_chunk_starts = noise_chunk_starts(margin_mask); % remove chunks that are not long enough according to the min margin
noise_chunk_stops = noise_chunk_stops(margin_mask);
noise_chunk_lens = noise_chunk_lens(margin_mask);
noise_chunk_mids = floor((noise_chunk_starts + noise_chunk_stops) / 2); % extract noise chunk midpoint indices. these are the candidate breakpoints

bp_inds = 1; % temporarily add first index to breakpoint list
for i = 1:length(noise_chunk_mids)
    if noise_chunk_lens(i) >= min_framelength_samples + 2*min_framemargin_samples % generate an all-noise frame
        if noise_chunk_starts(i) + min_framemargin_samples - bp_inds(end) + 1 >= min_framelength_samples
            bp_inds = [bp_inds; noise_chunk_starts(i) + min_framemargin_samples; noise_chunk_stops(i) - min_framemargin_samples];
        else
            bp_inds = [bp_inds; noise_chunk_stops(i) - min_framemargin_samples];
        end
    elseif noise_chunk_mids(i) - bp_inds(end) >= min_framelength_samples
        bp_inds = [bp_inds; noise_chunk_mids(i)]; % add to breakpoint list
    end
end
bp_inds(1) = []; % remove first index from breakpoint list

Nstarts = [1; bp_inds];
Nstops = [bp_inds-1; length(noise_mask)];
if Nstops(end)-Nstarts(end)+1 < min_framelength_samples % if the last frame is too short, combine it with the penultimate frame
    Nstarts = Nstarts(1:end-1);
    Nstops = Nstops(1:end-1);
end

end