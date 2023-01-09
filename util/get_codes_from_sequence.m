function [code, rcode, len, num_segments] = get_codes_from_sequence(seq)
% [code, rcode, len, num_segments] = get_codes_from_sequence(sequence)
% Inputs:
%     sequence = 1xlen or lenx1 vector of 1's and 0's (Manchester-encoded Barker/Gold/etc. sequence)
% Outputs:
%     code = 1xnum_segments vector of symbols to repeat (e.g. [1 0 1 0 ... ])
%     rcode = 1xnum_segments vector of run lengths of each symbol (normalized to unit sum)
%     len = length of sequence
%     num_segments = number of discrete segments (also length of code and rcode)

if ~all((seq == 1) | (seq == 0))
    error('Sequence may only consist of 1''s and 0''s!');
end

len = length(seq);
s_ind = 1; % sequence index
c_ind = 1; % code/rcode index
code = []; % initialize code
rcode = []; % initialize rcode
r_num = 1; % initialize run length counter

code(c_ind) = seq(s_ind); % insert first symbol
s_ind = s_ind + 1; % go to next symbol in sequence
while s_ind <= len
    if seq(s_ind) == seq(s_ind-1) % current symbol same as last symbol
        r_num = r_num + 1; % increment run length
    else
        rcode(c_ind) = r_num; % insert run length of last symbol
        c_ind = c_ind + 1; % increment code/rcode index
        code(c_ind) = seq(s_ind); % insert current symbol
        r_num = 1; % reset run length counter
    end
    s_ind = s_ind + 1; % go to next symbol 
end
rcode(c_ind) = r_num; % insert run length of last symbol
rcode = rcode ./ len; % normalize rcode for some reason

num_segments = length(rcode);

end