function [sequence, code, rcode, len, num_segments] = get_barker_codes(geom_str)
% [seq, code, rcode, len, num_segments] = get_barker_codes(geom_str)
% Inputs:
%     geom_str = string describing geometry of channel (choice of MB7, MB11, MB13, MB11i, and MB13i)
% Outputs:
%     sequence = 1xlen vector of 1's and 0's (Manchester-encoded Barker sequence)
%     code = 1xnum_segments vector of symbols to repeat (e.g. [1 0 1 0 ... ])
%     rcode = 1xnum_segments vector of run lengths of each symbol (normalized to unit sum)
%     len = length of sequence
%     num_segments = number of discrete segments (also length of code and rcode)

sequence = get_barker_sequence(geom_str);

[code, rcode, len, num_segments] = get_codes_from_sequence(sequence);

end