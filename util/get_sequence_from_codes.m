function seq = get_sequence_from_codes(code, rcode, len)
% sequence = get_sequence_from_codes(code, rcode)
% Inputs:
%     code = 1xm or mx1 vector of symbols to repeat (e.g. [1 0 1 0 ... ])
%     rcode = 1xm or mx1 vector of run lengths of each symbol (normalized to unit sum)
%     len = length of sequence (optional if rcode is given un-normalized)
% Outputs:
%     sequence = 1xlen vector of 1's and 0's (Manchester-encoded Barker/Gold/etc. sequence)

if nargin == 2
    len = sum(rcode);
end

if ~all((code == 1) | (code == 0))
    error('Sequence may only consist of 1''s and 0''s!');
end

if ~all(mod(rcode, 1) == 0) % rcode is not given as integer run lengths
   rcode = rcode .* len;
end

seq = []; % initialize sequence
for i = 1:length(code)
    seq = [seq, zeros(1, rcode(i))+code(i)];
end

end