function sequence = get_barker_sequence(geom_str)
% seq = get_barker_sequence(geom_str)
%   Input: geom_str = string describing geometry of channel (choice of MB7, MB11, MB13, MB11i, and MB13i)
%   Output: sequence = 1xlen vector of 1's and 0's (Manchester-encoded Barker sequence)

switch geom_str
    case 'MB7'
        sequence = [1 0 1 0 1 0 0 1 0 1 1 0 0 1];
    case 'MB7i'
        sequence = [0 1 0 1 0 1 1 0 1 0 0 1 1 0];
        
    case 'MB11'
        sequence = [1 0 1 0 1 0 0 1 0 1 0 1 1 0 0 1 0 1 1 0 0 1];
    case 'MB11i'
        sequence = [0 1 0 1 0 1 1 0 1 0 1 0 0 1 1 0 1 0 0 1 1 0];
        
    case 'MB13'
        sequence = [1 0 1 0 1 0 1 0 1 0 0 1 0 1 1 0 1 0 0 1 1 0 0 1 1 0];
    case 'MB13i'
        sequence = [0 1 0 1 0 1 0 1 0 1 1 0 1 0 0 1 0 1 1 0 0 1 1 0 0 1];
        
    otherwise
        error('Barker code string not recognized!');
end

end