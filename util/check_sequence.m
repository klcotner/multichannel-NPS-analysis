% helper function to verify that the input is a valid sequence
% a valid sequence must be a 1xn vector of 1's and 0's
% note: '../Util' must be on the search path (for check_sequence)
function tf = check_sequence(x)
    % x must be a numeric vector
    if isnumeric(x) && isvector(x)
        % x must be a row vector with >1 element
        if size(x,1)==1 && size(x,2)>1
            % all elements of x must be either 1 or 0
            tf = sum(arrayfun(@(xx) ~check_bool(xx), x))==0;
        else
            tf = false;
        end
    else
        tf = false;
    end
end