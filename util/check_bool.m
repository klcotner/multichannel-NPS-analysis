% helper function to identify a numeric or logical "boolean"
function tf = check_bool(x)
    % x can be a 1x1 logical
    if islogical(x) && isscalar(x)
        tf = true;
    % x can be a scalar with the value 1 or 0
    elseif isnumeric(x) && isscalar(x)
        tf = x==1 || x==0;
    % otherwise, x is not a boolean
    else
        tf = false;
    end
end