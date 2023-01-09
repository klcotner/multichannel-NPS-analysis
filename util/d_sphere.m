% Deblois & Bean, as derived by KLC 2022
% 
% `D_e`(double) = effective diameter [length unit]
% `dR` (double) = change in resistance caused by particle transit event [resistance unit]
% `R` (double) = baseline resistance [resistance unit]
% `L` (double) = lenth of channel [length unit]

function diameter = d_sphere(D_e, dR, R, L) ...
    
    diameter = ...
    (   ( (dR ./ R) .* D_e.^2 .* L ) ./ ...
        ( 1 + 0.8 .* (dR ./ R) .* (L ./ D_e) ) ...
    ) .^ (1/3);

end