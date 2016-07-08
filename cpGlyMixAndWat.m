function cp=cpGlyMixAndWat(x,temp)

% function cp=cpGlyMixAndWat(x,temp)
% This function gives the specific heat [J/kgK] of water and glycol/water
% mixtures given glycol content [%] and fluid temperature [degC] as input.
% Note that the temperature can be also passed as vector. In this case
% the function output will be a vector.

if x>0 % from Conde, M. "Thermophysical properties of brines"
    cp=(4.47642+0.60863*(x/100)-0.71497*(273.15./(temp+273.15))+...
        -1.93855*(x/100)*(273.15./(temp+273.15))+0.47873*(273.15./(temp+273.15)).^2)*1000;
    
elseif x==0  % if glycol%=0, switch to water equation (Furbo's)
    
    cp=4209.1-132.8*10^-2*temp+143.2*10^-4*temp.^2; % [J/kgK]
    
else
   error('Glycol % is not valid!');
end