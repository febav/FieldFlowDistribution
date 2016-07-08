function rho=densityGlyMixAndWat_2(x,temp)

% function rho=densityGlyMixAndWat(x,temp)
% This function computes the density [kg/m3] of water and glycol/water
% mixtures given glycol content [%] and fluid temperature [degC] as input.
% Note that the temperature can be also passed as vector. In this case
% the function output will be a vector.

if x==35   % density for glycol/water mixtures, DTU KEM (35%<x<50%)
    rho=1013-0.2682*temp+0.7225*x-0.00194*temp.^2-0.004964*temp*x; % [kg/m3]
    
elseif x==0     % if glycol%=0, switch to water equation (Furbo's)
    rho=1000.6-0.0128*temp.^1.76; % [kg/m3]
    
elseif x>0  % from Conde, M. "Thermophysical properties of brines"
    rho=508.41109-182.4082*(x/100)+965.76507*(273.15./(temp+273.15))+...
        +280.29104*(x/100)*(273.15./(temp+273.15))-472.2251*(273.15./(temp+273.15)).^2;
    
else
    error('Glycol % must be >=0!');
end