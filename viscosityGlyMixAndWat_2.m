function nu=viscosityGlyMixAndWat_2(x,temp)

% function nu=viscosityGlyMixAndWat(x,temp)
% This function gives the kinem. viscosity [m2/s] of water and glycol/water
% mixtures given glycol content [%] and fluid temperature [degC] as input.
% Note that the temperature can be also passed as vector. In this case
% the function output will be a vector.

if x==35 % sample from Hoeje Taastrup
    nu=((temp<=38).*(-1.449*10^-5*temp.^3+3.066*10^-3*temp.^2-0.2337*temp+7.289)...
        +((temp>38).*(180.294*temp.^-1.232)))./1000./...
        densityGlyMixAndWat_2(x,temp); % [m2/s]
    
elseif x==0  % if glycol%=0, switch to water equation (Kestin's)
    nu=1.002/1000*10.^((20-temp)./(temp+96).*(1.2378-1.303*10^-3*(20-temp)+...
        3.06*10^-6*(20-temp).^2+2.55*10^-8*(20-temp).^3))./...
        densityGlyMixAndWat_2(x,temp); % [m2/s]
    
elseif x>0  % from Conde, M. "Thermophysical properties of brines"
    nu=exp(-1.02798-10.03298*(x/100)-19.93497*273.15./(temp+273.15)+...
        14.65802*(x/100)*273.15./(temp+273.15)+14.6205*(273.15./(temp+273.15)).^2)./...
        densityGlyMixAndWat_2(x,temp); % [m2/s]
else
    error('Glycol % is not valid!');
end