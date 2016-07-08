function [Ycollector,Interpolated_Dp] = HT3508_Dp_Conde(Vx,Ty,density,InterpMethod)
% HT3508_Dp function calculates the presssure drop across a HT3508 collector
%     [Y,Interpolated_Dp] = HT3508_Dp(Vx,Ty,density,InterpMethod)
%  where 
%  - Vx, volume flow rate in [m3/s]
%  - Ty, fluid temperature in [degC]
%  - density, fluid density in [kg/m3]
%  - Ycollector is the resist.coeff [1/kg.m]
%  - Interpolated_Dp is the collector pressure drop in [Pa]
%  - InterpMethod can be 'linear', 'cubic' or 'spline'. See help interp2 for 
%    more info on the interpolation method
%
% The returned values are based on the data in InterpolationInfo.mat, which
% hence needs to be placed in the same folder as the function file to work
% correctly. InterpolationInfo.mat contains:
%  - Vtot_m3h_vec = x-axis of the tabled pressured drop data [m3/h]
%  - T_vec        = y-axis of the tabled pressured drop data [degC]
%  - DpMatrix     = tabled pressure drop data [kPa]

load InterpolationInfo_Conde.mat

if min(Ty)<-15 | max(Ty)>100 % check if temperature is inside the range
    warning('Temperature (T<-15degC or >100 degC');
end

Interpolated_Dp=interp2(Vtot_m3h_vec,T_vec,DpMatrix,Vx*3600,Ty,InterpMethod)*1000; % pressure drop [Pa], this is why *1000 is added
if sum(sum(isnan(Interpolated_Dp)))>0
    warning('HT3508_Dp function: some of the interpolated values are NaN. This is probably due to values outside range.')
end
Ycollector=Interpolated_Dp./(density.*Vx).^2;
