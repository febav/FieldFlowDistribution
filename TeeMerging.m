function [Y_Merg,Z_Merg] = TeeMerging(Q_s,Q_c,F_s,F_c,rho)
% TeeDivSide Computes the local loss coefficient for diverter in side
% passage
%   Z_CSi = TeeMerging(Q_s,Q_c,F_s,F_c,rho) = DeltaP/(rho*w_s^2*0.5)
%                   ___________________
%
%                   2s ->         <- 1s
%                   ______       ______
%                         |  |  |
%                         |  V  |
%                         |  c  |
%   Subscript used                   Variables
%   s1 = side passage 1              w  = fluid velocity
%   s2 = side passage 2              Q  = flow rate
%   c  = combined passage            F  = cross area
%                                    M = mass flow rate [kg/s]
%
%   This function computes the local loss coefficient Z_Merg for a converter in
%   side passage according to Idelchick "Handbook of Hydraulic Resistance",
%   CRC press, 3rd edition, 1994 (page 471). It also transforms the same 
%   coefficient into Y_Merg, which expresses the pressure drop as function
%   of the combined mass flow rate M_c, so that
%   Dp= Z_CSi*0.5*rho*w_s^2 = Y_MergTemp*M_c^2 = Y_Merg*M_s^2

% modify the main .m file so that Qs is the flow in each row with the row
% vector index = row number (Matlab numbering). 
Fratio=F_s./F_c;
Qratio=Q_s./Q_c;

A=(Fratio<=0.35)*1+...
(Fratio>0.35).*(Qratio<=0.4)*0.9.*(1-Qratio)+...
(Fratio>0.35).*(Qratio>0.4)*0.55;
Z_Merg=A.*(1+Fratio.^-2+3*(Qratio.^2-Qratio).*Fratio.^-2);    % [-] res.coef. by Idelchick
Y_Merg_temp=Z_Merg*0.5./(F_c.^2.*rho); % [1/m.kg] Dp=Y_MergTemp*M_c^2
Y_Merg=Y_Merg_temp.*Qratio.^-2;        % [1/m.kg] Dp=Y_Merg*M_s^2
end

