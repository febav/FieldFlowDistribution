function [Y_STAD,DpValve]=STAD_valve(Flow,rho,STADset,Nstad20)

% Nstad20 is the number of STAD 20 valves
V=Flow/rho*1000*3600; % [l/h]; convert Flow [kg/s] into [l/h]

% fitting of STAD20, for turns=[2.6-2.9] and V=[320-2500 l/h]
Kv20_1fun=@(turns,V_lh) -2.6837+2.17489*turns-1.78e-05*V_lh+...
                        +1.7e-05*turns.*V_lh-8.40274e-09*V_lh.^2;

% fitting of STAD20, for turns=[2-4] and V=[320-2500 l/h]
p00=21.9476;
p10=-31.43406;
p01=8.75506e-06;
p20=17.02211;
p30=-3.77447;
p40=0.307239;
Kv20_2fun=@(turns,V_lh) p00+p10*turns+p01*V_lh+p20*turns.^2+...
                        p30*turns.^3+p40*turns.^4;  
 
% fitting of STAD25, for turns=[2-4] and V=[320-2500 l/h]
Kv25fun=@(turns,V_lh) -1.712516+0.8661*turns-1.6384e-05*V_lh+...
                  +1.34175*turns.^2+1.434e-05*turns.*V_lh-0.22905*turns.^3;
if Nstad20>0
    cond=(STADset(1:Nstad20)<=2.9).*(STADset(1:Nstad20)>=2.6); % decide when to use different fitting curves for STAD20
    Kv20_1val=Kv20_1fun(STADset(1:Nstad20),V(1:Nstad20)).*cond;
    Kv20_2val=Kv20_2fun(STADset(1:Nstad20),V(1:Nstad20)).*(~cond);
    Kv20val=Kv20_1val+Kv20_2val;
    if Nstad20==length(STADset)
        DpValve=0.1*V.^2./(Kv20val.^2); % [Pa]
    else 
        Kv25val=Kv25fun(STADset(Nstad20+1:end),V(Nstad20+1:end));
        DpValve=0.1*V.^2./([Kv20val;Kv25val].^2); % [Pa]
    end  
else
    Kv25val=Kv25fun(STADset(Nstad20+1:end),V(Nstad20+1:end));
    DpValve=0.1*V.^2./(Kv25val.^2); % [Pa]
end
Y_STAD=DpValve./Flow.^2; % Y_STAD has same unit and ref as YR -> u can sum them up
end

