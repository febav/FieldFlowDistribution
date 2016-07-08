clear all
close all
clc

%% LAYOUT OF SOLAR COLLECTOR FIELD
% collector array
Np=24;     % number of hydraulic paths (collector rows)
Nsub1=12;  % number of rows in subfield
NCol=10;   % number of collectors per row

% supply and return pipes
Supply=load('SupplyPipesHT_constantD.txt'); % load supply pipes info
Return=load('ReturnPipesHT_constantD.txt'); % load return pipes info
DCi=Supply(:,3:4); % [m] diameter supply pipes
DCo=Return(:,3:4); % [m] diameter return pipes
LCi=Supply(:,1:2); % [m] length of supply segments
LCo=Return(:,1:2); % [m] length of return segments
ACi=DCi.^2*pi/4;
ACo=DCo.^2*pi/4;
ffrLaw=3; % friction factor: 1=Blasius 2=Colebrook 3=Haaland 4=Joseph&Jang
epLog=0.1e-3; % [m] surface roughness in distribution pipes

% hose and STAD valves
hoseDp=@(V_hose) 2.537e9*V_hose.^2-3.107e5*V_hose; % [Pa] if V_hose is in [m3/s]
setSTAD20=[2.40;2.43;2.45;2.47;2.50;2.52;2.55;2.57;2.60;2.68;2.73];
setSTAD25=[1.98;2.00;2.04;2.08;2.16;2.22;2.30;2.35;2.60;2.80;3.00;3.10;3.20];
STADset=[setSTAD20;setSTAD25];
valvesOnOff=0;      % 0 = no valves; 1 = valves are is use

% opstander
Dupop=0.032;        % [m] diameter of the pipe just before/after 1st/last hose
Aupop=Dupop^2*pi/4; 
Lupop=0.809;        % [m] length of //
Leq_OpIn=0.97;      % [m] equivalent length of the inlet opstander (assuming Din=32mm)
Leq_OpOut=1.73;     % [m] equivalent length of the outlet opstander (assuming Din=32mm)
Leq_OpInTot=Leq_OpIn+Lupop;   % [m] equivalent length
Leq_OpOutTot=Leq_OpOut+Lupop; % [m] equivalent length
Leq_OpOutMerg=2.98;           % [m] equivalent length vert pipe downstream Y-tee in outlet opstander
DSideDiv=0.0372;              % [m] diameter of the side pipe parting from the main supply pipe Do=42.4mm (series1)
ASideDiv=DSideDiv^2*pi/4;     % [m2] area of the pipe parting from the supply pipe
DSideConv=0.0545;             % [m] diameter of the pipe flowing into the main return pipe Do=60.3mm (series2)
ASideConv=DSideConv^2*pi/4;   % [m2] area of the pipe flowing into the return pipe
AMerg_s=8.553e-4*ones(Np,1);  % [m2] area of the side pipe merging in the symmetrical Y (Di=33mm)
AMerg_c=AMerg_s;              % [m2] area of the combined pipe resulting from the symmetrical Y
epGeb=0.01e-3;                % [m] surface roughness 2

% common supply/return pipes before branching
DCommIn=0.1071;         % [m] diameter of common supply pipe 
ACommIn=DCommIn^2*pi/4; % [m] area of common supply pipe
LCommIn=13.3;           % [m] overall length of common supply pipe
ZlocCommIn=1.8+0.18;    % [-] sum of Zloc's along //
LCommOut=13;            % [m] overall length of common return pipe
ZlocCommOut=3+1.43;     % [-] sum of Zloc's along //

%% FLUID PROPERTIES
Vtot_m3h_vector=[8 15 20 30 40 50]; % [m3/h] field flow rates to be tested
RMSD_vector=ones(1,length(Vtot_m3h_vector));
ReRow=RMSD_vector;
Tin=55;             % [degC] fluid (inlet) temperature
Tout=95;            % [degC] fluid (outlet) temp. (only 1st guess, if Tmode=1)
Tmode=1;            % =1 if temp=f(flow), !=1 Tout is fixed
Tm=(Tin+Tout)/2;    % [degC] mean fluid temp.
Tio_alongRow=linspace(Tin,Tout,NCol+1); % [degC] fluid temperature along at collector inlet/outlet
Tm_alongRow=[Tio_alongRow(1),Tio_alongRow(1:end-1)+(Tout-Tin)/NCol/2,Tio_alongRow(end)]; % [degC] fluid mean temperature profile along row + field in/outlet
TimoMap=repmat(Tm_alongRow,Np,1);
x=35.01;                                % [%] glycol content (=0 for water)
rhoMap=densityGlyMixAndWat_2(x,TimoMap);
rhoIn=rhoMap(1,1);                      % [kg/m3] fluid density at inlet
rhoOut=rhoMap(end,end);                 % [kg/m3] fluid density at outlet
nuIn=viscosityGlyMixAndWat_2(x,Tin);    % [m2/s] fluid kin.viscosity at inlet
nuOut=viscosityGlyMixAndWat_2(x,Tout);  % [m2/s] fluid kin.viscosity at outlet
nuROut=nuOut*ones(Np,1);                % [m2/s] kin.visc. at outlet of each row
nuYmerg=nuOut*ones(Nsub1,1);            % [m2/s] kin.visc. at Y-merging
nuCo=nuOut*ones(Nsub1,1);               % [m2/s] kin.visc. in return pipe segments
DpField_NoDh=zeros(1,length(Vtot_m3h_vector));  % [Pa] Dp_field (no hydrostatic)
MeasuredDp=DpField_NoDh;    % [Pa] Dp between manometers , including Dh
%% WEATHER and COLLECTOR EFFICIENCY
if Tmode==1
    ACol=12.6;  % [m2] HT3508 aperture area
    eta0=0.824; % [-] peak efficiency, based on aperture area, //
    a1=1.865;   % [W/m2K]  1st loss coef., //
    a2=0.018;   % [W/m2K2] 2nd loss coef., //
    Tamb=20;    % [degC] ambient temperature
    theta=0;    % [deg] incidence angle
    IAM=1-(tand(theta/2))^3.8; % iam
    cp=cpGlyMixAndWat(x,mean([Tin,Tout]));  % [J/kgK] specific heat
    Afield=ACol*NCol*Np; % [m2] field aperture area
    Gtot=(Vtot_m3h_vector/3600*rhoIn*cp*(Tout-Tin)+... % [W/m2] solar irradiance
        Afield*a1*(Tm-Tamb)+Afield*a2*(Tm-Tamb)^2)/(Afield*eta0*IAM);
end

%% Iteration parameters
itmax=40;       % max number of iteration
itend=itmax;
tolFlow=0.1e-4; % tolerance
Anal=fopen('Analysis.csv','w+');
fid=1;

% SAVE OUTPUT RESULT ON EXTERNAL FILE
outfile=fopen('FlowDistribution.txt','w');
fprintf(outfile,'Tin=%4.1f degC and Tout=%4.1f degC\r\n',Tin,Tout);
fprintf(outfile,'Flow[kg/h]\tDp[kPa]\trow1\trow2\trow3\trow4\trow5\trow6\trow7\trow8\trow9');
fprintf(outfile,'\trow10\trow11\trow12\trow13\trow14\trow15\trow16\trow17\trow18');
fprintf(outfile,'\trow19\trow20\trow21\trow22\trow23\trow24\titer#\r\n');

figure(2)
set(gcf,'DefaultAxesColorOrder',jet(length(Vtot_m3h_vector)));
legendFig2=cell(1,length(Vtot_m3h_vector));
hold on
if length(Vtot_m3h_vector)==1 % field distribution VS iteration#
    figure(1)
    set(gcf,'DefaultAxesColorOrder',jet(itmax/2)); % jet=dark blue->red
end

tic

for kk=1:length(Vtot_m3h_vector)
    
    % first guess on flow distribution in the field
    Vtot=Vtot_m3h_vector(kk)/3600;  % [m3/s] total volume rate in input
    Mtot=Vtot*rhoIn;                % [kg/s] total mass rate in input
    MoldR=Mtot/Np*ones(Np,1); % [kg/s] mass flow in rows (uniformity is assumed)
    MoldCi=Mtot/2*ones(Np,1); % MoldC(1)=MoldC(13)=Mtot/2 (other cells will be overwritten)
    for jj=2:Nsub1 
      MoldCi(jj)=MoldCi(jj-1)-MoldR(jj-1); % mass flow in the connecting pipes 
    end
    MoldCi(Nsub1+1:end)=MoldCi(1:Nsub1); % flow distribution in subfield1=subfield2
    MnewCi=MoldCi;
    
    Dp=1;
    it=1;
    BB = [Mtot zeros(1,Np-1)]';  % right hand side (RHS) column vector

    while it<=itmax     % matrix iteration till convergence of solution
        if length(Vtot_m3h_vector)==1   % if only 1 flow rate is analysed
            figure(1)                   % field distribution VS iteration#
            plot(1:Np,MoldR/Mtot*100)
            xlabel('Row #');ylabel('fraction of total flow \chi [%]')
            grid on; hold all
        end
        
        % update of pressure drop and mass flow distribution
        DpOld=Dp;  % [Pa] pressure distribution from the previous iteration
        MoldCo=MoldCi(1:Nsub1)+MoldCi(Nsub1+1:Np);
        Mmerg=MoldR(1:Nsub1)+MoldR(Nsub1+1:end);
        if Tmode==1
          % temperature profile along rows
          cp=cpGlyMixAndWat(x,mean([Tin,mean(TimoMap(:,end))])); % [J/kgK] specific heat
          Ra=-a2./MoldR/cp;
          Rb=(-a1+2*Tamb*a2)./MoldR/cp;
          Rc=(eta0*IAM*Gtot(kk)+a1*Tamb-a2*Tamb^2)./MoldR/cp;
          Rd=sqrt(4*Ra.*Rc-Rb.^2);
          Rk1=2*atan((Tin*2*Ra+Rb)./Rd)./Rd;
          TimoMap=(repmat(Rd,1,NCol+2).*tan(0.5*(repmat(Rk1.*Rd,1,NCol+2)+...
              Rd*ACol*[0,0.5:NCol-0.5,NCol]))-repmat(Rb,1,NCol+2))./(2*...
              repmat(Ra,1,NCol+2));
          rhoMap=densityGlyMixAndWat_2(x,TimoMap);
          nuROut=viscosityGlyMixAndWat_2(x,TimoMap(:,end));
          nuYmerg=(nuROut(1:Nsub1).*MoldR(1:Nsub1)+nuROut(Nsub1+1:Np).*MoldR(Nsub1+1:Np))...
              ./Mmerg;
          nuCo(Nsub1,1)=nuYmerg(Nsub1);
          for ij=1:(Nsub1-1)
            nuCo(Nsub1-ij)=(nuCo(Nsub1-ij+1)*MoldCo(Nsub1-ij+1)+nuYmerg(Nsub1-ij)...
                *Mmerg(Nsub1-ij))./(MoldCo(Nsub1-ij+1)+Mmerg(Nsub1-ij)); 
          end
        end
        VoldRIn=MoldR/rhoIn;  % [m3/s] volume flow rate in row
        VoldCi=MoldCi/rhoIn;  % [m3/s] vol.flow rate in conn. suppply pipes
        VoldCo=MoldCo/mean(rhoMap(:,end)); % [m3/s] vol.flow rate in conn. return pipes
        
        % Dp along row (collectors+hoses)
        VoldRmap=repmat(MoldR,1,NCol+2)./rhoMap; % vol.flow rate at inlet [m3/s]
        [YColl,DpColl]=HT3508_Dp_Conde(VoldRmap(:,2:end-1),TimoMap(:,2:end-1),rhoMap(:,2:end-1),'cubic');
        YColl_sum=sum(YColl,2); % resistance coeff. for all collectors of each row
        DpHose=hoseDp(VoldRmap(:,5))*(NCol+1); % [Pa] press.drop in hoses
        Yhoses=DpHose./MoldR.^2; % Y coeff for hoses in the same row
        
        % Dp in opstanders
        [Ystad,DpValve]=STAD_valve(MoldR,rhoIn,STADset,length(setSTAD20)); % [1/kg.m] STAD valves
        Ystad=Ystad*valvesOnOff;
        ReOpIn=Dupop./(nuIn*Aupop).*abs(MoldR/rhoIn);  % [-] Re in inlet opstander
        fOpIn=FrictionFactorFunc_matrix(ReOpIn,Dupop*ones(size(ReOpIn)),epGeb,ffrLaw); % [-] fr.f. //
        YOpIn=fOpIn.*Leq_OpInTot./Dupop*0.5/rhoIn./(Aupop.^2); % [1/kg.m] resistance coeff. //
        ReOpOut=Dupop./(nuROut*Aupop).*abs(MoldR./rhoMap(:,end));  % [-] Re in outlet opstander, before Y-merging
        fOpOut=FrictionFactorFunc_matrix(ReOpOut,Dupop*ones(size(ReOpOut)),epGeb,ffrLaw); % [-] fr.f. //
        YOpOut=fOpOut.*Leq_OpOutTot./Dupop*0.5./rhoMap(:,end)./(Aupop.^2); % [1/kg.m] resistance coeff. //
        ReOpOutMerg=DSideConv./(nuYmerg*ASideConv).*abs(Mmerg)./(0.5*(rhoMap(1:Nsub1,end)+rhoMap(Nsub1+1:Np,end)));  % [-] Re in outlet opstander, after Y-merging
        fOpOutMerg=FrictionFactorFunc_matrix(ReOpOutMerg,DSideConv*ones(size(ReOpOutMerg)),epLog,ffrLaw); % [-] fr.f. //
        YOpOutMerg=fOpOutMerg.*Leq_OpOutMerg./DSideConv*0.5./(0.5*(rhoMap(1:Nsub1,end)+rhoMap(Nsub1+1:Np,end)))./(ASideConv.^2); % [1/kg.m] res. coef. //
        
        % Dp in supply/return pipes
        ReCi=DCi./(nuIn*ACi).*abs(repmat(MoldCi,1,2)/rhoIn); % [-] Re in supply pipes
        ReCo=DCo./(repmat(nuCo,1,2).*ACo).*abs(repmat(MoldCo,1,2)/mean(rhoMap(:,end)));% [-] Re in return pipes
        fCi=FrictionFactorFunc_matrix(ReCi,DCi,epLog,ffrLaw); % [-] fr.f. in supply pipes
        fCo=FrictionFactorFunc_matrix(ReCo,DCo,epLog,ffrLaw); % [-] fr.f. in return pipes
        YCi=sum(fCi.*LCi./DCi*0.5/rhoIn./(ACi.^2),2); % [1/kg.m] res.coef. supply pipes.
        YCo=sum(fCo.*LCo./DCo*0.5/mean(rhoMap(:,end))./(ACo.^2),2); % [1/kg.m] res.coef. return pipes. YCi/o(2) refers to the segment between row 1 and 2, and so on
        YTot=YColl_sum+Yhoses+YOpIn+YOpOut+Ystad; % YT inizialized as Y(row components), as it includes at least the row Dp
        Yrow=YTot;
        [YTinSi,ZTinSi]=TeeDivSide(VoldRIn,VoldCi,ASideDiv,ACi(:,2),ACi(:,2),rhoIn,ReCi(:,2)); % diverter=in (side)
        [YTinSt,ZTinSt]=TeeDivSt(VoldRIn,VoldCi,ASideDiv,ACi(:,2),ACi(:,2),rhoIn,ReCi(:,2));   % diverter=in (straight)
        [YToutSi,ZToutSi,ZToutSiL]=TeeConvSide(Mmerg/mean(rhoMap(:,end)),VoldCo,ASideConv,ACo(:,2),ACo(:,2),mean(rhoMap(:,end)),ReCo(:,2)); % converter=out (side)
        [YToutSt,ZToutSt]=TeeConvSt(Mmerg/mean(rhoMap(:,end)),VoldCo,ASideConv,ACo(:,2),ACo(:,2),mean(rhoMap(:,end)),ReCo(:,2),ZToutSiL);   % converter=out (straight)  
        [YTMerg,ZTMerg]=TeeMerging(MoldR/mean(rhoMap(:,end)),repmat(Mmerg/mean(rhoMap(:,end)),2,1),AMerg_s,AMerg_c,mean(rhoMap(:,end)));   % merging Y (outlet)
        DpSupplyPipes=YCi.*MoldCi.^2;   % Dp in supply pipes segments
        DpReturnPipes=YCo.*MoldCo.^2;   % Dp in return pipes segments
        DpTeeInSt=YTinSt.*MoldCi.^2;    % Dp in straight passage of tee-inlet
        DpTeeOutSt=YToutSt.*MoldCo.^2;  % Dp in straight passage of tee-outlet
        YTot=YTot+YTMerg+...               % adding merging Y at row outlet
          (repmat(YToutSi.*MoldCo.^2,2,1)+...% adding outlet T-side
          YTinSi.*MoldCi.^2.+...             % adding inlet T-side
          repmat(YOpOutMerg.*Mmerg.^2,2,1))...% Logstor pipe between merging Y and return pipes
          ./MoldR.^2;                         % normalize the Dp to MR 

      for jj=1:Np   % for each jj-th row (see Fig. 1a)  
        if jj<=Nsub1            
          YTot(jj)=YTot(jj)+...
            (sum(DpSupplyPipes(1:jj))+...     % supply pipes          
            sum(DpReturnPipes(1:jj))+...      % return pipes
            (jj>1)*sum(DpTeeInSt(1:jj-1))+... % T-in (div) straight
            (jj>1)*sum(DpTeeOutSt(1:jj-1)))...% T-out(con) straight
            ./MoldR(jj)^2;                    % normalize Dp to MR             
        else
          YTot(jj)=YTot(jj)+...
            (sum(DpSupplyPipes(Nsub1+1:jj))+...           % supply pipes
            sum(DpReturnPipes(1:jj-Nsub1))+...            % return pipes
            (jj>Nsub1+1)*sum(DpTeeInSt(Nsub1+1:jj-1))+... % T-in (div) straight
            (jj>Nsub1+1)*sum(DpTeeOutSt(1:jj-Nsub1-1)))...% T-out (con) straight
            ./MoldR(jj)^2;                                % normalize Dp to MR  
        end
      end

      AA=diag(-YTot.*MoldR,0)+diag(YTot(1:end-1).*MoldR(1:end-1),-1);
      AA(1,:)=1;                        % coefficent matrix
      MnewR = AA\BB;                    % find solution vector
      MnewCi(1)=sum(MnewR(1:Nsub1));    % flow to subfield1
      MnewCi(Nsub1+1)=Mtot-MnewCi(1);	% flow to subfield2
      for jj=2:Nsub1
        MnewCi(jj)=MnewCi(jj-1)-MnewR(jj-1);% flow in connecting pipes (subfield 1)
        MnewCi(Nsub1+jj)=MnewCi(Nsub1+jj-1)-MnewR(Nsub1+jj-1); %  // (subfield 2)
      end
      
      % check convergence
      emaxFlow=max(abs(MnewR-MoldR)./MnewR);
      Dp=[MnewR(1:end-1).*diag(AA,-1);MnewR(end)*(-AA(end,end))];

      if length(Vtot_m3h_vector)==1   % if only 1 flow rate is analysed
          fprintf(1,'\nIntermediate flow distribution \n');
          fprintf(1,'  iter# = %2d     max error = %8.2e\n',it,emaxFlow);
          fprintf(1,  'Row  MnewR  MoldR  MnewC  MoldC\n');
          for j = 1:length(MoldR)
            fprintf(1,'%2i  %6.3f  %5.3f  %5.3f  %5.3f\n',j,MnewR(j),MoldR(j),MnewCi(j),MoldCi(j));
          end
      end

      if emaxFlow<tolFlow
        itend=it;
        it=itmax+1; 
      else 
        it=it+1;
        MoldR=(MnewR+MoldR)/2;
        MoldCi=MnewCi;
      end
    end
    if itend==itmax 
        warning('Hit max number of iterations!!! \n');
    end
    
    figure(2)
    plot(1:Np,MoldR/(Mtot/Np),'d-','markers',4)
    hold on
    legendFig2{kk}=['V=',num2str(Vtot_m3h_vector(kk)),' m^3/h'];
    
    % print flow distribution, DpTot and iter to file
    fprintf(outfile,'%9.1f\t%6.2f',Mtot*3600,Dp(1)/1000);
    fprintf(outfile,'\t%4.2f%%',MoldR/Mtot*100);
    fprintf(outfile,'\t%5.0f\r\n',itend);
    
    % compute parameters based on the last flow rates
    VRIn=MnewR/rhoIn;     % [m3/s] volume flow rate in different rows
    vR=VRIn./18/(0.0073^2*pi/4); % [m/s] average flow velocity in absorber
    VCiIn=MnewCi/rhoIn;    % [m3/s] volume flow rate in supply pipes
    MnewCo=MnewCi(1:Nsub1)+MnewCi(Nsub1+1:Np); % flow rate in return pipes
    VCoOut=MnewCo./mean(rhoMap(:,end)); % [m3/s] volume flow rate in return pipes
    vCi=repmat(VCiIn,1,2)./ACi; % [m/s] flow velocity conn. pipes (inlet)
    vCo=repmat(VCoOut,1,2)./ACo; % [m/s] flow velocity conn. pipes (outlet)
    DpR=(YColl_sum+Yhoses).*MnewR.^2;    % [Pa] press. drop along rows (coll+hose)
    DpTinSi=YTinSi.*MnewCi.^2; % [Pa] press.drop in diverging tees (side)
    DpTinSt=YTinSt.*MnewCi.^2; % [Pa] press.drop in diverging tees (straight)
    DpCi=YCi.*MnewCi.^2;  % [Pa] Pressure drop in inlet conn. pipe segment
    DpCo=YCo.*MnewCo.^2;  % [Pa] Pressure drop in outlet conn. pipe segment

    % MeasuredDp calculations  
    vCommIn=Mtot/rhoIn/ACommIn; % [m/s] velocity along the common forward pipe
    ReCommIn=DCommIn*vCommIn/nuIn; % [-] Re along //
    fCommIn=FrictionFactorFunc_matrix(ReCommIn,DCommIn,epLog,ffrLaw); % [-] fr.f. along //
    DpCommInFr=fCommIn*LCommIn/DCommIn*0.5*rhoIn*vCommIn^2; % [Pa] fricton Dp along // 
    DpCommInLoc=ZlocCommIn*0.5*rhoIn*vCommIn^2; % [Pa] local loss Dp //
    DpCommInTot=DpCommInFr+DpCommInLoc; % [Pa] total Dp //
    DpCommOutFr=fCo(1,1)*LCommOut/DCo(1,1)*0.5*mean(rhoMap(:,end))*vCo(1,1)^2; % [Pa] fricton Dp along common return pipe (from P02 to TB_h0)
    DpCommOutLoc=ZlocCommOut*0.5*mean(rhoMap(:,end))*vCo(1,1)^2; % [Pa] local loss Dp //
    DpCommOutTot=DpCommOutFr+DpCommOutLoc; % [Pa] total Dp //
    DhMan=1.365;    % [m] height difference manometers (H_return-H_supply)
    DpStatic=9.81*DhMan*densityGlyMixAndWat_2(35,mean([Tin,mean(TimoMap(:,end))]));  % [Pa] hydrostatic Dp due to manometers' Dh
    DpField_NoDh(kk)=Dp(1)+DpCommInTot+DpCommOutTot;  % [Pa] Overall press. drop
    MeasuredDp(kk)=DpField_NoDh(kk)+DpStatic; % [Pa] Dp between 2 manometer
    RMSD_vector(kk)=sqrt(1/Np*sum((MoldR/(Mtot/Np)-1).^2));
      
    % print summary results
    fprintf(fid,'\n\nResults for field flow rate [kg/h], %10.1f \n',Mtot*3600);
    if Tmode==1
      fprintf(fid,'Tmode=%d : row temp.profiles are f(V'',G_sol,T_amb)\n',Tmode);
    else
      fprintf(fid,'Tmode=%d : user defined outlet temperature\n',Tmode);
    end
    fprintf(fid,'Calculated DP across parallel [kPa]:  %7.2f\n',Dp(1)/1000);
    fprintf(fid,'RMSD of dimensionless flow distribution:  %6.6f\n',RMSD_vector(kk));
    fprintf(fid,'Total DP from PT7023 to PT7028 (with Dh) [kPa]: %7.2f\n',MeasuredDp(kk)/1000);
    fprintf(fid,'Tin [degC]: %4.1f    Tout [degC]: %5.1f -%5.1f\n',Tin,min(TimoMap(:,end)),max(TimoMap(:,end)));
    fprintf(fid,'Number of iterations to convergence:   %3d/%2d\n',itend,itmax);
    fprintf(fid,'Max relative error at convergence:   %8.2e\n',emaxFlow);
    fprintf(fid,'Row# MassFlow VolFlow AbsorbVel ReAbs Coll+Hose TeeSi.Dp\n');
    fprintf(fid,'     [kg/s]   [m^3/h]   [m/s]    [-]    [kPa]     [Pa]\n');
    for n=1:(length(MnewR)/2) %given the simmetry it is enough to plot 1subfield
        fprintf(fid,' %2i %7.2f  %7.3f %7.2f %8.0f %7.2f %8.0f\n', ...
                n,MnewR(n),VRIn(n)*3600,vR(n),...
                vR(n)*0.0073./viscosityGlyMixAndWat_2(x,TimoMap(n,round(NCol/2,0))),...
                DpR(n)/1000,DpTinSi(n));
    end
    fprintf(fid,' \n'); 
    fprintf(fid,' SUPPLY pipes (and pipe segments Seg.)\n');
    fprintf(fid,'Seg. Diam. Length MassFlow VolFlow Speed  Re#   fr.Darcy Fr.Loss TeeSt.Dp\n'); 
    fprintf(fid,'      [mm]   [m]   [kg/s]  [m^3/h] [m/s]  [-]     [-]     [Pa]     [Pa] \n');
    for n=1:length(VCiIn)        %given the simmetry it is enough to plot 1subfield
      for jj=1:2
        fprintf(fid,'%2i.%i  %3.0f  %5.1f %7.3f    %5.2f   %4.2f %6.0f  %6.4f %7.0f %6.0f\n', ...
          n,jj,DCi(n,jj)*1000,LCi(n,jj),MnewCi(n),VCiIn(n)*3600,vCi(n,jj),...
          ReCi(n,jj),fCi(n,jj),DpCi(n)*(jj==2),DpTinSt(n)*(jj==2));
      end
    end
    fprintf(fid,' RETURN \n');
    for n=1:length(VCoOut)
      for jj=1:2
        fprintf(fid,'%2i.%i  %3.0f  %5.1f %7.3f    %5.2f   %4.2f %6.0f  %6.4f %7.0f \n', ...
          n,jj,DCo(n,jj)*1000,LCo(n,jj),MnewCo(n),VCoOut(n)*3600,vCo(n,jj),ReCo(n,jj),fCo(n,jj),DpCo(n)*(jj==2));
      end
    end
    fprintf(fid,'\n');

    if Tmode==1
      if mod(kk,1)==0
        figure (10+kk)
        set(gcf,'DefaultAxesColorOrder',jet(Np));
        plot([0,0.5:NCol-0.5,NCol],TimoMap,'*-');
        xlabel('Collector #')
        ylabel('Fluid temperature [^oC]')
        title(['V=',num2str(Vtot_m3h_vector(kk)),'m^3/h'])
        legend(num2str((1:24)'),'Location','SouthEast')
        grid on
      end
    end
    PcR=(MoldR(1)-MoldR(12))/(MoldR(1));
end
toc
fclose('all');    % closing output text file

figure(2)
plot([1,Np],[1,1],'k--')
set(gcf,'color','white')
hlegend=legend(legendFig2{:},'Location','NorthEast');
set(hlegend,'FontSize',9);
set(gca,'FontName','Times New Roman');
set(gca,'defaultAxesFontSize',9)
xhandle=xlabel('Collector row','FontName','Times New Roman','FontSize',10);
yhandle=ylabel('Dimensionless flow rate, V'' [-]','FontName','Times New Roman','FontSize',10);
set(gca,'XLim',[1 Np])
set(gca,'XTick',1:Np);
grid on

figure(3)
plot(Vtot_m3h_vector,DpField_NoDh/1000,'r',...
     Vtot_m3h_vector,MeasuredDp/1000,'b');
legend('Field Dp (no Dh_{manometers})','Total Dp between manometers')
grid on
xlabel('Flow rate [m^3/h]')
ylabel('\Deltap [kPa]')

figure(4)
plot(Vtot_m3h_vector,RMSD_vector,'b-o');
grid on
xlabel('Flow rate [m^3/h]')
ylabel('RMSD [-]')