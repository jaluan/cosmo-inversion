function [gm] = forward_bedrockvJ1_mn4_BeAlNe(up,model,CNprop)

% function [gm] = forward_bedrockvJ1_mn4_(up,model,CNprop)
% Forward model for transient CN integration of bedrock samples 

% Inputs: 
%   up = model parameter vector (generated from bedrockMCvJ1_mn5.m), 
%   model = model setup (generated from bedrockMCvJ1_mn5.m), 
%   CNprop = structure with halflives and decay constants of radioactive 
%       nuclides (from getCNprop.m)

% Output: gm = predicted data vector

% Written by David Lundbek Egholm
% Modified by Jane Lund Andersen to also i) include [36Cl, 21Ne],
% ii) handle alternate climate inputs, and iii) add an (optional) slope or 
% step to the glaciation threshold (emergence)

%% model parameters

%generic glaciation parameters
d18Op = up(1);
T1 = up(2);
    
%% load and prepare climate curve data
dO = model.dO;
dOt = model.dOt;

if model.emergence == 1
 d18Osl=up(3); %slope of d18O threshold, permille/Myr
 dO=dO-d18Osl*dOt*1e-6;
elseif model.emergence == 2
 d18O2=up(3); %2nd d18O threshold, [permille]
 dO(dOt > model.dTTh*1e6)=dO(dOt > model.dTTh*1e6)-(d18O2-d18Op); %dTTh: timing of change [Myr]
end

%time steps 
model.dt = 1000; %yr

%% loop samples
for i=1:model.Nsnr
    n0 = (i-1)*model.Nsmp+model.Mmp; %parameter number start
    
    % Unpack exhumation parameters and pack time-depth vectors
    z1 = up(n0+1); %m
    dT2 = up(n0+2); %Myr
    dz2 = 10^up(n0+3); %m
    dT3 = up(n0+4); %Myr
    dz3 = 10^up(n0+5); %m
    dT4 = up(n0+6); %Myr
    dz4 = 10^up(n0+7); %m
    E5 = 10^up(n0+8); %m/Myr

    T2 = T1 + dT2; %Ma
    z2 = z1 + dz2; %m
    
    T3 = T2 + dT3; %Ma
    z3 = z2 + dz3; %m
    
    T4 = T3 + dT4; %Ma
    z4 = z3 + dz4; %m
    
    T5 = model.age; %this requires age > Tdg+dT2+dT3+dT4
    z5 = z4 + (model.age - T4)*E5; %m
    
    mT = [0,T1,T2,T3,T4,T5]*1e6; %Ma to a
    mz = [0,z1,z2,z3,z4,z5]; %m
    
    % Check starting condition of model and set flag ssbc
    if (z5 > model.z0) % Depth at start of model greater than max depth
        maxtime = interp1(mz,mT,model.z0);
        ssbc = 0; %starting concentrations set to 0 below
    else
        maxtime = model.age*1e6;
        ssbc = 1; %starting concentrations set from steady state erosion rate below
    end
    
    nt = ceil(maxtime/model.dt); %number of time steps
    time = maxtime*linspace(0,1,nt); %time vector
    
    
    %**************************************
    % Exhumation and exposure history
    %**************************************
    burial = interp1(mT,mz,time); %Depths at times in 'time'
    dOn = interp1(dOt,dO,time,'linear',dO(end)); %climate curve values at times in 'time'
    pfac = ones(nt,1);     %controls exposure, modified below
    pfac(dOn > d18Op) = 0; %No exposure when climate curve values above threshold (d18Op)
    pfac(time < 25e3) = 0; %correct exposure around deglaciation: set no
                           %exposure last 25 kyr
    pfac(time < T1*1e6) = 1; %then add exposure from time of deglaciation
    
    % plot(time,dOn), hold on, plot(time,pfac+4), plot(time,d18Op*ones(size(time)))
    % return

    %CN production parameters
    rho = model.data{i}.density;
    Lspal = model.data{i}.production.Lspal;
    
%     if ismember(1,model.data{i}.nuclides) %10Be, Note: ismember is very slow
        P10spal = model.data{i}.production.P10spal;
        P10m1 = model.data{i}.production.P10_m1;
        P10Lm1 = model.data{i}.production.P10_Lm1;
        P10m2 = model.data{i}.production.P10_m2;
        P10Lm2 = model.data{i}.production.P10_Lm2;
%     end
    
    
%     if ismember(2,model.data{i}.nuclides) %26Al, Note: ismember is very slow
        P26spal = model.data{i}.production.P26spal;
        P26m1 = model.data{i}.production.P26_m1;
        P26Lm1 = model.data{i}.production.P26_Lm1;
        P26m2 = model.data{i}.production.P26_m2;
        P26Lm2 = model.data{i}.production.P26_Lm2;
%     end
    
    
% %     if ismember(3,model.data{i}.nuclides) %36Cl, Note: ismember is very slow
%         P36spal = model.data{i}.production.P36spal;
%         P36m1 = model.data{i}.production.P36_m1;
%         P36Lm1 = model.data{i}.production.P36_Lm1;
%         P36m2 = model.data{i}.production.P36_m2;
%         P36Lm2 = model.data{i}.production.P36_Lm2;
%         P36th1 = model.data{i}.production.P36_th1;
%         P36Lth1 = model.data{i}.production.P36_Lth1;
%         P36th2 = model.data{i}.production.P36_th2;
%         P36Lth2 = model.data{i}.production.P36_Lth2;
%         P36th3 = model.data{i}.production.P36_th3;
%         P36Lth3 = model.data{i}.production.P36_Lth3;
%         P36eth1 = model.data{i}.production.P36_eth1;
%         P36Leth1 = model.data{i}.production.P36_Leth1;
%         P36eth2 = model.data{i}.production.P36_eth2;
%         P36Leth2 = model.data{i}.production.P36_Leth2;
%         P36theth1 = model.data{i}.production.P36_theth1;
%         P36Ltheth1 = model.data{i}.production.P36_Ltheth1;
%         P36theth2 = model.data{i}.production.P36_theth2;
%         P36Ltheth2 = model.data{i}.production.P36_Ltheth2;
%     end
%     
%     
%     if ismember(4,model.data{i}.nuclides) %21Ne, Note: ismember is very slow
        P21spal = model.data{i}.production.P21spal;
        P21m1 = model.data{i}.production.P21_m1;
        P21Lm1 = model.data{i}.production.P21_Lm1;
        P21m2 = model.data{i}.production.P21_m2;
        P21Lm2 = model.data{i}.production.P21_Lm2;
%     end
%     
    
    N10 = zeros(nt,1);
    N26 = zeros(nt,1);
%     N36 = zeros(nt,1);
    N21 = zeros(nt,1);
    
    if (ssbc == 0) %depth of sample at model start greater than maxdepth
        
        N10(nt) = 0.0;
        N26(nt) = 0.0;
%         N36(nt) = 0.0;
        N21(nt) = 0.0;

    else %assume steady state concentration at starting point
        
        erate = E5*1e-6;
        
%         if ismember(1,model.data{i}.nuclides) %10Be, Note: ismember is very slow
            fBe = CNprop.lambda_Be + rho*erate*100/Lspal; %spallation
            N10(nt) = P10spal*exp(-rho*100*burial(nt)/Lspal)/fBe;
            fBe = CNprop.lambda_Be + rho*erate*100/P10Lm2; %1st muon pathway
            N10(nt) = N10(nt) + P10m2*exp(-rho*100*burial(nt)/P10Lm2)/fBe;
            fBe = CNprop.lambda_Be + rho*erate*100/P10Lm1; %2nd muon pathway
            N10(nt) = N10(nt) + P10m1*exp(-rho*100*burial(nt)/P10Lm1)/fBe;
%         end
        
%         if ismember(2,model.data{i}.nuclides) %26Al, Note: ismember is very slow
            fAl = CNprop.lambda_Al + rho*erate*100/Lspal; %spallation
            N26(nt) = P26spal*exp(-rho*100*burial(nt)/Lspal)/fAl;
            fAl = CNprop.lambda_Al + rho*erate*100/P26Lm2; %1st muon pathway
            N26(nt) = N26(nt) + P26m2*exp(-rho*100*burial(nt)/P26Lm2)/fAl;
            fAl = CNprop.lambda_Al + rho*erate*100/P26Lm1; %2nd muon pathway
            N26(nt) = N26(nt) + P26m1*exp(-rho*100*burial(nt)/P26Lm1)/fAl;
%         end
        
% %         if ismember(3,model.data{i}.nuclides) %36Cl, Note: ismember is very slow
%             fCl = CNprop.lambda_Cl + rho*erate*100/Lspal; %spallation
%             N36(nt) = P36spal*exp(-rho*100*burial(nt)/Lspal)/fCl;
%             fCl = CNprop.lambda_Cl + rho*erate*100/P36Lm1; %1st muon pathway
%             N36(nt) = N36(nt) + P36m1*exp(-rho*100*burial(nt)/P36Lm1)/fCl;
%             fCl = CNprop.lambda_Cl + rho*erate*100/P36Lm2; %2nd muon pathway
%             N36(nt) = N36(nt) + P36m2*exp(-rho*100*burial(nt)/P36Lm2)/fCl;
%             %thermal and epithermal neutrons (36Cl)
%             fCl = CNprop.lambda_Cl + rho*erate*100/P36Lth1;
%             N36(nt) = N36(nt) + P36th1*exp(-rho*100*burial(nt)/P36Lth1)/fCl;
%             fCl = CNprop.lambda_Cl + rho*erate*100/P36Lth2;
%             N36(nt) = N36(nt) + P36th2*exp(-rho*100*burial(nt)/P36Lth2)/fCl;
%             fCl = CNprop.lambda_Cl + rho*erate*100/P36Lth3;
%             N36(nt) = N36(nt) + P36th3*exp(-rho*100*burial(nt)/P36Lth3)/fCl;
%             fCl = CNprop.lambda_Cl + rho*erate*100/P36Leth1;
%             N36(nt) = N36(nt) + P36eth1*exp(-rho*100*burial(nt)/P36Leth1)/fCl;
%             fCl = CNprop.lambda_Cl + rho*erate*100/P36Leth2;
%             N36(nt) = N36(nt) + P36eth2*exp(-rho*100*burial(nt)/P36Leth2)/fCl;
%             fCl = CNprop.lambda_Cl + rho*erate*100/P36Ltheth1;
%             N36(nt) = N36(nt) + P36theth1*exp(-rho*100*burial(nt)/P36Ltheth1)/fCl;
%             fCl = CNprop.lambda_Cl + rho*erate*100/P36Ltheth2;
%             N36(nt) = N36(nt) + P36theth2*exp(-rho*100*burial(nt)/P36Ltheth2)/fCl;
%         end
    
%         if ismember(4,model.data{i}.nuclides) %21Ne, Note: ismember is very slow
            fNe = rho*erate*100/Lspal; %spallation
            N21(nt) = P21spal*exp(-rho*100*burial(nt)/Lspal)/fNe;
            fNe = rho*erate*100/P21Lm2; %1st muon pathway
            N21(nt) = N21(nt) + P21m2*exp(-rho*100*burial(nt)/P21Lm2)/fNe;
            fNe = rho*erate*100/P21Lm1; %2nd muon pathway
            N21(nt) = N21(nt) + P21m1*exp(-rho*100*burial(nt)/P21Lm1)/fNe;
%         end
    end
    
        
    %integrate time
    for kk=(nt-1):-1:1
        
        dt = (time(kk+1)-time(kk)); %length of timestep
        pf = .5*(pfac(kk+1)+pfac(kk)); %exposure: 1=full exposure, 0=no exposure (glacial cover)
        bz = burial(kk); %depth at time kk
        erate = 100*(burial(kk+1)-burial(kk))/dt; %erosion rate within timestep
        
        % Be-10
%         if ismember(1,model.data{i}.nuclides) %10Be, Note: ismember is very slow
            %Calculate production profiles, spallation and muon pathways
            P10z_spal = pf*P10spal*exp(-rho*100*bz/Lspal);
            P10z_m1 = pf*P10m2*exp(-rho*100*bz/P10Lm2);
            P10z_m2 = pf*P10m1*exp(-rho*100*bz/P10Lm1);

            N10(kk) = N10(kk+1)*exp(-dt*CNprop.lambda_Be); %decay from previous step

            ff = CNprop.lambda_Be+rho*erate/Lspal; 
            N10(kk) = N10(kk) + P10z_spal*(1.0-exp(-ff*dt))/ff; %add spallation production (Dunai 2010 eq. 4.10)

            ff = CNprop.lambda_Be+rho*erate/P10Lm2;
            N10(kk) = N10(kk) + P10z_m1*(1.0-exp(-ff*dt))/ff; %add 1st muon pathway

            ff = CNprop.lambda_Be+rho*erate/P10Lm1;
            N10(kk) = N10(kk) + P10z_m2*(1.0-exp(-ff*dt))/ff; %add 2nd muon pathway
%         end
        
        % Al-26
%         if ismember(2,model.data{i}.nuclides) %26Al, Note: ismember is very slow
            %Calculate production profiles, spallation and muon pathways
            P26z_spal = pf*P26spal*exp(-rho*100*bz/Lspal);
            P26z_m1 = pf*P26m2*exp(-rho*100*bz/P26Lm2);
            P26z_m2 = pf*P26m1*exp(-rho*100*bz/P26Lm1);

            N26(kk) = N26(kk+1)*exp(-dt*CNprop.lambda_Al); %decay from previous step

            ff = CNprop.lambda_Al+rho*erate/Lspal;
            N26(kk) = N26(kk) + P26z_spal*(1.0-exp(-ff*dt))/ff; %add spallation production

            ff = CNprop.lambda_Al+rho*erate/P26Lm2;
            N26(kk) = N26(kk) + P26z_m1*(1.0-exp(-ff*dt))/ff; %add 1st muon pathway

            ff = CNprop.lambda_Al+rho*erate/P26Lm1;
            N26(kk) = N26(kk) + P26z_m2*(1.0-exp(-ff*dt))/ff; %add 2nd muon pathway
%         end
        
%         % Cl-36   
%         if ismember(3,model.data{i}.nuclides) %36Cl, Note: ismember is very slow
%             %Calculate production profiles, spallation, muon, thermal and epithermal pathways
%             P36z_spal = pf*P36spal*exp(-rho*100*bz/Lspal);
%             P36z_m1 = pf*P36m1*exp(-rho*100*bz/P36Lm1);
%             P36z_m2 = pf*P36m2*exp(-rho*100*bz/P36Lm2);
%             P36z_th1 = pf*P36th1*exp(-rho*100*bz/P36Lth1);
%             P36z_th2 = pf*P36th2*exp(-rho*100*bz/P36Lth2);
%             P36z_th3 = pf*P36th3*exp(-rho*100*bz/P36Lth3);
%             P36z_eth1 = pf*P36eth1*exp(-rho*100*bz/P36Leth1);
%             P36z_eth2 = pf*P36eth2*exp(-rho*100*bz/P36Leth2);
%             P36z_theth1 = pf*P36theth1*exp(-rho*100*bz/P36Ltheth1);
%             P36z_theth2 = pf*P36theth2*exp(-rho*100*bz/P36Ltheth2);
% 
%             N36(kk) = N36(kk+1)*exp(-dt*CNprop.lambda_Cl); %decay from previous step
% 
%             ff = CNprop.lambda_Cl+rho*erate/Lspal;
%             N36(kk) = N36(kk) + P36z_spal*(1.0-exp(-ff*dt))/ff; %add spallation production
% 
%             ff = CNprop.lambda_Cl+rho*erate/P36Lm1;
%             N36(kk) = N36(kk) + P36z_m1*(1.0-exp(-ff*dt))/ff; %add 1st muon pathway
% 
%             ff = CNprop.lambda_Cl+rho*erate/P36Lm2;
%             N36(kk) = N36(kk) + P36z_m2*(1.0-exp(-ff*dt))/ff; %add 2nd muon pathway
% 
%             ff = CNprop.lambda_Cl+rho*erate/P36Lth1;
%             N36(kk) = N36(kk) + P36z_th1*(1.0-exp(-ff*dt))/ff; %add thermal 1 production
% 
%             ff = CNprop.lambda_Cl+rho*erate/P36Lth2;
%             N36(kk) = N36(kk) + P36z_th2*(1.0-exp(-ff*dt))/ff; %add thermal 2 production
% 
%             ff = CNprop.lambda_Cl+rho*erate/P36Lth3;
%             N36(kk) = N36(kk) + P36z_th3*(1.0-exp(-ff*dt))/ff; %add thermal 3 production
% 
%             ff = CNprop.lambda_Cl+rho*erate/P36Leth1;
%             N36(kk) = N36(kk) + P36z_eth1*(1.0-exp(-ff*dt))/ff; %add epithermal 1 production
% 
%             ff = CNprop.lambda_Cl+rho*erate/P36Leth2;
%             N36(kk) = N36(kk) + P36z_eth2*(1.0-exp(-ff*dt))/ff; %add epithermal 2 production
% 
%             ff = CNprop.lambda_Cl+rho*erate/P36Ltheth1;
%             N36(kk) = N36(kk) + P36z_theth1*(1.0-exp(-ff*dt))/ff; %add thermal-epithermal 1 production
% 
%             ff = CNprop.lambda_Cl+rho*erate/P36Ltheth2;
%             N36(kk) = N36(kk) + P36z_theth2*(1.0-exp(-ff*dt))/ff; %add thermal-epithermal 2 production
%         end
%         
%         % Ne-21  
%         if ismember(4,model.data{i}.nuclides) %21Ne, Note: ismember is very slow
%             %Calculate production profiles, spallation and muon pathways
            P21z_spal = pf*P21spal*exp(-rho*100*bz/Lspal);
            P21z_m1 = pf*P21m1*exp(-rho*100*bz/P21Lm1);
            P21z_m2 = pf*P21m2*exp(-rho*100*bz/P21Lm2);

            N21(kk) = N21(kk+1); %No decay - 21Ne is stable

            ff = rho*erate/Lspal;
            N21(kk) = N21(kk) + P21z_spal*(1.0-exp(-ff*dt))/ff; %add spallation production

            ff = rho*erate/P21Lm2;
            N21(kk) = N21(kk) + P21z_m1*(1.0-exp(-ff*dt))/ff; %add 1st muon pathway

            ff = rho*erate/P21Lm1;
            N21(kk) = N21(kk) + P21z_m2*(1.0-exp(-ff*dt))/ff; %add 2nd muon pathway
%         end
    end

    % ouput modelled data    
    for j=1:model.data{i}.Nnuc %loop nuclides
        nuclide=model.data{i}.nuclides(j); %get nuclide identification
            if nuclide == 1 %10Be
                gm((i-1)*model.Nds+j) = N10(1);
            elseif nuclide == 2 %26Al
                gm((i-1)*model.Nds+j) = N26(1);
%             elseif nuclide == 3 %36Cl
%                 gm((i-1)*model.Nds+j) = N36(1);
            elseif nuclide == 4 %21Ne
                gm((i-1)*model.Nds+j) = N21(1);
            end
    end

end