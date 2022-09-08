function compile_mn_data()

% This code reads sample data from an excelfile, calculates production
% parameters using 'st' scaling, and saves a data structure for use in
% bedrockMC codes.
% 
% Note 1: At present this and subsequent codes are set up to handle 10Be 
% and 26Al (mandatory), while 36Cl and 21Ne are optional
% 
% Note 2: Production rates are presently calculated using the 'St'
% scaling scheme. This is primarily due to the lack of an open source 
% framework for calculating the 'LSDn' production rates for all four 
% nuclides in a computationally consistent way. 'LSDn' scaling is 
% better, especially at high latitudes.
% 
% Input excelfile format:
% Each sample is presented with a line for each nuclide measured. 
% The first line contains sample site information (Sample name and number,
% lat (dd), lon (dd), elevation (m), topographic shielding factor (0-1), 
% sample thickness (cm) and density (g/cm3), depth to top of sample (cm), 
% year of sampling (yr), estimated timing of last deglaciation and 
% uncertainty (ka), number of nuclides measured (this determines how many 
% lines of excel-file to read before next sample occurs). 
% Each line has a nuclide identificator (1=10Be, 2=26Al, 3=36Cl,4=21Ne), 
% the cosmogenic nuclide concentration (for 36Cl the radiogenic component 
% must first be subtracted) in at/g and related uncertainty.
% 36Cl lines must have full chemistry in the following order: 1. fractional
% volumetric water-content (unitless) 2. %CO2 in bulk rock, 
% 3-12. bulk rock composition (weight% oxides)in following order:
% [Na2O,MgO,Al2O3,SiO2,P2O5,K2O,CaO,TiO2,MnO,Fe2O3], 
% 13-20. bulk rock composition (ppm): [Cl,B,Sm,Gd,U,Th,Cr,Li]
% 21-25. Target composition [K2O(%),CaO(%),TiO2(%),Fe2O3(%),Cl(ppm)]

% This script cites the following literature:
% Balco 2017: Production rate calculations for cosmic-ray-muon-produced 
% 10Be and 26Al benchmarked against geological calibration data. Quaternary
% Geochronology 39 (150-173).
% Balco & Shuster 2009: Balco, G. and Shuster, D.L., 2009. Production rate
% of cosmogenic 21Ne in quartz estimated from 10Be, 26Al, and 21Ne 
% concentrations in slowly eroding Antarctic bedrock surfaces. Earth and 
% Planetary Science Letters, 281(1-2), pp.48-58.
% Borchers et al., 2016: Borchers, B., Marrero, S., Balco, G., Caffee, M., 
% Goehring, B., Lifton, N., Nishiizumi, K., Phillips, F., Schaefer, J. and 
% Stone, J., 2016. Geological calibration of spallation production rates in 
% the CRONUS-Earth project. Quaternary Geochronology, 31, pp.188-198.
% Fernandez-Mosquera 2010: ?
% Marrero, S. M., Phillips, F. M., Borchers, B., Lifton, N., Aumer, R., & 
% Balco, G. (2016). Cosmogenic nuclide systematics and the CRONUScalc 
% program. Quaternary Geochronology, 31, 160-187.

% Codes used for calculation of production rates retrieved from:
% https://bitbucket.org/cronusearth/cronus-calc/src/master/ in july 2021
% Website: https://cronus.cosmogenicnuclides.rocks/2.1/
% See Marrero et al., 2016 for details

% Written by Jane Lund Andersen, Aarhus University, 2019-2022

close all; clear

% Add paths to functions called in this script
addpath('Functions','Functions/cl36_v2_1','Functions/CronusCalc') 

%% Read data info from Excelfile
excelfile = 'data/Input_mn.xlsx';
sheet = 'MDML1'; %'MDML1' = Heimefrontfjella, 'MDML2' = Jutulstraumen
[num,text,~] = xlsread(excelfile,sheet);

ns = max(num(:,1)); %number of samples/models in sheet

ij=1; %row number of first sample
for i=1:ns %loop through all samples
    sample{i}.Nnuc = num(ij,12); %Number of nuclides measured for sample i
    sample{i}.nuclides = num(ij:ij+sample{i}.Nnuc-1,13); %nuclide identifications for sample i
    
    %Sample specific data
    sample{i}.name = text(ij+1,1);
    sample{i}.lat=num(ij,2);
    sample{i}.lon=num(ij,3);
    sample{i}.elev=num(ij,4);
    % Calculate atm pressure (hPa) from elevation
        p = ERA40atm(sample{i}.lat,sample{i}.lon,sample{i}.elev);    
        sample{i}.pressure = p; 
    sample{i}.shield=num(ij,5); %Topographic shielding
    sample{i}.thick=num(ij,6); %Sample thickness (cm)
    sample{i}.density = num(ij,7); %(g/cm3)
    rho = sample{i}.density; %Density of rock sample (g/cm3) for internal use in function
    sample{i}.depth=num(ij,8); %surface samples = 0
    sample{i}.sampleyr = num(ij,9);
    sample{i}.minDgla = (num(ij,10)-num(ij,11))*1e-3; %last deglaciation minimum estimate -> [Myr]
    sample{i}.maxDgla = (num(ij,10)+num(ij,11))*1e-3; %last deglaciation maximum estimate -> [Myr]
    
    %% Site-specific production parameters
    
    % Define depths below surface z/rho [cm/(g/cm3)]=[g/cm^2] for fitting of production profiles
    D_m = 100; %Max depth (m), changed below
    z_m = linspace(0,10,100);
    z_D = D_m*z_m.^3/10*rho; %denser depth-grid near surface
    
    % Spallation attenuation length calculated from CronusCalc functions 
    % based on site cutoff rigidity, not considering terrain shielding [g/cm2]
    % Lspal=attenuationlength(sample{i}.lat,sample{i}.lon,sample{i}.elev,p);
    Lspal = 155; %Alternative, constant value [g/cm2]
    sample{i}.production.Lspal=Lspal; %[g/cm2]
    
    % Thickness correction
    sf_spal = exp(-sample{i}.thick/2*rho/Lspal); %Factor to correct production 
    % rate for thickness of sample, sets surface production = production midway 
    % through sample. Make sure this is not already factored in to site-specific 
    % production rates
    
    maxZ = 1200; %maxdepth (g/cm2) for muon-production profile used for fitting
    % of exponentials below. If this depth is very large, exponential terms will
    % be dominated by fast muon production, which isn't ideal. 1200 g/cm2=4.5m
    % with rho ~2.65-2.7, Test effect of this choice


%     nuclide specific data
        for j=1:sample{i}.Nnuc %loop over nuclides measured for sample i
            nuclide=num(ij+j-1,13); %get nuclide identification
            switch nuclide %calculate relevant production parameters only
                case 1 %10Be
                sample{i}.N10 = num(ij+j-1,14);
                sample{i}.dN10 = num(ij+j-1,15);
                
                % Spallation surface production in atoms / g qtz / year
                % 4.01 = Be10-spallation at SLHL from Borchers et al.,
                % 2016, Table 7, for the st scaling framework
                sample{i}.production.P10spal = 4.01.*stone2000(sample{i}.lat,p,1); 
                sample{i}.production.P10spal = sample{i}.production.P10spal*sf_spal*sample{i}.shield; %sample thickness and topographic shielding correction
    
                %Muon production following Balco 2017, 
                %Muon production parameters, from BCO fit, Model 1A, alpha=1; f_star*f_C*f_D
                mc10.k_neg = 0.00191 .* 0.704 .* 0.1828; %summary cross-section for negative muon capture (at/muon)
                mc10.sigma0 = 0.280e-30; %x-section for fast muon production at 1 Gev
                mc10.Natoms = 2.006e22; %Oxygen atoms pr gram Quartz
                
                % Fit muon production profile calculated with P_mu_total_alpha1.m 
                % with two exponential terms following Balco 2017, Eq. 7/Fig. 15
                p10_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,2,mc10,0);
                sample{i}.production.P10_Lm1=p10_muons.L(1); %attenuation, first exponential term
                sample{i}.production.P10_Lm2=p10_muons.L(2); %attenuation, second exponential term
                sample{i}.production.P10_m1 = p10_muons.P(1); %production, first exponential term
                sample{i}.production.P10_m2 = p10_muons.P(2); %production, second exponential term
                shield_fac10_m1 = exp(-sample{i}.thick/2*rho/p10_muons.L(1)); %sample thickness correction
                shield_fac10_m2 = exp(-sample{i}.thick/2*rho/p10_muons.L(2)); %sample thickness correction
                sample{i}.production.P10_m1 = sample{i}.production.P10_m1*shield_fac10_m1*sample{i}.shield; % production, first exponential term, corrected for sample thickness and topographic shielding
                sample{i}.production.P10_m2 = sample{i}.production.P10_m2*shield_fac10_m2*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding
                
                case 2 %26Al
                sample{i}.N26 = num(ij+j-1,14);
                sample{i}.dN26 = num(ij+j-1,15);
                
                % Spallation surface production in atoms / g qtz / year
                % 27.93 = Al26-spallation at SLHL from Borchers et al.,
                % 2016, Table 7, for the st scaling framework
                sample{i}.production.P26spal = 27.93.*stone2000(sample{i}.lat,p,1);
                sample{i}.production.P26spal = sample{i}.production.P26spal*sf_spal*sample{i}.shield; %sample thickness and topographic shielding correction
                
                %Muon production following Balco 2017, 
                %Muon production parameters, from BCO fit, Model 1A, alpha=1; f_star*f_C*f_D
                mc26.k_neg = 0.0133 .* 0.296 .* 0.6559; %summary cross-section for negative muon capture (at/muon)
                mc26.sigma0 = 3.89e-30; %x-section for fast muon production at 1 Gev
                mc26.Natoms = 1.003e22; %Si atoms pr gram Quartz
                              
                % Fit muon production profile calculated with P_mu_total_alpha1.m 
                % with two exponential terms following Balco 2017, Eq. 7/Fig. 15
                p26_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,2,mc26,0);
                sample{i}.production.P26_Lm1=p26_muons.L(1); %attenuation, first exponential term
                sample{i}.production.P26_Lm2=p26_muons.L(2); %attenuation, second exponential term
                sample{i}.production.P26_m1 = p26_muons.P(1); %production, first exponential term
                sample{i}.production.P26_m2 = p26_muons.P(2); %production, second exponential term
                shield_fac26_m1 = exp(-sample{i}.thick/2*rho/p26_muons.L(1)); %sample thickness correction
                shield_fac26_m2 = exp(-sample{i}.thick/2*rho/p26_muons.L(2)); %sample thickness correction
                sample{i}.production.P26_m1 = sample{i}.production.P26_m1*shield_fac26_m1*sample{i}.shield; %production, first exponential term, corrected for sample thickness and topographic shielding
                sample{i}.production.P26_m2 = sample{i}.production.P26_m2*shield_fac26_m2*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding
                
                case 3 %36Cl
                sample{i}.N36=num(ij+j-1,14);
                sample{i}.dN36=num(ij+j-1,15);
                sample{i}.chem=num(ij+j-1,16:40);
                
                %Assemble Cl36 input data and calculate production parameters
                Cl36chem = [sample{i}.N36 zeros(1,2) sample{i}.chem(1) ...
                    sample{i}.density sample{i}.thick sample{i}.lat ...
                    sample{i}.lon,sample{i}.elev p sample{i}.shield Lspal...
                    sample{i}.chem(2:end) sample{i}.depth sample{i}.sampleyr];
                scaling_model = 'st';
                %calc production rate from spallation, muons, thermal and epithermal neutrons
                Cl36prod=get36clProd(Cl36chem,z_D(z_D<1000),scaling_model); 

                %Spallation 36Cl (at/g/yr)
                sample{i}.production.P36spal=Cl36prod.P36s(1); %takes only surface value, topographic shielding correction already applied
                sample{i}.production.P36spal = sample{i}.production.P36spal*sf_spal; %apply sample thickness correction
                
                %Muons 36Cl (at/kg/yr)
                p36_muons = fit_P_exp(z_D(z_D<1000),Cl36prod.P36nm,Cl36prod.P36fm,2,0); %Fit muon production profile with two exponentials
                sample{i}.production.P36_Lm1=p36_muons.L(1);
                sample{i}.production.P36_Lm2=p36_muons.L(2);
                sample{i}.production.P36_m1 = p36_muons.P(1);
                sample{i}.production.P36_m2 = p36_muons.P(2);
                shield_fac36_m1 = exp(-sample{i}.thick/2*rho/p36_muons.L(1));
                shield_fac36_m2 = exp(-sample{i}.thick/2*rho/p36_muons.L(2));
                sample{i}.production.P36_m1 = p36_muons.P(1)*shield_fac36_m1; %sample thickness correction
                sample{i}.production.P36_m2 = p36_muons.P(2)*shield_fac36_m2; %sample thickness correction
            %     plot(Cl36prod.P36nm,z_D(z_D<1000)/rho,'.-',Cl36prod.P36fm,z_D(z_D<1000)/rho,'.-'), set(gca,'ydir','reverse'),legend('Negative muons','Fast muons')

                % Thermal neutron production 36Cl (at/kg/yr) - 4 terms 
                % Term 1-3
                sample{i}.production.P36_Lth1=Cl36prod.Lth1;
                sample{i}.production.P36_Lth2=Cl36prod.Lth2; 
                sample{i}.production.P36_Lth3=Cl36prod.Lth3; 
                shield_fac36_th1 = exp(-sample{i}.thick/2*rho/Cl36prod.Lth1);
                shield_fac36_th2 = exp(-sample{i}.thick/2*rho/Cl36prod.Lth2);
                shield_fac36_th3 = exp(-sample{i}.thick/2*rho/Cl36prod.Lth3);
                sample{i}.production.P36_th1 = Cl36prod.P36th1(1)*shield_fac36_th1; %sample thickness correction
                sample{i}.production.P36_th2 = Cl36prod.P36th2(1)*shield_fac36_th2; %sample thickness correction
                sample{i}.production.P36_th3 = Cl36prod.P36th3(1)*shield_fac36_th3; %sample thickness correction
            %     plot(Cl36prod.P36th3,z_D(z_D<1000)/rho,'.-'), set(gca,'ydir','reverse','xscale','log'),
            %     hold on, plot(Cl36prod.P36th3(1).*exp(-z_D(z_D<1000)./Cl36prod.Lth3),z_D(z_D<1000)/rho,'.-')

                %Epithermal 36Cl (at/kg/yr) - 3 terms% Term 1-2
                sample{i}.production.P36_Leth1=Cl36prod.Lth1; %same lengthscale as thermal1
                sample{i}.production.P36_Leth2=Cl36prod.Lth2;%same lengthscale as thermal2
                sample{i}.production.P36_eth1 = Cl36prod.P36eth1(1)*shield_fac36_th1; %sample thickness correction
                sample{i}.production.P36_eth2 = Cl36prod.P36eth2(1)*shield_fac36_th2; %sample thickness correction
            %     plot(Cl36prod.P36eth1,z_D(z_D<1000)/rho,'.-',Cl36prod.P36eth2,z_D(z_D<1000)/rho,'.-',Cl36prod.P36eth3,z_D(z_D<1000)/rho,'.-',Cl36prod.P36eth,z_D(z_D<1000)/rho,'k-'), set(gca,'ydir','reverse'),legend('Epithermal1','Epithermal2','Epithermal3','Epithermal'), grid on
            %     hold on, plot(Cl36prod.P36eth2(1).*exp(-z_D(z_D<1000)./Cl36prod.Lth2),z_D(z_D<1000)/rho,'.-')        

                % Non-exponential thermal+epithermal production terms, fit with two exponentials
                p36_theth = fit_P_exp(z_D(z_D<1000),Cl36prod.P36eth3+Cl36prod.P36th4,0,2,0);
                sample{i}.production.P36_Ltheth1=p36_theth.L(1);
                sample{i}.production.P36_Ltheth2=p36_theth.L(2);
                shield_fac36_theth1 = exp(-sample{i}.thick/2*rho/p36_theth.L(1));
                shield_fac36_theth2 = exp(-sample{i}.thick/2*rho/p36_theth.L(2));
                sample{i}.production.P36_theth1 = p36_theth.P(1)*shield_fac36_theth1; %sample thickness correction
                sample{i}.production.P36_theth2 = p36_theth.P(2)*shield_fac36_theth2; %sample thickness correction
            %     plot(Cl36prod.P36th1,z_D(z_D<1000)/rho,'.-',Cl36prod.P36th2,z_D(z_D<1000)/rho,'.-',Cl36prod.P36th3,z_D(z_D<1000)/rho,'.-',Cl36prod.P36th4,z_D(z_D<1000)/rho,'.-',Cl36prod.P36th,z_D(z_D<1000)/rho,'k-'), set(gca,'ydir','reverse'),legend('Thermal1','Thermal2','Thermal3','Thermal4','Thermal')
            %     plot(Cl36prod.P36th3,z_D(z_D<1000)/rho,'.-'), set(gca,'ydir','reverse','xscale','log'),%legend('Thermal1','Thermal2','Thermal3','Thermal4','Thermal')
 
                case 4 %21Ne
                sample{i}.N21=num(ij+j-1,14);
                sample{i}.dN21=num(ij+j-1,15);
                
                % 21Ne/10Be production ratio from Balco and Shuster, 2009:
                % 4.08 +- 0.37 at /g qtz /yr, multiplied w. 10Be spal. prod
                sample{i}.production.P21spal = 4.08*4.01.*stone2000(sample{i}.lat,p,1);
                sample{i}.production.P21spal = sample{i}.production.P21spal*sf_spal*sample{i}.shield; %sample thickness and topographic shielding correction
                
                % Muon production following Balco 2017, 
                % Muon production parameters, from BCO fit, Model 1A, alpha=1; f_star*f_C*f_D
                % Muons 21Ne. These are not calibrated, taken from Fernandez-Mosquera 2010
                mc21.Natoms = 1.0003e22; %Si atoms pr gram Quartz
                mc21.k_neg = 0.296.*0.6559.*0.0029; %summary cross-section for negative muon capture (at/muon)
                mc21.sigma190 = 0.79e-27; %x-section for fast muon production at 1 Gev
                
                % Fit muon production profile calculated with P_mu_total_alpha1.m 
                % with two exponential terms following Balco 2017, Eq. 7/Fig. 15
                p21_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,2,mc21,0);
                sample{i}.production.P21_Lm1=p21_muons.L(1); %attenuation, first exponential term
                sample{i}.production.P21_Lm2=p21_muons.L(2); %attenuation, second exponential term
                sample{i}.production.P21_m1 = p21_muons.P(1); %production, first exponential term
                sample{i}.production.P21_m2 = p21_muons.P(2); %production, second exponential term
                shield_fac21_m1 = exp(-sample{i}.thick/2*rho/p21_muons.L(1)); %sample thickness correction
                shield_fac21_m2 = exp(-sample{i}.thick/2*rho/p21_muons.L(2)); %sample thickness correction
                sample{i}.production.P21_m1 = sample{i}.production.P21_m1*shield_fac21_m1*sample{i}.shield; %production, first exponential term, corrected for sample thickness and topographic shielding 
                sample{i}.production.P21_m2 = sample{i}.production.P21_m2*shield_fac21_m2*sample{i}.shield; %production, second exponential term, corrected for sample thickness and topographic shielding 
                
            end
        end
        
        %% Calculate nuclide ratios
        if isfield(sample{i},'N26')
            sample{i}.r2610 = sample{i}.N26./sample{i}.N10;
            sample{i}.dr2610 = sample{i}.r2610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN26./sample{i}.N26).^2);
        end       
        if isfield(sample{i},'N21')
            sample{i}.r1021 = sample{i}.N10./sample{i}.N21;
            sample{i}.dr1021 = sample{i}.r1021.*sqrt((sample{i}.dN21./sample{i}.N21).^2+(sample{i}.dN10./sample{i}.N10).^2);
            sample{i}.r2621 = sample{i}.N26./sample{i}.N21;
            sample{i}.dr2621 = sample{i}.r2621.*sqrt((sample{i}.dN21./sample{i}.N21).^2+(sample{i}.dN10./sample{i}.N10).^2);
        end
        if isfield(sample{i},'N36')
            sample{i}.r3610 = sample{i}.N36./sample{i}.N10;
            sample{i}.dr3610 = sample{i}.r3610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN36./sample{i}.N36).^2);
            sample{i}.r3626 = sample{i}.N36./sample{i}.N26;
            sample{i}.dr3626 = sample{i}.r3626.*sqrt((sample{i}.dN26./sample{i}.N26).^2+(sample{i}.dN36./sample{i}.N36).^2);
        end
        
    ij=ij+sample{i}.Nnuc; %Look for next sample in this row
end


save(['./data/',sheet,'_mn_data.mat'],'sample','excelfile')