function [CNprop] = getCNprop()

% Compilation of basic cosmogenic nuclide properties
% Written by Jane Lund Andersen, Aarhus University

% Half-lives and decay constants
TBe = 1.387e6;
TAl = 0.705e6;
TCl = 0.301e6; %Dunai 2010 citing Holden 1990: 301 +-2 ka
CNprop.lambda_Be = log(2)/TBe;
CNprop.lambda_Al = log(2)/TAl;
CNprop.lambda_Cl = log(2)/TCl;

% Production at Sea-level high-latitude (SLHL) using Stone (2000) scaling
CNprop.PBe0 = 4.01; % Borchers et al., 2016
CNprop.PAl0 = 27.93;% Borchers et al., 2016
CNprop.PNe0 = 4.01*4.08; % 21Ne/10Be production ratio from Balco and Shuster, 2009:
                        % 4.08 +- 0.37 at /g qtz /yr, multiplied w. 10Be spal. prod

%% Other values, no longer in use:
% CNprop.pr_fm_Be = 0.005; %fast muons (of total p)
% CNprop.pr_fm_Al = 0.006;
% CNprop.pr_nmc_Be = 0.015;
% CNprop.pr_nmc_Al = 0.018;
% CNprop.Lspal = 150; %g/cm^2 %Spallation attenuation
% CNprop.Lnmc = 1500; %negative muon capture attenuation
% CNprop.Lfm = 4320; %fast muon attenuation
% CNprop.rho = 2.65; %density
% CNprop.pr_fm_Ne = 0.018; % N21 muon props (3.6% total)
% CNprop.pr_nmc_Ne = 0.018;