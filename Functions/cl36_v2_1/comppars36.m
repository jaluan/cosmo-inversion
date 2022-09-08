%
% cp=comppars36(physpars,samppars,scalefactors,maxdepth)
%
% Generates a complete set of parameters for Chlorine 36 dating
% from the basic physics parameters, sample specific parameters, 
% and scale factors.
%
% maxdepth is an optional parameter.  The default value is 2500
% g/cm^2, or about 10 meters.  Computing muon production profiles
% down this far is extremely expensive.  For surface samples, a max
% depth of 25 g/cm^2 (or about 10cm)  is more appropriate.  
%
function cp=comppars36(pp,sp36,sf36,maxdepth)
%
% First, set maxdepth to a default value if not supplied.
%
if (nargin < 4)
  maxdepth=2500;
end
%
% Either way, keep track of maxdepth so we can make sure we're safe
% in prodz().
%
cp.maxdepth=maxdepth;
%
% Setup Lambdafe.  
%
cp.Lambdafe=sp36.Lambdafe;
cp.ls=sp36.ls;
cp.rb=sp36.rb;
%
% Ni are the atoms per gram for C, Na, Mg, Al, Si, P, K, Ca, Ti,
% Mn, Fe, followed by Cl, B, Sm, Gd, U, Th, Cr, Li.  Note that C starts in 
% entry 3 of the array to leave room for H(2) and O(1). 
% Target elements will be calculated in cp.tni in a similar fashion.
% Target elements in order are K, Ca, Ti, Fe, Cl.  
%
cp.ni=zeros(21,1);
for i=3:13
  cp.ni(i)=sp36.ci(i)*pp.table(i,9)*1.0e19;
end
for i=14:21
  cp.ni(i)=sp36.ci(i)*pp.table(i,9)*1.0e15;
end

cp.ni(2)=sp36.qavg*0.111*pp.ag/sp36.rb;
cp.ni(1)=cp.ni(2:13)'*pp.table(2:13,10);

cp.tni=zeros(5,1);
for i=1:3
  cp.tni(i)=sp36.tci(i)*pp.table((i+8),9)*1.0e19;
end
cp.tni(4)=sp36.tci(4)*pp.table(13,9)*1.0e19;
cp.tni(5)=sp36.tci(5)*pp.table(14,9)*1.0e15;
%
% Fraction of O, H, C, Na, Mg, etc.
%
cp.fi=zeros(21,1);
cp.fi=cp.ni .* pp.table(:,1)/ ...
    pp.ag;
% Do I need to add something here?
%
% Ass is the average atomic weight of the subsurface material.
%
cp.Ass=pp.ag/sum(cp.ni);
%
% Sigma_sc,ss=sumproduct(Ni,ssc_i)/1.0e24
% Eqn 3.22? in Gosse & Phillips
%
cp.Sigmascss=cp.ni'*pp.table(:,3)/ ...
    1.0e24;
%
% Sigma_th,ss=sumproduct(Ni,ssth_i)/1.0e24
% Eqn 3.6 in Gosse & Phillips
%

cp.Sigmathss=cp.ni'*pp.table(:,4)/ ...
    1.0e24;
%
% I_eff=sumproduct(Ni,Ia_i)/1.0e24;
% I is the effective resonance integral for 
% absorption of epithermal neutrons
% Eqn 3.9 in G&P
%
cp.Ieff=cp.ni'*pp.table(:,5)/1.0e24;
%
% X_ss=sumproduct(Ni,xi,ssc_i)/sumproduct(Ni,ssc_i)
%
cp.Xss=...
    sum((cp.ni.*pp.table(:,2)).*pp.table(:,3))/...
    (cp.ni'*pp.table(:,3));
%
% Sigma_eth,ss=Xss*(Ieff+Sigma_sc,ss)
%
cp.Sigmaethss=cp.Xss*...
    (cp.Ieff+cp.Sigmascss);
%
% Lambda_f,t=(4.3/3.3)*Lambda_f,a
%
%cp.Lambdaft=(4.3/3.3)*pp.Lambdafe;
%
% Lambda_eth,a=1/Sigma_eth,a
%
cp.Lambdaetha=1/pp.Sigmaetha;
%
% Lambda_eth,ss=1/Sigma_eth,ss
%
cp.Lambdaethss=1/cp.Sigmaethss;
%
% Lambda_th,ss=1/Sigma_th,ss
%
cp.Lambdathss=1/cp.Sigmathss;
%
% Added new variables for individual spallation rxns
%
cp.PsCa=cp.tni(2)*pp.PsCa0*(40.078/6.022E+23); % () is a conversion to at 36Cl/at Ca;
cp.PsK=cp.tni(1)*pp.PsK0*(39.0983/6.022E+23); % () is a conversion to at 36Cl/at K;
cp.PsTi=cp.tni(3)*pp.PsTi0;
cp.PsFe=cp.tni(4)*pp.PsFe0;
%
% B=sumproduct(xi,Ni,ssc_i)/1.0e24;
% D34 on Parameters page in CHLOE
% B is the scattering rate parameter
%
cp.B=...
    sum((cp.ni.*pp.table(:,2)).*pp.table(:,3))/1.0e24;
%
% p(Eth)_ss=exp(-Ieff/B)
%
cp.pEthss=exp(-cp.Ieff/cp.B);
%
% Reth=sqrt(Ass/14.5)
% Reth is the production rate of epithermal neutrons from 
% fast neutrons in the air at land/atmosphere boundary
% P. 1498 in G&P
%
cp.Reth=sqrt(cp.Ass/14.5);
%
% Rth=pEthss/pEtha
% Eqn 3.20 in G&P
%
cp.Rth=cp.pEthss/pp.pEtha;
%
% D_eth,a=1/(3*Sigma_sc,a*(1-2/(3*Aa)))
% Eqn 3.16 in G&P
%
cp.Detha=1/(3*pp.Sigmasca*...
                        (1-2/(3*pp.Aa)));
cp.Da=cp.Detha;
%
% D_th,ss=(3*Sigma_sc,ss*(1-2/(3*Ass)))^(-1)
% Eqn 3.33 in G&P
%
cp.Dthss=(3*cp.Sigmascss*(1-2/(3*cp.Ass)))^(-1);
cp.Dethss=cp.Dthss;
%
% Lambda_th,a=1/Sigma_th,a
%
cp.Lambdatha=1/pp.Sigmatha;
%
% L_eth,a=1/sqrt(3*SSc_a*Seth_a)
% epithermal diffusion length (eqn 3.21? in G&P)
cp.Letha=1/sqrt(3*pp.Sigmasca* ...
                            pp.Sigmaetha);
%
% L_eth,ss=1/sqrt(Sigma_eth,ss*(3*Sigma_sc,ss))
%
cp.Lethss=1/sqrt(cp.Sigmaethss*...
                           (3*cp.Sigmascss));
%
% L_th,ss=sqrt(D_th,ss/Sigma_th,ss)
% thermal neutron diffusion length (eqn 3.34? in G&P)
%
cp.Lthss=sqrt(cp.Dthss/ ...
                          cp.Sigmathss);
%
% f_th=((s_th,cl*Ncl)/Sigma_th,ss)/1.0e24;
% Eqn 3.32 in G&P
%
cp.fth=((pp.table(14,4)*cp.tni(5))/...
    (cp.Sigmathss))/1.0e24;
%
% f_eth=((la_cl*Ncl)/Ieff)/1.0e24;
%
cp.feth=((pp.table(14,5)*cp.tni(5))/...
    cp.Ieff)/1.0e24;
%
% Phistaretha=Pf0/(Sigmaetha-Detha/Lambdafe^2)
% Eqn 3.26 in G&P
%
cp.Phistaretha=pp.Pf0/...
    (pp.Sigmaetha-cp.Detha/ ...
     cp.Lambdafe^2);
%
% Phistartha=p(Eth)_a*Sigma_eth,a*PhiStar_eth,a/(Sigma_th,a-Da/Lf_a^2)
% Eqn 3.38 in G&P
%
cp.Phistartha=pp.pEtha*pp.Sigmaetha*...
    cp.Phistaretha/...
    (pp.Sigmatha-cp.Da/cp.Lambdafe^2);
%
% Phistarethss=Pf0*Reth/(Sigmaethss-Dethss/Lambda_fa^2)
% Eqn 3.26 in G&P
%
cp.Phistarethss=pp.Pf0*cp.Reth/...
    (cp.Sigmaethss-cp.Dethss/ ...
     cp.Lambdafe^2);
%
% Phistarthss=pEtha*Sigmaethss*PhiStarethss*Rth/(Sigmathss-Dthss/Lambda_fa^2);
% Eqn 3.38 in G&P
%
cp.Phistarthss=pp.pEtha*cp.Sigmaethss*...
    cp.Phistarethss*cp.Rth/...
    (cp.Sigmathss-cp.Dthss/cp.Lambdafe^2);
%
% DeltaPhistareth=Phistaretha-Phistarethss
% Eqn 3.29 in G&P
%
cp.DeltaPhistareth=cp.Phistaretha-...
    cp.Phistarethss;
cp.DeltaPhistarethss=cp.DeltaPhistareth;
cp.DeltaPhistaretha=-cp.DeltaPhistareth;
%
% DeltaPhistarth=Phistartha-Phistarthss
% Eqn 3.3.41 in G&P
%
cp.DeltaPhistarth=cp.Phistartha-...
    cp.Phistarthss;
cp.DeltaPhistarthss=cp.DeltaPhistarth;
cp.DeltaPhistartha=-cp.DeltaPhistarth;
%
% DeltaPhistarstareth=Phistarethss-Detha*Phistaretha/Dethss
% Eqn 3.30 in G&P
%
cp.DeltaPhistarstareth=cp.Phistarethss-...
    cp.Detha*cp.Phistaretha/ ...
    cp.Dethss;
%
% FDeltaPhistarethss=((Detha*DelPhiStareth/Letha)-
%                     (Dethss*DelPhiDblStareth/Lfa))/
%                     ((Detha/Letha)+(Dethss/Lethss))
% Eqn 3.28 in G&P
%
cp.FDeltaPhistarethss=...
    (cp.Detha*cp.DeltaPhistareth/cp.Letha-...
     (cp.Dethss*cp.DeltaPhistarstareth/cp.Lambdafe))/...
    ((cp.Detha/cp.Letha)+...
     (cp.Dethss/cp.Lethss));
cp.FDeltaPhistareth=cp.FDeltaPhistarethss;
%
% FDeltaPhistaretha=((Dethss*(-DeltaPhistareth)/Lethss)-(Dethss*DeltaPhistarstareth/Lambdafe))/...
%  ((Detha/Letha)+Dethss/Lethss))
% Eqn 3.28 in G&P
%
cp.FDeltaPhistaretha=...
    ((cp.Dethss*(-cp.DeltaPhistareth)/cp.Lethss)-...
     (cp.Dethss*cp.DeltaPhistarstareth/cp.Lambdafe))/...
    ((cp.Detha/cp.Letha)+...
     (cp.Dethss/cp.Lethss));
%
% Eqn 3.39 in G&P
%
cp.SFDeltaPhistaretha=pp.pEtha*pp.Sigmaetha*cp.FDeltaPhistaretha/...
    (pp.Sigmatha-cp.Da/cp.Letha^2);
%
%Eqn 3.39 in G&P
%
cp.SFDeltaPhistarethss=cp.Sigmaethss*pp.pEtha*cp.FDeltaPhistareth*cp.Rth/...
    (cp.Sigmathss-cp.Dthss/cp.Lethss^2);
%
% Eqn 3.42 in G&P
%
cp.DeltaSFDeltaPhistareth=cp.SFDeltaPhistaretha- ...
    cp.SFDeltaPhistarethss;
cp.DeltaSFDeltaPhistaretha=-cp.DeltaSFDeltaPhistareth;
cp.DeltaSFDeltaPhistarethss=cp.DeltaSFDeltaPhistareth;
%
% cp.SFDeltaPhistartha.  Note correction to Eqn 3.40 in G&P
%
%
% Corrected version used below.
%
cp.SFDeltaPhistartha=...
    (cp.Da*(cp.Phistartha/cp.Lambdafe-cp.SFDeltaPhistaretha/cp.Letha)-...
     cp.Dthss*(cp.Phistarthss/cp.Lambdafe+cp.SFDeltaPhistarethss/cp.Lethss)+...
     (cp.Dthss/cp.Lthss)*(-cp.DeltaPhistarth-cp.DeltaSFDeltaPhistareth))/...
    (cp.Dthss/cp.Lthss+cp.Da/pp.La);
%
% SFDeltaPhistarthss.
%
%
% Corrected version of G&P eqn 3.40 used below.
%
cp.SFDeltaPhistarthss=...
    (cp.Da*(cp.Phistartha/cp.Lambdafe-cp.SFDeltaPhistaretha/cp.Letha)-...
     cp.Dthss*(cp.Phistarthss/cp.Lambdafe+cp.SFDeltaPhistarethss/cp.Lethss)+...
     (cp.Da/pp.La)*(cp.DeltaPhistarthss+cp.DeltaSFDeltaPhistarethss))/...
    (cp.Dthss/cp.Lthss+cp.Da/pp.La);
%
% Zs
%
%cp.Zs=sp36.ls*sp36.rb;
%
%---------------------------------------------------------------------
%
% The following are modifications for Sato/Heisinger muons.  
%
%
% Pmunf
%
% old formulation: 
%  cp.Pmunf=0.0000058*pp.Yf*pp.Phimuf0*exp(-cp.Zs/pp.Lambdamu);
%
%  New Formulation uses Greg's Heisinger code to calculate the fluxes at 
% a vector of depths to create a table that can be referred to later
% Actually get the full data from P_mu_total -- needed later
    % Use middle of sample
%     mu = P_mu_total((sample.thick.*sample.rho./2),sample.pressure,mconsts,'yes');
% RcEst = 14.9.*((cos(abs(d2r(sp36.latitude)))).^4); %Dipolar estimate for these purposes
% The RcEst above was originally used and has been replaced by Nat Lifton's
% new model.  He  fit the function below to trajectory-traced cutoffs for the
% modern dipolar field strength, so it should be more accurate overall.

dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];

%   Trajectory-traced dipolar estimate for these purposes
   RcEst = (dd(1)*cos(d2r(sp36.latitude)) + ...
       dd(2)*(cos(d2r(sp36.latitude))).^2 + ...
       dd(3)*(cos(d2r(sp36.latitude))).^3 + ...
       dd(4)*(cos(d2r(sp36.latitude))).^4 + ...
       dd(5)*(cos(d2r(sp36.latitude))).^5 + ...
       dd(6)*(cos(d2r(sp36.latitude))).^6);
   
     %     % %Use mean Solar Modulation Parameter (SPhiInf)
cp.depthvector=[0:10:maxdepth];
%the following assignments of pp are so the code doesn't crash.  We don't
%use the output values, so they are meaningless. If they are included this
%way, we don't have to recode the muon code for each nuclide. Muons are
%time-independent as coded.

flux=muonfluxsato(cp.depthvector,sf36.P,RcEst,pp.SPhiInf,pp,'yes');
%store the output fluxes that we need
cp.negflux=flux.R;
cp.totalflux=flux.phi;
%Also store the muons production rates from the code
cp.muon36(1,:)=cp.depthvector;

%store muon scaling factor
cp.SFmufast=flux.SFmufast;
cp.SFmuslow=flux.SFmuslow;

% Need to calculate the chemical compound factor
% Start by calculating the capture probability for target elements
% (numerator)

PcapturetargetK=cp.tni(1)* pp.table(9,12);
PcapturetargetCa=cp.tni(2)* pp.table(10,12);
%for Ti & Fe
%PcapturetargetTi=cp.tni(3)*consts.table(11,12);
%PcapturetargetFe=cp.tni(4)*consts.table(13,12);

%calculate the overall probability of capture (denominator)
Pcapturebulk=cp.ni.*pp.table(:,12);
%divide these two in the correct direction
cp.FcompoundK=PcapturetargetK/sum(Pcapturebulk);
cp.FcompoundCa=PcapturetargetCa/sum(Pcapturebulk);

% For Ti & Fe
%FcompoundTi=PcapturetargetTi/sum(Pcapturebulk);
%FcompoundFe=PcapturetargetFe/sum(Pcapturebulk);
% negative muon capture

P_negK = cp.negflux.*pp.k_negpartial36K.*cp.FcompoundK.*pp.fstar36K;
P_negCa = cp.negflux.*pp.k_negpartial36Ca.*cp.FcompoundCa.*pp.fstar36Ca;

%calculating fast muon contribution to chlorine (individually for Ca and K)
  z=cp.depthvector;
% Balco original - from Heisinger fast muon paper Sea Level
% Beta = 0.846 - 0.015 .* log((z./100)+1) + 0.003139 .* (log((z./100)+1).^2);

% % For Beacon heights
% aalpha = 0.75;
% Beta =  0.842344 - 0.0107398 log((z./100)+1) + 0.00240182 log((z./100)+1)^2
% 
% aalpha = 0.85;
% Beta =  0.888695 - 0.00716992 log((z./100)+1) + 0.00169676 log((z./100)+1)^2
% 
aalpha = 1.0;
Beta = 1.0;
% 
% alpha = 1.15;
% Beta =  1.18129 + 0.00903804 log((z./100)+1) - 0.00273586 lnlog((z./100)+1)^2
% 
% aalpha = 1.30;
% Beta =  1.45283+0.04615 lnlog((z./100)+1) - 0.0153481 lnlog((z./100)+1)^2 + 0.000672339 lnlog((z./100)+1)^3
% % 

Ebar = 7.6 + 321.7.*(1 - exp(-8.059e-6.*z)) + 50.7.*(1-exp(-5.05e-7.*z));

P_fastK = cp.totalflux.*Beta.*(Ebar.^aalpha).*pp.sigma036K.*cp.tni(1);
P_fastCa = cp.totalflux.*Beta.*(Ebar.^aalpha).*pp.sigma036Ca.*cp.tni(2);
cp.P_fasttotal=P_fastK+P_fastCa;

%The following lines are for adding in Fe & Ti when we're ready
% P-negTi = R.*consts.k_neg_Ti*FcompoundTi;
% P-negFe = R.*consts.k_neg_Fe*FcompoundFe;
%Put in these as place holders until we can actually include them
% P_negTi=0;
% P_negFe=0;
% ProdmuTi=P_negTi;
% ProdmuFe=P_negFe;

%sum parts of the muon production

%Prodmu(jkl)=ProdmuK+ProdmuCa+P_negTi+P_negFe;
cp.Prodmu=P_negK+P_negCa+cp.P_fasttotal;
cp.muon36(2,:)=P_negCa+P_fastCa;
cp.muon36(3,:)=P_negK+P_fastK;
cp.muon36(4,:)=P_negCa;
cp.muon36(5,:)=P_negK;
cp.muon36(6,:)=P_fastCa;
cp.muon36(7,:)=P_fastK;
%
%  New Formulation uses Greg's Heisinger table from above to find the 
%  surface production of neutrons for that particular elevation.  
%first, interpolate and make sure we have the surface values

cp.negflux0=interpolate(cp.depthvector,cp.negflux,0);
cp.totalflux0=interpolate(cp.depthvector,cp.totalflux,0);

%calculate neutron production from muons
cp.Pmun0=pp.Ys*(cp.negflux0)+0.0000058*(cp.totalflux0);
%
% NCa and NK are useful temporary variables.
%
NCa=cp.tni(2);
NK=cp.tni(1);
CU=sp36.ci(18);
CTh=sp36.ci(19);
NCl=cp.tni(5);
%
% Pnsf
%
cp.Pnsf=0.429*CU;
%
% X
%
% Depends on fi (cp.fi), Si (table(:,6), and Yn_iu (table(:,7)
%
cp.X=sum(pp.table(:,6).*cp.fi.*pp.table(:,7))/sum(pp.table(:,6).*cp.fi);
%
% Y
%
cp.Y=sum(pp.table(:,6).*cp.fi.*pp.table(:,8))/sum(pp.table(:,6).*cp.fi);
%
% Pnan
%
cp.Pnan=cp.X*CU+cp.Y*CTh;
%
% Pethr
%
cp.Pethr=(cp.Pnan+cp.Pnsf)*(1-cp.pEthss);
%
% Pthr
%
cp.Pthr=(cp.Pnan+cp.Pnsf)*cp.pEthss;
%
% N36r         Radiogenic production in atoms/gram.
%
cp.N36r=cp.Pthr*cp.fth/pp.lambda36Cl+cp.Pethr*cp.feth/pp.lambda36Cl;
%
% N36m         Measured Chlorine 36.
%
cp.N36m=sp36.concentration36;
%
% N36c         Cosmogenic production in atoms/gram.
%
cp.N36c=cp.N36m-cp.N36r-sp36.inheritance36;

