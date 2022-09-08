% Calculate 36Cl production rates based on CronusCalc functions

% input: 
% z: depths of interest [g/cm^2]
% scaling_model: two-letter code, e.g. 'st','sa','de','du','li'
% sampledata: bulk rock and target chemistry of sample, structure as follows:
%
%1.     Sample cosmogenic 36-Cl concentration (atoms of 36-Cl/g of target)
%2.     Inheritance (atoms 36-Cl/g of target)  
%3.     erosion-rate epsilon (g/(cm^2*kyr))
%4.     fractional volumetric water-content (unitless) 
%5.     bulk density (g/cm^3)
%6.     sample thickness (cm)
%7.     Latitude (decimal degrees, -90(S) to +90(N))
%8.     Longitude (decimal degrees, 0-360 degrees east)
%9.     Elevation (meters)
%10.    Pressure (hPa)                Both 9 and 10 must be present!
%11.    Shielding factor for terrain, snow, etc. (unitless)
%12.    Effective attenuation length -Lambdafe (g/cm^2)
%13.    % CO2                        Rock
%14.    % Na2O                       Rock
%15.    % MgO                        Rock
%16.    % Al2O3                      Rock
%17.    % SiO2                       Rock
%18.    % P2O5                       Rock
%19.    % K2O                        Rock
%20.    % CaO                        Rock
%21.    % TiO2                       Rock
%22.    % MnO                        Rock
%23.    % Fe2O3                      Rock
%24.    Cl (ppm)                     Rock
%25.    B (ppm)                      Rock
%26.    Sm (ppm)                     Rock
%27.    Gd (ppm)                     Rock
%28.    U (ppm)                      Rock
%29.    Th (ppm)                     Rock
%30.    Cr (ppm)                     Rock
%31.    Li (ppm)                     Rock
%32.	Target element %K2O          Target
%33.    Target element %CaO          Target
%34.    Target element %TiO2         Target
%35.    Target element %Fe2O3        Target
%36.    Target element Cl (ppm)      Target
%37.    Depth to top of sample (g/cm^2)
%38.    Year sampled (e.g. 2010)
% 
% output: 
% out: structure containing averaged production for each production pathway
% within sample thickness (out.P36s, out.P36m, out.P36th, out.P36eth)

% Jane LA, April 2019

function out=get36clProd(sampledata,z,scaling_model)

% Make sampledata and uncertainties column vectors if they aren't already.
  
if (size(sampledata,1)==1)
  sampledata=sampledata';
end
  
%Check input data length
if (length(sampledata) ~= 38)
    error('sampledata has wrong size!');
end

% Setup the physical parameters
pp=physpars(scaling_model);

% Extract the sample parameters from the sampledatavector
sp=samppars36(sampledata);
  
% Get the scale factors
sf=scalefacs36(sp,scaling_model);

% We need an absolute maximum age for several purposes, including
% detecting saturated samples and setting the maximum depth for comppars.
maxage=2000;               % 2Ma > 6 half lives              
  
% Figure out the maximum possible depth at which we'll ever need a
% production rate.  This is depthtotop + maxage * erosion +
% thickness * density + a safety factor.
maxdepth=sp.depthtotop+maxage*sp.epsilon+sp.ls*sp.rb+1000;

%check for maximum depth: muon formulation is only good down to 2e5 g/cm2.
% This will likely never happen, but if it does, Matlab will error out saying we have not
% supplied times and plotprod and we will receive that error and be able to help the user.
if maxdepth > 2e5
  fprintf(1,'Maximum sample depth (%f) exceeds muon formulation of 2e5 g/cm2. \n Options: Lower the erosion rate, lower the maxage in file cl36age, or change muon formulation in muonfluxsato.m',[maxdepth])
  warning('This sample exceeds muon maximum depth. Try lowering the erosion rate of the sample or decreasing sample depth.');
  return;
end
  
% Computed parameters.
cp=comppars36(pp,sp,sf,maxdepth);
  
% Get contemporary depth production rates in atoms/g 
sf.currentsf=getcurrentsf(sf,0,scaling_model,'cl');
  
% % interpolate production across the thickness of the sample (n points)
% n=10;
% z=linspace(thickness/n/2,(n-1)*thickness/n+thickness/n/2,n);
% z=0;

% Compute production at each depth
out36=prodz36(z,pp,sf,cp);

% average results for each production pathway
out.P36s=out36.Prods;
out.P36nm=out36.Prodnm;
out.P36fm=out36.Prodfm;
out.P36th=out36.Prodth;
out.P36th1=out36.Prodth1;
out.P36th2=out36.Prodth2;
out.P36th3=out36.Prodth3;
out.P36th4=out36.Prodth4;
out.Lth1=cp.Lambdafe;
out.Lth2=cp.Lethss;
out.Lth3=cp.Lthss;
out.P36eth=out36.Prodeth;
out.P36eth1=out36.Prodeth1;
out.P36eth2=out36.Prodeth2;
out.P36eth3=out36.Prodeth3;
