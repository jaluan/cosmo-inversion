function out = fit_P_exp(z,Pnm,Pfm,n,plotFlag)

% This function generates an exponential approximation to subsurface
% production rates over a finite depth range. 
% 
% Syntax: out = fit_Pmu_exp(z,Pnm,Pfm,n,plotFlag);
%
% input args are:
%   depth vector x
%   Pnm: negative muon production profile
%   Pfm: fast muon production profile
%   n number of exponentials used to fit (maximum is 3)
%   plotFlag 1 to plot fitting results; 0 to not (default 1)
%
%   out has out.P surface production rates (atoms/g/yr)
%           and out.L e-folding lengths (g/cm2)
%   so then P_mu(z) = sum(out.P.*exp(-z./out.L));
%   out also has out.actP and out.predP, which are actual (i.e., calculated
%   with model 1A) and predicted production rates, out.z (depths of above), 
%   and out.scatter, which is std(out.predP./out.actP). 
%
% Written by Greg Balco
% Berkeley Geochronology Center
% October, 2016
% 
% Modified to fit input production profiles, JLA - April, 2019

% Checks

if nargin < 5; plotFlag = 1; end;

if n > 3; error('fit_Pmu_exp: n can''t be greater than 3');end;

% Define z vector
fitx = z;
y = Pnm+Pfm;

% Redundant assignments to output arg
out.z = fitx;
out.actP = y;

% Do fitting

% Always do 1-order fit to get starting point. 
fity = log(y);
pf1 = polyfit(fitx,fity,1);
P1 = exp(pf1(2)); 
L1 = -1./pf1(1);

if n == 1;
    fity = log(y);
    pf1 = polyfit(fitx,fity,1);
    out.P = P1; 
    out.L = L1;
    out.predP = exp(polyval(pf1,out.z));
    ts = 'Single exponential';
elseif n == 2;
    % Starting guess
    x0 = [P1/2 P1/2 L1*1.5 L1./1.5];
    xopt = fminsearch(@(x) sum(((x(1).*exp(-fitx./x(3)) + x(2).*exp(-fitx./x(4)))-y).^2),x0);
    out.P = xopt([1 2]);
    out.L = xopt([3 4]);
    out.predP = out.P(1).*exp(-fitx./out.L(1)) + out.P(2).*exp(-fitx./out.L(2));
    ts = 'Two exponentials';
elseif n == 3;
    % Starting guess
    x0 = [P1/4 P1/2 P1/4 L1*2 L1 L1./2];
    %xopt = fmincon(@(x) sum(((x(1).*exp(-fitx./x(4)) + x(2).*exp(-fitx./x(5)) + x(3).*exp(-fitx./x(6)))-y).^2),x0,[],[],[],[],[0 0 0 0 0 0],[P1*2 P1*2 P1*2 Inf Inf Inf]);
    xopt = fminsearch(@(x) sum(((x(1).*exp(-fitx./x(4)) + x(2).*exp(-fitx./x(5)) + x(3).*exp(-fitx./x(6)))-y).^2),x0);
    out.P = xopt([1 2 3]);
    out.L = xopt([4 5 6]);
    out.predP = out.P(1).*exp(-fitx./out.L(1)) + out.P(2).*exp(-fitx./out.L(2)) + out.P(3).*exp(-fitx./out.L(3));
    ts = 'Three exponentials';
end;

% Do plotting

if plotFlag == 1;
    figure;
    plot(out.actP,out.z,'go','markerfacecolor','g');
    hold on;
    plot(out.predP,out.z,'k');
    title(ts);
    xlabel('Pmu (atoms/g/yr)');
    set(gca,'ydir','reverse','xscale','log');
    ylabel('Depth (g/cm2)');
end;
    

out.scatter = std(out.predP./out.actP);
