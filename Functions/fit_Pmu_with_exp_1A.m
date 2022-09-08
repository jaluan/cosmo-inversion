function out = fit_Pmu_with_exp_1A(p,mindepth,maxdepth,n,mc,plotFlag)

% This function generates an exponential approximation to subsurface
% production rates due to muons over a finite depth range. 
% 
% Syntax: out = fit_Pmu_with_exp_1A(elv,mindepth,maxdepth,n,mc);
%
% Uses "model 1A" with alpha = 1
%
% input args are:
%   p site pressure (hPa)
%   mindepth, maxdepth range of depth needed (g/cm2)
%   n number of exponentials used to fit (maximum is 3)
%   mc muon cross-sections structure - example (for Be-10):
%       mc.Natoms = 2.006e22;
%       mc.k_neg = 0.00191 .* 0.704 .* 0.1828; 
%       mc.sigma0 = 0.280e-30; 
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

% Checks

if nargin < 6; plotFlag = 1; end;

if n > 3; error('fit_Pmu_with_exp_1A: n can''t be greater than 3');end;

% Define z vector
fitx = linspace(mindepth,maxdepth,100); % Not sure how many to use...
y = zeros(size(fitx));

% % Compute "actual" production rates using model 1A
% for a = 1:length(fitx);
%     y(a) = P_mu_total_alpha1(fitx(a),p,mc);
% end;
y = P_mu_total_alpha1(fitx,p,mc); %JLA, March 12th 2019 - ok?

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
