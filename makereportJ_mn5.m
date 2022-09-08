function makereportJ_mn5(snr,CCurve,emergence,ss)

% function makereportJ_mn5(snr,CCurve,emergence,ss)

% This function plots: 
% a. Measured vs. modelled Be-Al, 
% b. histogram of d18O-threshold parameter and c. last deglaciation pm,
% d. 2D-histogram of exhumation pathways of sample to surface through time,
% (e-h). If the model contains 36Cl or 21Ne, the measured vs modelled
% concentrations compared to 10Be and 26Al will be shown,
% i. Distribution of all model parameters with separate coloured curves for 
% each walker.
% j. parameter values vs model-time through burn-in (gray) and post-burn in
% (colours, rejected models are grey). Separate colours for different walkers.

% Inputs: 
%  - snr: scalar or vector containing sample numbers of samples in input data
%       file generated from compile_mn_data.m.
%  - CCurve: Climate curve, 1: d18Ocurves.mat, 2:zachos.mat,
%       3: milankovitch.mat, 4-6:iceDMC.mat
%  - emergence: Threshold to climate curve exposure parameter/'d18O threshold'. 
%       0: no slope to threshold, 1=slope to threshold (only tested with d18O 
%       curve), 2=two-step threshold changing at time dTTh (hardcoded for now)
%  - ss: season, specific to MDML samples, only for data loading, can be
%       removed.

% Note that the inversion must have been generated first by use of bedrockMC-code

% Written by David Lundbek Egholm, Aarhus University
% Modified by Jane Lund Andersen to also i) include [36Cl, 21Ne],
% ii) handle alternate climate inputs, and iii) add an (optional) slope or 
% step to the glaciation threshold (emergence)

close all;
set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

addpath('Functions','Functions/export_fig')
CNprop = getCNprop(); %Cosmogenic properties

%% load MC results
str = num2str(snr(1));
if length(snr) > 1 %if multiple samples were modelled together
    for i=2:length(snr)
        str = [str,['-',num2str(snr(i))]];
    end
end

% define model name depending on input
if ss == 1
    season='MDML1';
elseif ss==2
    season='MDML2';
end

if CCurve == 1
    if emergence == 0
        mname = ['models/MDML/MC_sample_',str,'_mn5_d18Ocurves_',season,'.mat'];
    elseif emergence == 1
        mname = ['models/MDML/MC_sample_',str,'_mn5_d18Ocurves_',season,'_emergence.mat'];
    elseif emergence == 2
        mname = ['models/MDML/MC_sample_',str,'_mn5_d18Ocurves_',season,'_twostep.mat'];
    end
elseif CCurve == 2
    if emergence == 0
        mname = ['models/MDML/MC_sample_',str,'_mn5_zachos_',season,'.mat'];
    elseif emergence == 1
        mname = ['models/MDML/MC_sample_',str,'_mn5_zachos_',season,'_emergence.mat'];
    end
elseif CCurve == 3
    mname = ['models/MDML/MC_sample_',str,'_mn5_milankovitch_',season,'.mat'];      
elseif CCurve == 4
    mname = ['models/MDML/MC_sample_',str,'_mn5_IceDMC1_',season,'.mat'];
elseif CCurve == 5
    mname = ['models/MDML/MC_sample_',str,'_mn5_IceDMC2_',season,'.mat'];
elseif CCurve == 6
    mname = ['models/MDML/MC_sample_',str,'_mn5_IceDMC3_',season,'.mat'];
else
    warning('ClimateCurve not implemented')    
end

load(mname,'model');

%% colors
col1 = 0.9*[1,1,1];
col2 = 0.6*[1,1,1];
col3 = 0.3*[1,1,1];

map = colormap;
nc = length(map);
xc = linspace(0,1,nc);
rc = map(:,1);
gc = map(:,2);
bc = map(:,3);
wcol=zeros(model.Nwalk,3);
for i=1:model.Nwalk
    ri = interp1(xc,rc,(i-1)/model.Nwalk);
    gi = interp1(xc,gc,(i-1)/model.Nwalk);
    bi = interp1(xc,bc,(i-1)/model.Nwalk);
    wcol(i,:) = [ri,gi,bi];
end

%% figures showing walkers *************

set(gcf,'units','normalized','position',[.1,.3,.4,.6]);
set(gca,'position',[0,0,1,1],'visible','off');
set(gca,'xlim',[0,1],'ylim',[0,1]);
set(gcf,'Name','Walker information');
% ax0 = gca;

[np,~] = numSubplots(model.Nmp);

%plot models parameters as a function of model time
for i = 1:model.Nmp
    subplot(np(1),np(2),i); hold on; box on; grid on;
    xlabel('model nr.');
    ylabel(model.mp{i}.name);
    set(gca,'xlim',[0,length(model.walker{1}.status)],'ylim',[model.mp{i}.vmin,model.mp{i}.vmax]);
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == -2); %rejected during burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',col1);
        I = find(model.walker{nw}.status == 0); %rejected after burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',col2);
        I = find(model.walker{nw}.status == -1); %accepted during burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',col3);
        I = find(model.walker{nw}.status == 1); %accepted after burnin
        plot(I,model.walker{nw}.up(i,I),'.','color',wcol(nw,:));
    end
end

print('temp3.pdf','-dpdf','-fillpage'); %print temporary pdf, appended below

%% Parameter distributions ************
figure()
set(gcf,'units','normalized','position',[.3,.2,.4,.6]);
set(gcf,'Name','Parameter distributions');

[np,~] = numSubplots(model.Nmp);

%loop parameters
for i=1:model.Nmp
    
    subplot(np(1),np(2),i); 
    hold on; box on; grid on;
    ylabel('Frequency');
    xlabel(model.mp{i}.name);
    set(gca,'xlim',[model.mp{i}.vmin, ...
                    model.mp{i}.vmax]);
    uval = [];
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        if (~isempty(I))
            [f,xi] = ksdensity(model.walker{nw}.up(i,I));
            line(xi,f,'color',wcol(nw,:));
            uval = [uval(:);model.walker{nw}.up(i,I)'];
        end
    end
    if (~isempty(uval))
        [f,xi] = ksdensity(uval); %compute kernel density
    end
    line(xi,f,'color','k','linewidth',2);

end

print('temp2.pdf','-dpdf','-fillpage'); %print temporary pdf, appended below

%% set plot margins
lpm = 0.1;
ddm = 0.02;

%% 36Cl and 21Ne model-data fit figures 
%36Cl figures
if (any(ismember(model.data{1,1}.nuclides,3))) %if 36Cl in model
    % initiate figure
    figure()
    set(gcf,'units','normalized','position',[.3,.2,.4,.6]);
    set(gcf,'Name','36Cl model-data fit');
    ax0 = gca;
    set(ax0,'position',[0,0,1,1]);
    set(gca,'visible','off');
    dpx = (1-2*lpm-(model.Nsnr-1)*ddm)/model.Nsnr;
    
    for ns = 1:model.Nsnr %loop samples

        axes('position',[lpm+(ns-1)*(ddm+dpx),0.7,dpx,0.25]);
        hold on; box on; grid on;
        title(model.data{ns}.name);

        n0 = (ns-1)*model.Nds;
        
        % retrieve modelled nuclide distributions
        Be = [];
        Al = [];
        Cl = [];
        for nw = 1:model.Nwalk
            I = find(model.walker{nw}.status == 1);
            Be = [Be(:);model.walker{nw}.gm(n0+1,I)'];
            Al = [Al(:);model.walker{nw}.gm(n0+2,I)'];
            Cl = [Cl(:);model.walker{nw}.gm(n0+3,I)'];
        end
        %calculate ratios
        ClBe = Cl./Be;
        ClAl = Cl./Al;
        
        Beint = linspace(min(Be),max(Be),40)';
        ClBeint = linspace(0.5*min(ClBe),1.5*max(ClBe),40)';
        Alint = linspace(min(Al),max(Al),40)';
        ClAlint = linspace(0.5*min(ClAl),1.5*max(ClAl),40)';

        %Cl-Be banana
        xlabel(['Normalized NBe (atoms/g), sample ',num2str(ns)]);
        if (ns == 1), ylabel('N36Cl/N10Be');
        else, set(gca,'yticklabel',[]);
        end

        N = hist3([ClBe,Be],{ClBeint Beint});
        N = N/sum(N(:));
        
        %Normalize to SLHL using surface production rate (spallation+muons)
        nfac = CNprop.PBe0/(model.data{ns}.production.P10spal +...
            model.data{ns}.production.P10_m1 + model.data{ns}.production.P10_m2);

        [X,Y] = meshgrid(Beint*nfac,ClBeint);
        contour(X,Y,N,40);
    
        errorbar(model.data{ns}.N10*nfac,model.data{ns}.r3610,model.data{ns}.dN10*nfac,'horizontal','.k');
        errorbar(model.data{ns}.N10*nfac,model.data{ns}.r3610,model.data{ns}.dr3610,'vertical','.k');
        
        %Cl-Al banana
        axes('position',[lpm+(ns-1)*(ddm+dpx),0.075,dpx,0.25]);
        hold on; box on; grid on;
        
        xlabel(['Normalized NAl (atoms/g), sample ',num2str(ns)]);
        if (ns == 1), ylabel('N36Cl/N26Al');
        else, set(gca,'yticklabel',[]);
        end

        N = hist3([ClAl,Al],{ClAlint Alint});
        N = N/sum(N(:));
        
        %Normalize to SLHL using surface production rate (spallation+muons)
        nfac = CNprop.PAl0/(model.data{ns}.production.P26spal + ...
            model.data{ns}.production.P26_m1 + model.data{ns}.production.P26_m2);
        
        [X,Y] = meshgrid(Alint*nfac,ClAlint);
        contour(X,Y,N,40);
        
        model.data{ns}.r3626 = model.data{ns}.N36./model.data{ns}.N26;
        model.data{ns}.dr3626 = model.data{ns}.r3626.*sqrt((model.data{ns}.dN26./...
        model.data{ns}.N26).^2+(model.data{ns}.dN36./model.data{ns}.N36).^2);
    
        errorbar(model.data{ns}.N26*nfac,model.data{ns}.r3626,model.data{ns}.dN26*nfac,'horizontal','.k');
        errorbar(model.data{ns}.N26*nfac,model.data{ns}.r3626,model.data{ns}.dr3626,'vertical','.k');
       
        % Histogram/kernel density with modelled concentrations
        axes('position',[lpm+(ns-1)*(ddm+dpx),0.4,dpx,0.2]);
        hold on; box on; grid on;
        xlabel('N36 (atoms/g)'); ylabel ('Probability')
        histogram(Cl,'Normalization','probability')

        % Measured data on top
        errorbar(model.data{ns}.N36,.015,model.data{ns}.dN36,...
        'horizontal','.k','Linewidth',1.5,'Color',[.7 .2 .2]);
    
    end
    print('temp1.pdf','-dpdf','-fillpage'); %print temporary pdf, appended below
end

%21Ne figures
if (any(ismember(model.data{1,1}.nuclides,4))) %if 21Ne in model
    % initiate figure
    figure()
    set(gcf,'units','normalized','position',[.3,.2,.4,.6]);
    set(gcf,'Name','21Ne model-data fit');
    ax0 = gca;
    set(ax0,'position',[0,0,1,1]);
    set(gca,'visible','off');

    dpx = (1-2*lpm-(model.Nsnr-1)*ddm)/model.Nsnr;
    
    for ns = 1:model.Nsnr %loop samples

        axes('position',[lpm+(ns-1)*(ddm+dpx),0.7,dpx,0.25]);
        hold on; box on; grid on;
        title(model.data{ns}.name);

        n0 = (ns-1)*model.Nds;

        % retrieve modelled nuclide distributions
        Be = [];
        Ne = [];
        Al = [];
        for nw = 1:model.Nwalk
            I = find(model.walker{nw}.status == 1);
            Be = [Be(:);model.walker{nw}.gm(n0+1,I)'];
            Al = [Al(:);model.walker{nw}.gm(n0+2,I)'];
            Ne = [Ne(:);model.walker{nw}.gm(n0+model.Nds,I)'];
        end
        BeNe = Be./Ne; %calculate ratios
        AlNe = Al./Ne;

        Neint = linspace(min(Ne),max(Ne),40)';
        BeNeint = linspace(0.5*min(BeNe),1.5*max(BeNe),40)';
        AlNeint = linspace(0.5*min(AlNe),1.5*max(AlNe),40)';

        %Ne-Be banana
        xlabel(['Normalized NNe (atoms/g), sample ',num2str(ns)]);
        if (ns == 1), ylabel('N10Be/N21Ne');
        else, set(gca,'yticklabel',[]);
        end

        N = hist3([BeNe,Ne],{BeNeint Neint});
        N = N/sum(N(:));
        
        %Normalize to SLHL using surface production rate (spallation+muons)
        nfac = CNprop.PNe0/(model.data{ns}.production.P21spal + ...
            model.data{ns}.production.P21_m1 + model.data{ns}.production.P21_m2);
        [X,Y] = meshgrid(Neint*nfac,BeNeint);
        contour(X,Y,N,40);

        errorbar(model.data{ns}.N21*nfac,model.data{ns}.r1021,model.data{ns}.dN21*nfac,'horizontal','.k');
        errorbar(model.data{ns}.N21*nfac,model.data{ns}.r1021,model.data{ns}.dr1021,'vertical','.k');
        
        %Ne-Al banana
        axes('position',[lpm+(ns-1)*(ddm+dpx),0.075,dpx,0.25]);
        hold on; box on; grid on;
        
        xlabel(['Normalized NNe (atoms/g), sample ',num2str(ns)]);
        if (ns == 1), ylabel('N26Al/N21Ne');
        else, set(gca,'yticklabel',[]);
        end

        N = hist3([AlNe,Ne],{AlNeint Neint});
        N = N/sum(N(:));
        
        %Normalize to SLHL using surface production rate (spallation+muons)
        nfac = CNprop.PNe0/(model.data{ns}.production.P21spal + ...
            model.data{ns}.production.P21_m1 + model.data{ns}.production.P21_m2);
        
        [X,Y] = meshgrid(Neint*nfac,AlNeint);
        contour(X,Y,N,40);
        
        model.data{ns}.r2621 = model.data{ns}.N26./model.data{ns}.N21;
        model.data{ns}.dr2621 = model.data{ns}.r2621.*sqrt((model.data{ns}.dN26./...
        model.data{ns}.N26).^2+(model.data{ns}.dN21./model.data{ns}.N21).^2);
    
        errorbar(model.data{ns}.N21*nfac,model.data{ns}.r2621,model.data{ns}.dN21*nfac,'horizontal','.k');
        errorbar(model.data{ns}.N21*nfac,model.data{ns}.r2621,model.data{ns}.dr2621,'vertical','.k');
        
        % Histogram/kernel density with modelled concentrations
        axes('position',[lpm+(ns-1)*(ddm+dpx),0.4,dpx,0.2]);
        hold on; box on; grid on;
        xlabel('N21 (atoms/g)'); ylabel ('Probability')
        histogram(Ne,'Normalization','probability')

        % Measured data on top
        errorbar(model.data{ns}.N21,.015,model.data{ns}.dN21,...
        'horizontal','.k','Linewidth',1.5,'Color',[.7 .2 .2]);
    
    end
    print('temp0.pdf','-dpdf','-fillpage'); %print temporary pdf, appended below
end

%% ******* main report figure *************
% Initiate figure
figure;
set(gcf,'papertype','a4');
set(gcf,'units','centimeters','position',[5,5,21,29.7]);
set(gcf,'Name','Report');
ax0 = gca;
set(ax0,'position',[0,0,1,1]);
set(gca,'visible','off');
text(0.05,0.97,['File: ',mname],'HorizontalAlignment','left','fontsize',14);

dpx = (1-2*lpm-(model.Nsnr-1)*ddm)/model.Nsnr;

% Al-Be banana
for ns = 1:model.Nsnr %loop samples

    axes('position',[lpm+(ns-1)*(ddm+dpx),0.675,dpx,0.25]);
    hold on; box on; grid on;
    set(gca,'ylim',[3,7.5]);
    
    %title
    Sname=['$',model.data{1,ns}.name{1},'$'];
    title(Sname);
    
    n0 = (ns-1)*model.Nds(ns);
    
    Be = [];
    Al = [];
    
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        Be = [Be(:);model.walker{nw}.gm(n0+1,I)'];
        Al = [Al(:);model.walker{nw}.gm(n0+2,I)'];
    end
    AlBe = Al./Be;
    
    Beint = linspace(min(Be),max(Be),40)';
    %Alint = linspace(min(Al),max(Al),40)';
    AlBeint = linspace(0.5*min(AlBe),1.5*max(AlBe),40)';
    
    xlabel(['Normalized NBe (atoms/g), sample ',num2str(ns)]);
    if (ns == 1), ylabel('NAl/NBe');
    else, set(gca,'yticklabel',[]);
    end
    
    N = hist3([AlBe,Be],{AlBeint Beint});
    N = N/sum(N(:));

    nfac = CNprop.PBe0/(model.data{ns}.production.P10spal +...
        model.data{ns}.production.P10_m1 + model.data{ns}.production.P10_m2);
    %nfac2 = CNprop.PAl0/(model.data{ns}.production.P26spal +...
    %    model.data{ns}.production.P26_m1 + model.data{ns}.production.P26_m2);
    [X,Y] = meshgrid(Beint*nfac,AlBeint);
    contour(X,Y,N,40);
    
    errorbar(model.data{ns}.N10*nfac,model.data{ns}.r2610,model.data{ns}.dN10*nfac,'horizontal','.k');
    errorbar(model.data{ns}.N10*nfac,model.data{ns}.r2610,model.data{ns}.dr2610,'vertical','.k');
    
    % Histogram/kernel density with modelled N26 concentrations
    axes('position',[0.65,0.8,0.2,0.1]);
    hold on; box on; grid on;
    xlabel('N26 (atoms/g)'); %ylabel ('Probability')
    histogram(Al,'Normalization','probability')

    % Measured data on top
    errorbar(model.data{ns}.N26,.015,model.data{ns}.dN26,...
    'horizontal','.k','Linewidth',1.5,'Color',[.7 .2 .2]);
end

% glaciation history parameter distibutions

dpx = (1-2*lpm-ddm)/2;

%loop parameters
for i=1:2
    
    axes('position',[lpm+(i-1)*(ddm+dpx),0.37,dpx,0.25]);
    hold on; box on; grid on;
    
    if (i == 1), ylabel('Frequency'); end
    xlabel(model.mp{i}.name);
    set(gca,'xlim',[model.mp{i}.vmin, ...
                    model.mp{i}.vmax]);
    set(gca,'yticklabel',[]);
    uval = [];
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        if (~isempty(I))
            [f,xi] = ksdensity(model.walker{nw}.up(i,I));
            line(xi,f,'color',wcol(nw,:));
            uval = [uval(:);model.walker{nw}.up(i,I)'];
        end
    end
    if (~isempty(uval))
        [f,xi] = ksdensity(uval); %compute kernel density
    end
    line(xi,f,'color','k','linewidth',2);

end

% Draw 2D-histogram of exhumation
Maxz = model.z0;
Maxt = model.age;
z0 = model.z0;
Nz = 100;
Nt = 200;
zint = linspace(0,Maxz,Nz);
tint = linspace(0,Maxt,Nt);
[tbin,zbin]=meshgrid(tint,zint);

dpx = (1-2*lpm-(model.Nsnr-1)*ddm)/model.Nsnr;

map = colormap;
map(1,:) = [1,1,1];
colormap(map);

for ns = 1:model.Nsnr

    %subplot(np(1),np(2),ns); 
    axes('position',[lpm+(ns-1)*(ddm+dpx),0.05,dpx,0.25]);
    hold on; box on; set(gca,'layer','top');
    
    set(gca,'ydir','reverse');
    xlabel('Time (Ma)');
    if (ns == 1) 
        ylabel('Burial depth');
    else
        set(gca,'yticklabel',[]);
    end
    
    histgrid = zeros(size(tbin));
    dzi = zint(2)-zint(1);
    
    %title
    Sname=['$',model.data{1,ns}.name{1},'$'];
    title(Sname);
    
    n0 = model.Mmp + (ns-1)*model.Nsmp;
    
    for nw=1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        for i=1:length(I)
            
            T1 = model.walker{nw}.u(2,I(i));
            z1 = model.walker{nw}.u(n0+1,I(i));
            dT2 = model.walker{nw}.u(n0+2,I(i));
            dz2 = 10^model.walker{nw}.u(n0+3,I(i));
            dT3 = model.walker{nw}.u(n0+4,I(i));
            dz3 = 10^model.walker{nw}.u(n0+5,I(i));
            dT4 = model.walker{nw}.u(n0+6,I(i));
            dz4 = 10^model.walker{nw}.u(n0+7,I(i));
            E5 = 10^model.walker{nw}.u(n0+8,I(i));
            
            T2 = T1 + dT2;
            z2 = z1 + dz2;
            T3 = T2 + dT3;
            z3 = z2 + dz3;
            T4 = T3 + dT4;
            z4 = z3 + dz4;
            T5 = model.age; %this requires age > Tdg+dT2+dT3+dT4
            z5 = z4 + (model.age - T4)*E5;
            
            
            Tm = [0,T1,T2,T3,T4,T5];
            zm = [0,z1,z2,z3,z4,z5];
            
            zinterp = interp1(Tm,zm,tint,'linear',2*z0);        
            izs = 1+floor(zinterp/dzi); %bin index at each tsfine
            for itsfine=1:Nt
                if (izs(itsfine)<=Nz)
                    histgrid(izs(itsfine),itsfine)=histgrid(izs(itsfine),itsfine)+ 1;
                end
            end
        end
    end
    
    
    contourf(tbin,zbin,(histgrid+0.5).^0.25,50, ...
                     'linestyle','none');
    
    set(gca,'ylim',[0 3])            
end

%% Name of output file
iDMCn = ["Spector";"Pollard";"deBoer"];

if model.emergence==0
    if CCurve <=3 || CCurve >=7
        mname = ['models/MDML/reports/Report_sample_',str,'_',model.ClimCurve,'_',season,'_mn5.pdf'];
    else
        mname = ['models/MDML/reports/Report_sample_',str,'_',model.ClimCurve,char(iDMCn(CCurve-3)),'_',season,'_mn5.pdf'];
    end
elseif model.emergence==1
    mname = ['models/MDML/reports/Report_sample_',str,'_',model.ClimCurve,'_',season,'_emergence_mn5.pdf'];
elseif model.emergence==2
    mname = ['models/MDML/reports/Report_sample_',str,'_',model.ClimCurve,'_',season,'_twostep_mn5.pdf'];
end

%% save figures to one pdf and delete temporary files
print(mname,'-dpdf','-fillpage');

if (any(ismember(model.data{1,1}.nuclides,3))) %if 36Cl in model
    append_pdfs(mname, 'temp1.pdf')
    delete('temp1.pdf')
end
if (any(ismember(model.data{1,1}.nuclides,4))) %if 21Ne in model
    append_pdfs(mname, 'temp0.pdf')
    delete('temp0.pdf')
end

append_pdfs(mname, 'temp2.pdf', 'temp3.pdf')
delete('temp2.pdf','temp3.pdf')
%% Set interpreter to default
set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default'); 
set(groot, 'defaultLegendInterpreter','default');