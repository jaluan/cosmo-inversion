function compile_results_mn5(ss)

% This function compiles all results for a study area and saves compiled 
% results in the model output folder

% Inputs: 
% ss: season, 'MDML1': Heimefrontfjella, 'MDML2': Jutulstraumen.

% Written by David Lundbek Egholm, Aarhus University
% Modified by Jane Lund Andersen to also i) include [36Cl, 21Ne],
% ii) handle alternate climate inputs, and iii) add an (optional) slope or 
% step to the glaciation threshold (emergence)

close all;
addpath('Functions','ClimateCurves')

set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

set(gcf,'units','centimeters','position',[10,5,20,30]);
set(gcf,'paperunits','centimeters','paperposition',[1,5,20,30]);
set(gca,'position',[0,0,1,1],'visible','off');
set(gca,'xlim',[0,1],'ylim',[0,1]);

pcol1 = [.7,.7,.9];
prt = [2.5,25,50,75,97.5]; %percentile values

%% get relationship between dO and burial in past 1Myr
load('burialTimeZachos.mat','bt','dO_i')

%% get data
load(['data/' ss '_mn_data.mat'], 'sample');
Ns = length(sample); %Option to specify subset of samples only
List = 1:Ns; %Option to specify sample #

%% Output parameters
%times for measuring erosion
TE1 = 1.0; %erosion since TE1/depth at TE1;
TE2 = 2.6;
TE3 = 5.3;
TE4 = 15.0;

%depth for recording times
ZT1 = 1.0; %Time spent within ZT1 m of surface
ZT2 = 2.6;
ZT3 = 0.25;

di = 0.25; %distance for plotting
maxres = 100; % maximum residual

%% ******* initialise figures *************
ax1=subplot(8,1,1);
hold on; box on;
set(gca,'ylim',[3.0,5.0]);
set(gca,'xlim',[0,Ns+1]);
title('d18O threshold');
ylabel('d18O (permil)');
grid on

ax2=subplot(8,1,2);
hold on; box on;
set(gca,'ylim',[0,10]);
set(gca,'xlim',[0,Ns+1]);
title(['Erosion since ' num2str(TE1) ' myr']);
ylabel('E1 (m)');
grid on

ax3=subplot(8,1,3);
hold on; box on;
set(gca,'ylim',[0,10]);
set(gca,'xlim',[0,Ns+1]);
title(['Erosion since ' num2str(TE2) ' myr']);
ylabel('E2 (m)');
grid on

ax4=subplot(8,1,4);
hold on; box on;
set(gca,'ylim',[0,20]);
set(gca,'xlim',[0,Ns+1]);
title(['Erosion since ' num2str(TE3) ' myr']);
ylabel('E3 (m)');
grid on

ax5=subplot(8,1,5);
hold on; box on;
set(gca,'ylim',[0,5.0]);
set(gca,'xlim',[0,Ns+1]);
title(['Erosion since ' num2str(TE4) ' myr']);
ylabel('E4 (m)');
grid on

ax6=subplot(8,1,6);
hold on; box on;
set(gca,'ylim',[0,10.0]);
set(gca,'xlim',[0,Ns+1]);
title(['Time of ' num2str(ZT1) ' m erosion']);
ylabel('T1 (Ma)');
grid on

ax7=subplot(8,1,7);
hold on; box on;
set(gca,'ylim',[0,10.0]);
set(gca,'xlim',[0,Ns+1]);
title(['Time of ' num2str(ZT2) ' m erosion']);
ylabel('T2 (Ma)');
grid on

ax8=subplot(8,1,8);
hold on; box on;
set(gca,'ylim',[0,10.0]);
set(gca,'xlim',[0,Ns+1]);
title(['Time of ' num2str(ZT3) ' m erosion']);
ylabel('T3 (Ma)');
grid on

%% loop bedrock samples
for i=1:Ns
    
    %sample number
    sn = List(i);
    
    %load MC results
    mname = ['models/MDML/MC_sample_',num2str(sn),'_mn5_zachos_',ss,'.mat'];
    load(mname,'model');

    axes(ax1);
    %extract total erosion data
    d18O = [];
    res = [];
    E1 = [];
    E2 = [];
    E3 = [];
    E4 = [];
    TZ1 = [];
    TZ2 = [];
    TZ3 = [];
    b1My = [];

    n0 = 2;
    
    for nw = 1:model.Nwalk %loop over walkers
        I = find((model.walker{nw}.status > 0)&(model.walker{nw}.restot < maxres)); %accepted walkers with residual < maxres
        if (~isempty(I))
            d18Ow = model.walker{nw}.u(1,I);
            resw = model.walker{nw}.restot(I); 
                        
            d18O = [d18O(:);d18Ow(:)];
            res = [res(:);resw(:)];
            
            E1w = [];
            E2w = [];
            E3w = [];
            E4w = [];
            TZ1w = [];
            TZ2w = [];
            TZ3w = [];
            b1Myw = [];
            
            for k=1:length(I)
                        
                T1 = model.walker{nw}.u(2,I(k));
                z1 = model.walker{nw}.u(n0+1,I(k));
                dT2 = model.walker{nw}.u(n0+2,I(k));
                dz2 = 10^model.walker{nw}.u(n0+3,I(k));
                dT3 = model.walker{nw}.u(n0+4,I(k));
                dz3 = 10^model.walker{nw}.u(n0+5,I(k));
                dT4 = model.walker{nw}.u(n0+6,I(k));
                dz4 = 10^model.walker{nw}.u(n0+7,I(k));
                E5 = 10^model.walker{nw}.u(n0+8,I(k));
            
                T2 = T1 + dT2;
                z2 = z1 + dz2;
                T3 = T2 + dT3;
                z3 = z2 + dz3;
                T4 = T3 + dT4;
                z4 = z3 + dz4;
                T5 = model.age; %this requires age > Tdg+dT2+dT3
                z5 = z4 + (model.age - T4)*E5;
                
                Tm = [0,T1,T2,T3,T4,T5]; %Ma
                zm = [0,z1,z2,z3,z4,z5]; %m

                E1w(k) = interp1(Tm,zm,TE1); %interpolate time-depth vector at time TE1
                E2w(k) = interp1(Tm,zm,TE2); %interpolate time-depth vector at time TE2
                E3w(k) = interp1(Tm,zm,TE3); %interpolate time-depth vector at time TE3
                E4w(k) = z5; %interp1(Tm,zm,TE4); %Depth at beginning of model run
                
                if (z5 > ZT1)
                    TZ1w(k) = interp1(zm,Tm,ZT1); %interpolate depth-time vectors at depth ZT1
                else
                    TZ1w(k) = (ZT1 - z5)/E5 + T5; %Calculate time if ZT1 is below max sample depth
                end
                if (z5 > ZT2)
                    TZ2w(k) = interp1(zm,Tm,ZT2); %interpolate depth-time vectors at depth ZT2
                else
                    TZ2w(k) = (ZT2 - z5)/E5 + T5; %Calculate time if ZT2 is below max sample depth
                end
                if (z5 > ZT3)
                    TZ3w(k) = interp1(zm,Tm,ZT3);  %interpolate depth-time vectors at depth ZT3
                else
                    TZ3w(k) = (ZT3 - z5)/E5 + T5; %Calculate time if ZT3 is below max sample depth
                end
                
                b1Myw(k) = interp1(dO_i,bt,d18Ow(k)); %Find ice burial since 1 Ma (%) from glaciation threshold value
                
            end
            
            %save things
            E1 = [E1(:);E1w(:)];
            E2 = [E2(:);E2w(:)];
            E3 = [E3(:);E3w(:)];
            E4 = [E4(:);E4w(:)];
            TZ1 = [TZ1(:);TZ1w(:)];
            TZ2 = [TZ2(:);TZ2w(:)];
            TZ3 = [TZ3(:);TZ3w(:)];
            b1My = [b1My(:);b1Myw(:)];
            
        end
    end
    
    %Calculate percentiles
    d18Oprc = prctile(d18O,prt);
    E1prc = prctile(E1,prt);
    E2prc = prctile(E2,prt);
    E3prc = prctile(E3,prt);
    E4prc = prctile(E4,prt);
    TZ1prc = prctile(TZ1,prt);
    TZ2prc = prctile(TZ2,prt);
    TZ3prc = prctile(TZ3,prt);
    b1Myprc = prctile(b1My,prt);
    
    %Box-and-whisker plots
    axes(ax1); myboxwhisker(i,di,d18Oprc,pcol1);
    axes(ax2); myboxwhisker(i,di,E1prc,pcol1);
    axes(ax3); myboxwhisker(i,di,E2prc,pcol1);
    axes(ax4); myboxwhisker(i,di,E3prc,pcol1);
    axes(ax5); myboxwhisker(i,di,E4prc,pcol1);
    axes(ax6); myboxwhisker(i,di,TZ1prc,pcol1);
    axes(ax7); myboxwhisker(i,di,TZ2prc,pcol1);
    axes(ax8); myboxwhisker(i,di,TZ3prc,pcol1);
    
    %save things
    sample{List(i)}.d18O = d18O; 
    sample{List(i)}.d18Oprc = d18Oprc;
    sample{List(i)}.TE1 = TE1;  
    sample{List(i)}.E1 = E1; 
    sample{List(i)}.E1prc = E1prc; 
    sample{List(i)}.TE2 = TE2;  
    sample{List(i)}.E2 = E2; 
    sample{List(i)}.E2prc = E2prc; 
    sample{List(i)}.TE3 = TE3;  
    sample{List(i)}.E3 = E3;
    sample{List(i)}.E3prc = E3prc;
    sample{List(i)}.TE4 = TE4;  
    sample{List(i)}.E4 = E4; 
    sample{List(i)}.E4prc = E4prc;
    sample{List(i)}.ZT1 = ZT1;   
    sample{List(i)}.TZ1 = TZ1; 
    sample{List(i)}.TZ1prc = TZ1prc; 
    sample{List(i)}.ZT2 = ZT2;   
    sample{List(i)}.TZ2 = TZ2; 
    sample{List(i)}.TZ2prc = TZ2prc; 
    sample{List(i)}.ZT3 = ZT3;   
    sample{List(i)}.TZ3 = TZ3; 
    sample{List(i)}.TZ3prc = TZ3prc; 
    sample{List(i)}.b1My = b1My; 
    sample{List(i)}.b1Myprc = b1Myprc;
    sample{List(i)}.residual = res;  

end

%% save to file
save(['models/MDML/compiled_results_' ss '_zachos_acc.mat'],'sample','-v7.3'); %-v7.3 for variables larger tham 2GB

%% return to default interpreters
set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default'); 
set(groot, 'defaultLegendInterpreter','default');