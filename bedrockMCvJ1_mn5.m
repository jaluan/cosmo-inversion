function bedrockMCvJ1_mn5(snr,CCurve,emergence,ss)

% function bedrockMCvJ1_mn5(snr,CCurve,emergence,ss)
% Inverse MCMC model exploring complex exposure and exhumation histories of 
% samples with measured cosmogenic nuclides from (formerly) glaciated
% landscapes. Delineates the ensemble of parameter values with lowest 
% data misfit.

% Inputs: 
% - snr: scalar or vector containing sample numbers of samples in input data
%    file generated from compile_mn_data.m. If multiple numbers are given, the
%    exposure history of the samples will be assumed to be the same, while
%    exhumation may vary. Note that the function retrieves deglaciation timing
%    from the first sample only.
% - CCurve: Climate curve input controlled by number: 
%    1: d18Ocurves.mat (Lisiecki & Raymo 2005), 
%      Data downloaded from https://lorraine-lisiecki.com/stack.html
%    2: zachos.mat (Zachos et al., 2001; default in manuscript),
%    3: milankovitch.mat. Data described in Laskar 2011 (La2004) and 
%      downloaded here: https://biocycle.atmos.colostate.edu/shiny/Milankovitch/
%    4-6:iceDMC.mat, Ice-sheet model output data generated for DML
%    localities using the ICE-D: DMC interface: http://dmc.ice-d.org by
%    Perry Spector. 4=Spector et al., 2018, 5=Pollard & DeConto 2009, 
%    6=deBoer et al., 2014
% - emergence: Threshold to climate curve exposure parameter/'d18O threshold'. 
%    0: no slope to threshold (use as default), 
%    1=slope to threshold (only tested with d18O-curve by Lisiecki&Raymo), 
%    2=two-step threshold changing at time dTTh (hardcoded for now)
% - ss: season, specific to DML samples, only used for data loading, can be
%       removed: 1=MDML1, Heimefrontfjella, 2=MDML2, Jutulstraumen.

% Written by David Lundbek Egholm, Aarhus University
% Modified by Jane Lund Andersen to also i) include [36Cl, 21Ne],
% ii) handle alternate climate inputs, and iii) add a slope or step to the 
% glaciation threshold (emergence)

close all;
addpath('ClimateCurves','Functions')

if ss == 1
    season='MDML1';
elseif ss==2
    season='MDML2';
end

%**** load data *****
load(['./data/',season,'_mn_data.mat'],'sample')

%Cosmogenic halflives
CNprop = getCNprop();

%%%%%% Model setup %%%%%%%%%%

%calculate number of samples input
model.Nsnr = length(snr);

%Number of free depth points
model.Nfree = 3;

%number of sample specific parameters
model.Nsmp = 2*model.Nfree + 2;

% Switch ClimateCurve
model.CCurve=CCurve;
if CCurve == 1
    model.ClimCurve = 'd18Ocurves';
    model.age = 5.0; %Max time (Ma), this curve goes back to 5.32 Ma
    model.mp{1}.name = 'd18O threshold';
    model.mp{1}.vmin = 3.0;
    model.mp{1}.vmax = 5.0;
elseif CCurve == 2
    model.ClimCurve = 'zachos';
    model.age = 15.0; % Max time (Ma), this curve goes back to 67 Ma
    model.mp{1}.name = 'd18O threshold';
    model.mp{1}.vmin = 0;
    model.mp{1}.vmax = 4.98;
elseif CCurve == 3
    model.ClimCurve = 'milankovitch';
    model.age = 15.0; %Max time (Ma), this curve goes back to 20 Ma
    model.mp{1}.name = 'Glaciation threshold';
    model.mp{1}.vmin = 0.0;
    model.mp{1}.vmax = 1.0;       
elseif CCurve >= 4 && CCurve <= 6
    model.ClimCurve = 'iceDMC';
    model.age = 5.0; %Max time (Ma), these curves go back to 5 Ma
    model.mp{1}.name = 'Glaciation threshold';
    model.mp{1}.vmin = 0.0;
    model.mp{1}.vmax = 1.0;   
else
    warning('ClimateCurve not implemented')    
end
    
%initialize
model.z0 = 20; %Max depth
model.Temp = 1.0; %'temperature'; >1 to increase nuclide uncertainties, if=1 no significance, overwritten below
model.Mmp = 2; %number of generic model parameters, is overwritten if emergence>0

%Last deglaciation timing, set wide bounds if this is not well-known. If no
%glacial cover during the LGM is expected, vmax/maxDgla should be = 25 kyr
model.mp{2}.name = 'Time of deglaciation (Tdg)'; 
model.mp{2}.vmin = sample{snr(1)}.minDgla; %Ma, 5e-3
model.mp{2}.vmax = sample{snr(1)}.maxDgla; %Ma, 18e-3

model.emergence=emergence; %0=no slope to threshold
if model.emergence == 1 %1=slope to threshold, only tested with d18O curves
    model.Mmp = 3; %number of generic model parameters overwritten
    model.mp{3}.name = 'slope of d18O threshold';
    model.mp{3}.vmin = -30.0; %permille/Myr
    model.mp{3}.vmax = 30.0;
elseif model.emergence == 2 %2=two-step threshold changing at dTTh, only tested with d18O curves
    model.Mmp = 3; %number of generic model parameters overwritten
    model.mp{3}.name = '2nd d18O threshold';
    model.mp{3}.vmin = 3.0; %[permille]
    model.mp{3}.vmax = 5.0;
    model.dTTh=1; %[Ma] - timing of second threshold, change to variable?
end

%loop samples to set model parameter boundaries for exhumation and get sample info
for i=1:model.Nsnr
    model.mp{model.Mmp+(i-1)*model.Nsmp+1}.name =  ['Z at Tdg, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+1}.vmin = 0;
    model.mp{model.Mmp+(i-1)*model.Nsmp+1}.vmax = 0.1; %m
    model.mp{model.Mmp+(i-1)*model.Nsmp+2}.name =  ['dT2, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+2}.vmin = 0;
    model.mp{model.Mmp+(i-1)*model.Nsmp+2}.vmax = 4.0; %Myr, increase if model.age is long
    model.mp{model.Mmp+(i-1)*model.Nsmp+3}.name =  ['dz2, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+3}.vmin = -1.3;% -1=0.1 m, -1.3=0.05 m
    model.mp{model.Mmp+(i-1)*model.Nsmp+3}.vmax = 1; %log 10 m
    model.mp{model.Mmp+(i-1)*model.Nsmp+4}.name =  ['dT3, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+4}.vmin = 0;
    model.mp{model.Mmp+(i-1)*model.Nsmp+4}.vmax = 5.0; %Myr, increase if model.age is long 
    model.mp{model.Mmp+(i-1)*model.Nsmp+5}.name =  ['dz3, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+5}.vmin = -1.3;% -1=0.1 m, -1.3=0.05 m
    model.mp{model.Mmp+(i-1)*model.Nsmp+5}.vmax = 1; %log 10 m
    model.mp{model.Mmp+(i-1)*model.Nsmp+6}.name =  ['dT4, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+6}.vmin = 0;
    model.mp{model.Mmp+(i-1)*model.Nsmp+6}.vmax = 5.0; %Myr, increase if model.age is long 
    model.mp{model.Mmp+(i-1)*model.Nsmp+7}.name =  ['dz4, sample ',num2str(i)];
    model.mp{model.Mmp+(i-1)*model.Nsmp+7}.vmin = -1.3;% -1=0.1 m, -1.3=0.05 m
    model.mp{model.Mmp+(i-1)*model.Nsmp+7}.vmax = 1.3; %log 10 m
    model.mp{model.Mmp+(i-1)*model.Nsmp+8}.name =  ['E5, sample ',num2str(i)]; %Initial erosion rate
    model.mp{model.Mmp+(i-1)*model.Nsmp+8}.vmin = -2; %log10 -2 =0.01m/Myr
    model.mp{model.Mmp+(i-1)*model.Nsmp+8}.vmax = 2; %2=100 m/Myr

    if (model.mp{model.Mmp+(i-1)*model.Nsmp+2}.vmax + ...
            model.mp{model.Mmp+(i-1)*model.Nsmp+4}.vmax + ...
            model.mp{model.Mmp+(i-1)*model.Nsmp+6}.vmax) >= (model.age - model.mp{2}.vmax)
        error('Check times: dT2+dT3+dT4 should not exceed model.age-Tdg')
    end
    
    %add sample info to models
    model.data{i} = sample{snr(i)};
    
    %number of data per sample (nuclides)
    model.Nds(i) = sample{snr(i)}.Nnuc;

end

%save some general MCMC paramers
model.Nwalk = 10; % number of walkers, reduce for testing purposes
model.burnin = 40e3; % length of burn-in phase
model.Nmod = 400e3; % target number of accepted models per walker
model.Nmax = 100e4; % maximum number of models, stops inversion from running indefinitely if something is wrong

%number of total model parameters
model.Nmp = model.Mmp + model.Nsnr*model.Nsmp;

%resample model parameters into vectors
for i=1:model.Nmp
    umin(i) = model.mp{i}.vmin;
    umax(i) = model.mp{i}.vmax;
    du0(i) = model.mp{i}.vmax - model.mp{i}.vmin;
end
umin = umin(:);
umax = umax(:);

% Load climate curve and check input data depending on Climate Curve used
ccname=[model.ClimCurve,'.mat'];
load(ccname)
if CCurve == 1 %'d18Ocurves by Lisiecki & Raymo, 2005'
    I = find(Age < model.age*1e6);
    model.dO = d18O_4ky(I);
    model.dOt = Age(I)';
    if model.age>5.3
        error('Make sure model.age is < 5.32 Ma')
    end
    if model.mp{1}.vmin<2.65 || model.mp{1}.vmax>5.08
        error('Check d18O boundaries: Should be between 2.65 and 5.08')
    end
elseif CCurve == 2 %'zachos by Zachos et al., 2001'
    I = find(t < model.age*1e6);
    model.dO = dO_4ky(I);
    model.dOt = t(I);
    if model.age>67
        error('Make sure model.age is < 67 Ma')
    end
    if model.mp{1}.vmin<-0.38 || model.mp{1}.vmax>4.98
        error('Check d18O boundaries: Should be between -0.38 and 4.98')
    end
elseif CCurve == 3 %'milankovitch'
    I = find(mil20Ma(:,1) < model.age*1e3);
    model.dO = mil20Ma(I,7);
    model.dOt = mil20Ma(I,1)*1e3;
    if model.age>20
        error('Make sure model.age is < 20 Ma')
    end
    if model.mp{1}.vmin<0 || model.mp{1}.vmax>1
        error('Check d18O boundaries: Should be between 0 and 1')
    end    
elseif CCurve == 4 %'iceDMC Spector et al., 2018'
    I = find(MABSpector(:,1) < model.age*1e3);
    model.dO = MABSpector(I,2);
    model.dOt = MABSpector(I,1)*1e3;
    if model.age>5
        error('Make sure model.age is < 5 Ma')
    end    
    if model.mp{1}.vmin<0 || model.mp{1}.vmax>1
        error('Check d18O boundaries: Should be between 0 and 1')
    end  
elseif CCurve == 5 %'iceDMC: Pollard & DeConto 2009'
    I = find(MABPollard(:,1) < model.age*1e3);
    model.dO = MABPollard(I,2);
    model.dOt = MABPollard(I,1)*1e3;
    if model.age>5
        error('Make sure model.age is < 5 Ma')
    end    
    if model.mp{1}.vmin<0 || model.mp{1}.vmax>1
        error('Check d18O boundaries: Should be between 0 and 1')
    end  
elseif CCurve == 6 %'iceDMC: de Boer et al. 2014'
    I = find(MABdeBoer(:,1) < model.age*1e3);
    model.dO = MABdeBoer(I,2);
    model.dOt = MABdeBoer(I,1)*1e3;
    if model.age>5
        error('Make sure model.age is < 5 Ma')
    end    
    if model.mp{1}.vmin<0 || model.mp{1}.vmax>1
        error('Check d18O boundaries: Should be between 0 and 1')
    end  
end

%data and covariance
for i=1:model.Nsnr %loop samples
    for j=1:sample{snr(i)}.Nnuc %loop nuclides
        nuclide=sample{snr(i)}.nuclides(j); %get nuclide identification
            if nuclide == 1 %10Be
                    dobs((i-1)*model.Nds+j) = model.data{i}.N10;
                    sigd((i-1)*model.Nds+j) = model.data{i}.dN10;
            elseif nuclide == 2 %26Al
                    dobs((i-1)*model.Nds+j) = model.data{i}.N26;
                    sigd((i-1)*model.Nds+j) = model.data{i}.dN26;
            elseif nuclide == 3 %36Cl
                    dobs((i-1)*model.Nds+j) = model.data{i}.N36;
                    sigd((i-1)*model.Nds+j) = model.data{i}.dN36;
            elseif nuclide == 4 %21Ne
                    dobs((i-1)*model.Nds+j) = model.data{i}.N21;
                    sigd((i-1)*model.Nds+j) = model.data{i}.dN21;
            else
                    warning('Nuclide not implemented')
            end
    end
end

Cobs = model.Temp*diag(sigd.^2);
Cobsinv = inv(Cobs);

%%%%%% Initiate random model generation %%%%%%%

%Initialize random sequence
rng('default');

%set walker starting positions
if model.Nwalk > 1 %multiple walkers
    for i=1:model.Nmp
        wini(i,:) = 0.8*(randperm(model.Nwalk)-1)/(model.Nwalk-1)+0.1;
    end
else %only one walker
    wini = randi([0,4],model.Nmp,model.Nwalk)/4;
end

%loop walkers
for nw = 1:model.Nwalk

    %walker starting position - for initial parameter vector
    for i = 1:model.Nmp
        u(i) = (1-wini(i,nw))*model.mp{i}.vmin + wini(i,nw)*model.mp{i}.vmax;
    end     
 
    %initialize
    minres = 1e20; %not used
    res_current = 1e20; %current residual - first model run
    restot = 0; %not used
    acount = 0; %acceptance count
    bcount = 0; %burn-in count
    rcount = 0; %reject count
    accrat = 0; %acceptance rate
    status = zeros(model.Nmax,1); %model status
    erosion_rec = zeros(model.Nmax,1); %not used
    up_rec = zeros(model.Nmp,model.Nmax); %proposed parameters
    u_rec = zeros(model.Nmp,model.Nmax);
    N10_rec = zeros(model.Nmax,1); %not used
    N26_rec = zeros(model.Nmax,1); %not used
    restot_rec = zeros(model.Nmax,1);
    accrat_rec = zeros(model.Nmax,1);
    k_rec = zeros(model.Nmax,1);
    duR = zeros(model.Nmp,1); %not used
   
    accfac = 1e-2;
    ktarget = 0.025; %not used
    
    %run models
    mi = 0; %model iteration
    k = 0.01; %initial step length

    while ((mi < model.Nmax)&&(acount < model.Nmod))
        
        mi = mi + 1; %update model iteration

        disp(['nw = ',num2str(nw),'/',num2str(model.Nwalk),' mi = ',num2str(mi),'/',num2str(model.Nmax),' bcount = ',num2str(bcount),' acount = ',num2str(acount),' accrat = ',num2str(accrat),' k = ',num2str(k)]);
    
       %***** step length ******
       
       %set acceptance ratio
        if (mi > 100)
            accrat = (sum(abs(status((mi-100):(mi-1))))+1)/100; %update based on status of last 100 models
        elseif (bcount < model.burnin)
            accrat = 0.1; %first 100 models of burnin
        else 
            accrat = 0.3; %only gets here if model.burnin is <100?
        end
        
        %steplength
        if (bcount < 0.5*model.burnin) %first half of burn-in
            
            k = k*((1-accfac) + accfac*accrat/0.1);
            
            %take larger steps if too few models are accepted intially
            if ((mi > 100)&&(bcount < 2)) k = 1.0; %limited to 0.5 below
            elseif ((mi > 10)&&(bcount < 2)) k = 0.01*mi;
            end
            
        elseif (bcount < model.burnin) %second half of burn-in

            k = k*((1-accfac) + accfac*accrat/0.2);
           
            if ((mi > 100)&&(bcount < 2)) k = 1.0; %Does not enter here unless burn-in is <200?
            elseif ((mi > 10)&&(bcount < 2)) k = 0.01*mi; %and here if burn-in is <20?
            end
            
        
        elseif (acount < model.Nmod) %after burn-in

            k = k*((1-accfac) + accfac*accrat/0.3);    
            
        end  
        
        if (k > 0.5) k = 0.5; end %steplength limited to 0.5
        
        %model 'temperature' gradually reduced from 11 to 1 during burn-in
        if (bcount < model.burnin) model.Temp = 1.0 + 10.0*(model.burnin-bcount)/model.burnin;
        else model.Temp = 1.0;
        end
        
        %********* propose new parameters ************
        
         %random step
        du = 0.5*randn(model.Nmp,1).*du0(:);
        
        %proposed model
        up = u(:) + k*du(:);
     
        
        %retake if outside parameter boundaries
        while (any(up(:) < umin(:))||any(up(:) > umax(:))) %||((up(2)+up(model.Mmp+2)+up(model.Mmp+4)+up(model.Mmp+6)) > model.age)
        
            %random step
            du = 0.5*randn(model.Nmp,1).*du0(:);

            %proposed model
            up = u(:) + k*du(:);
            
        end
        
        %********** Forward model *****************
        % the ismember(#,NNuc)implementation within one forward-model was 
        % very slow, so I have changed to this clunkier solution with 
        % multiple forward-models for now
        if model.Nds(1) == 2 %10Be, 26Al
            [gmp] = forward_bedrockvJ1_mn4_BeAl(up,model,CNprop);  
        elseif model.Nds(1) == 3 && model.data{1}.nuclides(3) == 3 %10Be, 26Al, 36Cl
            [gmp] = forward_bedrockvJ1_mn4_BeAlCl(up,model,CNprop);
        elseif model.Nds(1) == 3 && model.data{1}.nuclides(3) == 4 %10Be, 26Al, 21Ne
            [gmp] = forward_bedrockvJ1_mn4_BeAlNe(up,model,CNprop);
        else
            [gmp] = forward_bedrockvJ1_mn4_BeAlClNe(up,model,CNprop); %10Be, 26Al, 36Cl, 21Ne
        end
            
        %Acceptance critieria
        %restot = (dobs(:)-gmp(:))'*Cobsinv*(dobs(:)-gmp(:));
        restot = (dobs(:)-gmp(:))'*Cobsinv*(dobs(:)-gmp(:))/model.Temp;
        rfrac = exp(-0.5*restot)/exp(-0.5*res_current);
        alpha = rand(1);
    
    
        %if model is accepted
        if ((alpha < rfrac)||(mi == 1))
        
            u = up;
            gm = gmp;
            res_current = restot;
        
            %accepted after burn-in
            if (bcount > model.burnin)
        
                status(mi) = 1;
                acount = acount + 1;
                
            %if accepted during burn-in    
            else
        
                status(mi) = -1;
                bcount = bcount + 1;
            
            end
        
        %rejected
        else
        
            status(mi) = 0;
            rcount = rcount + 1;                 
    
        end
       
        %save things
        up_rec(:,mi) = up(:);
        u_rec(:,mi) = u(:);
        gm_rec(:,mi) = gm(:);
        restot_rec(mi) = res_current;
        accrat_rec(mi) = accrat;
        k_rec(mi) = k;
                
    end

    %change status for models that are rejected during burn-in
    Imin=find(status == 1,1); %find index of first accepted model after burn-in
    I = find(status(1:Imin) == 0); %find indices of rejected models in burn-in
    status(I) = -2; %set status of these models to -2
    
    %save things
    model.walker{nw}.status = status(1:mi);
    model.walker{nw}.up = up_rec(:,1:mi);
    model.walker{nw}.u = u_rec(:,1:mi);
    model.walker{nw}.gm = gm_rec(:,1:mi);
    model.walker{nw}.restot = restot_rec(1:mi);
    model.walker{nw}.acount = acount;
    model.walker{nw}.bcount = bcount; 
    model.walker{nw}.rcount = rcount;
    model.walker{nw}.accrate = accrat_rec(1:mi);
    model.walker{nw}.kstep = k_rec(1:mi);

end

%%save output
%sample numbers
str = num2str(snr(1));
if length(snr) > 1
    for i=2:length(snr)
        str = [str,['-',num2str(snr(i))]];
    end
end
%filename
if model.emergence==0
    if CCurve <=3 || CCurve >=7
        savefile = ['models/MDML/MC_sample_',str,'_mn5_',model.ClimCurve,'_',season,'.mat'];
    else 
        savefile = ['models/MDML/MC_sample_',str,'_mn5_',model.ClimCurve,num2str(CCurve-3),'_',season,'.mat'];
    end
elseif model.emergence==1
    savefile = ['models/MDML/MC_sample_',str,'_mn5_',model.ClimCurve,'_',season,'_emergence.mat'];
elseif model.emergence==2
    savefile = ['models/MDML/MC_sample_',str,'_mn5_',model.ClimCurve,'_',season,'_twostep.mat'];
end
%save
save(savefile,'model');
