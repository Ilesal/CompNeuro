%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time Series Analysis: frequency-domain directional measures of connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% In this report I will analyse local-field potentials (LFP) recorded from one Parkinson’s Disease patient in the On and Off medication state from four brain areas: 
% subthalamic nucleus (STN), globus pallidus internus (GPI), thalamus (Th) and premotor cortex (PM). 
% Using frequency-domain directional measures of connectivity such as Spectral Granger Causality (SGC) and Partial Directed Coherence (PDC) I will investigate the directional influences between these areas.

% Code Tested on Matlab R2021a 


close all
clear all

%% Add fieldtrip to your path  

% addpath(genpath('~YOURPAHTH/fieldtrip_YOURVERSION/'))
% ft_detaults

%% load dataset

load data_PD_OFF.mat %load the dataset from condition OFF
data = dataoff;  

%% load second dataset

% load data_PD_ON.mat %load the dataset from condition ON
% data = dataon;  

%% Plot the four timeseries from the first trial of the recording

% The simulated data consists of 4 channels, 500 trials of 1 second length, sampled at 200Hz.
% You can easily visualize the data for example in the first trial using

figure
plot(data.time{1}, data.trial{1}) %plot the first trial
legend(data.label) %add legend
xlabel('time (s)') %add label on the x-axis
%title('OFF Condition') %uncomment to add title
%title('ON Condition')



%% Obtaining frequency-domain directional measures of connectivity

%% Approach 1: Computation of the multivariate autoregressive model (MVAR)  
%
% To be able to compute spectrally resolved Granger causality, or other frequency-domain directional measures of connectivity, 
% we have to fit an autoregressive model to the data. This is done using the ft_mvaranalysis function.
% 
% For the actual computation of the autoregressive coefficients FieldTrip makes use of an implementation from third party toolboxes. 
% At present ft_mvaranalysis supports the biosig and bsmart toolboxes for these computations.
% 
% In this assessment we will use the bsmart toolbox. 
% The relevant functions have been included in the FieldTrip release in the fieldtrip/external/bsmart directory.

%Parametric Route
cfg         = [];
cfg.order   = 5;   % model order: 5 time-lags
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, data); %this variable contains a description of the data in terms of a multivariate autoregressive model

 
% Computing connectivity measures requires first to estimate the spectral transfer function

cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata); % The resulting mfreq data structure contains the pairwise transfer function 
% between the 4 channels for 101 frequencies.


%% Compute the spectral Granger-causality index from 
% the spectral transfer function 

cfg           = [];
cfg.method    = 'granger';
mgranger       = ft_connectivityanalysis(cfg, mfreq); 


% Visualisation: Directions of influence go 
% from row to column, as indicated in the figure
cfg           = [];
cfg.parameter = 'grangerspctrm';
cfg.jackknife  = 'yes'  
cfg.zlim      = [0 1.5];
figure('Name','MVAR Granger causality')
ft_connectivityplot(cfg, mgranger);
 
%Alternatively, to see the frequencies involved:

figure
for row=1:4
for col=1:4
  subplot(4,4,(row-1)*4+col);
  plot(mgranger.freq, squeeze(mgranger.grangerspctrm(row,col,:)))
  ylim([0 1.5])
end
end
 

%% Compute the Partial Directed Coherence (PDC) index - Parametric (MVAR)
% from the spectral transfer function 

cfg           = [];
cfg.method    = 'pdc';
mpdc       = ft_connectivityanalysis(cfg, mfreq);  

% Visualisation
cfg           = [];
cfg.parameter = 'pdcspctrm';
cfg.jackknife  = 'yes'  
cfg.zlim      = [0 1.5];
figure('Name','MVAR PDC')
ft_connectivityplot(cfg, mpdc);

%to get a more detailed description of the freqeuncies involved
figure
for row=1:4
for col=1:4
  subplot(4,4,(row-1)*4+col);
  plot(mpdc.freq, squeeze(mpdc.pdcspctrm(row,col,:)))
  ylim([0 1.5])
end
end
 


%% Approach 2: Multitaper frequency decomposition method - Non-parametric
% 
%  It is also possible to compute the spectral transfer function using 
%  non-parametric spectral factorization of the cross-spectral density matrix. 
%  For this, we need a Fourier-type decomposition of the data. Here we will
%  use multitaper transformation (not covered in class but it's similar to
%  wavelet transformation, with basis functions (tapers) based on
%  discrete prolate spheroidal sequences (DPSS), also known as the Slepian
%  sequence.

cfg           = [];
cfg.method    = 'mtmfft'; % implements multitaper frequency transformation
cfg.taper     = 'dpss';  
cfg.output    = 'fourier';  %return the complex Fourier-spectra
cfg.tapsmofrq = 2;  %the amount of spectral smoothing through
%                    multi-tapering. Note that 2 Hz smoothing means
%                    plus-minus 2 Hz, i.e. a 4 Hz smoothing box.                  

freq          = ft_freqanalysis(cfg, data);  %NON_PARAMETRIC         



% Spectral Granger causality index
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq);    

% Visualisation
cfg           = [];
cfg.parameter = 'grangerspctrm';
cfg.zlim      = [0 1.5];
figure('Name','Taper Granger causality')
ft_connectivityplot(cfg, granger);

%Alternatively, to get a better understanding:

figure
for row=1:4
for col=1:4
  subplot(4,4,(row-1)*4+col);
  plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)))
  ylim([0 1.5])
end
end
 
 

%% Compute the Partial Directed Coherence (PDC) directly from the
% frequency-transformed data - Non-parametric

cfg           = [];
cfg.method    = 'pdc';
pdc       = ft_connectivityanalysis(cfg, freq);    

% Visualisation
cfg           = [];
cfg.parameter = 'pdcspctrm';
cfg.zlim      = [0 1];
figure('Name','Taper PDC')
ft_connectivityplot(cfg, pdc); colormap(1-hot)

%to get a more detailed description of the freqeuncies involved
figure
for row=1:4
for col=1:4
  subplot(4,4,(row-1)*4+col);
  plot(pdc.freq, squeeze(pdc.pdcspctrm(row,col,:)))
  ylim([0 1.5])
end
end


%% Compare visually Granger-causality indexes obtained with both methods: MVAR and Multitaper transformation
% Write your own code

% Visualisation
cfg           = [];
cfg.parameter = 'grangerspctrm';
cfg.zlim      = [0 1.5];
figure('Name','Granger causality')
ft_connectivityplot(cfg, mgranger,granger);


%Alternatively, to see the frequencies in the x-axis: 

figure
for row=1:4
for col=1:4
  subplot(4,4,(row-1)*4+col);
  plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)),'k')  
  hold on
  plot(mgranger.freq, squeeze(mgranger.grangerspctrm(row,col,:)),'m')  
  ylim([0 1.5])
  
end
end
legend('Granger: MVAR', 'Granger: MT')



%% Compare visually PDC indexes obtained with both methods: MVAR and Multitaper transformation


% Visualisation
cfg           = [];
cfg.parameter = 'pdcspctrm';
cfg.zlim      = [0 1.5];
figure('Name','PDC')
ft_connectivityplot(cfg, mpdc, pdc);


%Alternatively, to see the frequencies in the x-axis: 

figure
for row=1:4
for col=1:4
  subplot(4,4,(row-1)*4+col);
  plot(pdc.freq, squeeze(pdc.pdcspctrm(row,col,:)),'m')  
  hold on
  plot(mpdc.freq, squeeze(mpdc.pdcspctrm(row,col,:)),'k')  
  ylim([0 1.5])
  
end
end
legend('PDC: MVAR', 'PDC: MT')



%% Compute the coherence  (additional analysis)


cfg          = [];
cfg.method   = 'coh';
coh          = ft_connectivityanalysis(cfg, freq); %non-parametric
cohm         = ft_connectivityanalysis(cfg, mfreq); %parametric MVAR

% visualisation 

cfg           = [];
cfg.parameter = 'cohspctrm';
cfg.zlim      = [0 1];
ft_connectivityplot(cfg, coh, cohm); %plot them together


%Plot the parametric and non-parametric coherence together to see the
%frequencies:
figure,plot(cohm.freq,squeeze(cohm.cohspctrm(1,2,:)),'b'); hold on
plot(cohm.freq,squeeze(cohm.cohspctrm(1,3,:)),'c');
plot(cohm.freq,squeeze(cohm.cohspctrm(1,4,:)),'r');
plot(cohm.freq,squeeze(cohm.cohspctrm(2,3,:)),'k');
plot(cohm.freq,squeeze(cohm.cohspctrm(2,4,:)),'g');
plot(cohm.freq,squeeze(cohm.cohspctrm(3,4,:)),'m');
plot(coh.freq,squeeze(coh.cohspctrm(1,2,:)),'b');
plot(coh.freq,squeeze(coh.cohspctrm(1,3,:)),'c');
plot(coh.freq,squeeze(coh.cohspctrm(1,4,:)),'r');
plot(coh.freq,squeeze(coh.cohspctrm(2,3,:)),'k');
plot(coh.freq,squeeze(coh.cohspctrm(2,4,:)),'g');
plot(coh.freq,squeeze(coh.cohspctrm(3,4,:)),'m');
legend('1->2','1->3','1->4','2->3','2->4','3->4');
%title('Coherence Spectrum OFF-Therapy');
%title('Coherence Spectrum ON-Therapy');
xlabel('Frequency');
ylabel('Coherence');

 
%% Phase-slope index (PSI) (additional analysis)

%Parametric PSI (with MVAR)
cfg           = [];
cfg.method    = 'psi';
cfg.bandwidth = 4;
psi1 = ft_connectivityanalysis(cfg, mfreq);

%Visualise
 
%Visualise
figure;plot(psi1.freq,squeeze(psi1.psispctrm(1,2,:))); hold on;
plot(psi1.freq,squeeze(psi1.psispctrm(1,3,:)));
plot(psi1.freq,squeeze(psi1.psispctrm(1,4,:)));
plot(psi1.freq,squeeze(psi1.psispctrm(2,1,:)));  
plot(psi1.freq,squeeze(psi1.psispctrm(2,3,:)));
plot(psi1.freq,squeeze(psi1.psispctrm(2,4,:)));
plot(psi1.freq,squeeze(psi1.psispctrm(3,1,:)));
plot(psi1.freq,squeeze(psi1.psispctrm(3,2,:)));
plot(psi1.freq,squeeze(psi1.psispctrm(3,4,:)));
plot(psi1.freq,squeeze(psi1.psispctrm(4,1,:)));
plot(psi1.freq,squeeze(psi1.psispctrm(4,2,:)));
plot(psi1.freq,squeeze(psi1.psispctrm(4,3,:)));
legend('1->2','1->3','1->4','2->1','2->3','2->4','3->1','3->2','3->4','4->1','4->2','4->3');
% title('Phase Slope Index (PSI) Parametric OFF-Therapy'); 
% title('Phase Slope Index (PSI) Parametric ON-Therapy'); 
xlabel('Frequency');
ylabel('Phase Slope');


%%PSI Non-Parametric analysis

% cfg           = [];
% cfg.method    = 'psi';
% cfg.bandwidth = 4;
% psi2 = ft_connectivityanalysis(cfg, freq);

%Visualise

% figure;plot(psi2.freq,squeeze(psi2.psispctrm(1,2,:))); hold on;
% plot(psi2.freq,squeeze(psi2.psispctrm(1,3,:)));
% plot(psi2.freq,squeeze(psi2.psispctrm(1,4,:)));
% plot(psi2.freq,squeeze(psi2.psispctrm(2,1,:)));  
% plot(psi2.freq,squeeze(psi2.psispctrm(2,3,:)));
% plot(psi2.freq,squeeze(psi2.psispctrm(2,4,:)));
% plot(psi2.freq,squeeze(psi2.psispctrm(3,1,:)));
% plot(psi2.freq,squeeze(psi2.psispctrm(3,2,:)));
% plot(psi2.freq,squeeze(psi2.psispctrm(3,4,:)));
% plot(psi2.freq,squeeze(psi2.psispctrm(4,1,:)));
% plot(psi2.freq,squeeze(psi2.psispctrm(4,2,:)));
% plot(psi2.freq,squeeze(psi2.psispctrm(4,3,:)));
% legend('1->2','1->3','1->4','2->1','2->3','2->4','3->1','3->2','3->4','4->1','4->2','4->3');
% %title('Phase Slope Index (PSI) Non-Parametric OFF-Therapy'); 
% %title('Phase Slope Index (PSI) Non-Parametric ON-Therapy'); 
% xlabel('Frequency');
% ylabel('Phase Slope');


%% Spectral analysis (additional analysis)
%calculate the fourier coefficients (non-parametric derivation of power)

cfg = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss'; 
cfg.output    = 'fourier';
cfg.foilim    = [0 60]; %choose the frequency ranges
cfg.tapsmofrq = 5;  
freq          = ft_freqanalysis(cfg, data);
fd            = ft_freqdescriptives(cfg, freq); %freqdescriptives calculates the power spectrum  


figure;plot(fd.freq, fd.powspctrm);
title('Non-Parametric Power Spectrum');
legend(data.label)

%% Obtain the autocorrelation function for each of the four time series

% Concatenate data epochs
X = []; 

for n = 1:500
X = cat(2, X, data.trial{n});
end


% Estimate ACF and visualise
for nchan = 1 : 4
    figure;autocorr(X(nchan,:));
    
%     %uncomment the following lines to add title
%     capt = sprintf('Signal %d OFF-Therapy', nchan);
%     capt = sprintf('Signal %d ON-Therapy', nchan);
%     title(capt); 
    
end




