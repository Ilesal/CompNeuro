% Code Tested on Matlab R2021a 
% Mini-project in Cortical Modelling
% First Part

%%  Implement a perfect integrator and fire neuron
 
close all
clear all


% Define the unit-less parameters

V_thresh=1; % membrane threshold
Rm=1;       % membrane resistance
tau=1;      % time constant
dt=0.01;    % timestep
E=0;        % resting potential
V_reset=-1; % membrane reset
Ie=4;       % current injected
V_0=0;      % initial membrane potential
t_stop=2.5; % max time of simualation, 2.5 unit of time


% Initialise the variables that will chnage with time

T= 0:dt:t_stop;             % define the time vector
V_hat = zeros(size(T));     % inizialise the membrane potential variable 
V_hat(1) = V_0;             % The first value of V_hat is equal to the initial membrane potential value
S=zeros(size(T));           % inizialise a spike array discretised with the same size as the time vector T


for i=2:length(T)           % start from index 2 to compute the next step as a function of the previous step    
    if V_hat(i-1)<V_thresh  % if the membrane potential value is lower than threshold, go ahead with the differential equation in the next step
     V_hat(i) = V_hat(i-1)+Ie*dt; %he value of V_hat at index t is given by the value of V at the previous timestep plus the current times the timestep
    else
        V_hat(i)= V_reset;  % when the threshold is crossed (and the neuron fire a spike) reset the membrane potential to -1        
        S(i)=1;      % store the spike at the particular timestep. If there is a spike store 1, otherwise 0
    end
end

%Find the spikes
n_spikes=find(S)*0.01;

%Plot the membrane potential as a function of time  
figure, plot(T,V_hat);  
xline([n_spikes],'-r','Linewidth',1) %draw a red vertical line for each spike
yline(V_thresh,'--m'); %draw an horizontal line for the threshold
yline(V_reset,'--c'); %draw an horizontal line for the membrane reset value
ylim([-1.5 1.5])  %change y-axis limit
text(0.28,1.07,'Threshold','Color','magenta','FontSize',12) %add line description
text(0.28,-1.07,'Reset','Color','cyan','FontSize',12) %add line description
xlabel('Time'); %add a x-axis label
ylabel('Membrane Potential'); %add a y-axis label
title('Perfect Integrate-Fire Neuron'); %add the tile

%Measure the ISI and CV 

isi=diff(n_spikes); %calculate interspike interval
cv=std(isi)/mean(isi); %calculate the coefficient of variation, to quantify the variability

figure, histogram(isi,20); %Plot the ISI histogram 
xlabel('Interspike Interval'); %add x-label
ylabel('Frequency'); %add ylabel
text(0.75,3.5,['ISI = ' num2str(mean(isi))],'FontSize',12) %add ISI mean value on the plot
text(0.75,3.3,['CV = ' num2str(cv)],'FontSize',12) %add CV value on the plot
xline((2/Ie),'--r'); %plot with a red vertical dotted line the exact solution ISI=2/Ie

firing_rate= length(n_spikes)/max(T); %quantify the firing rate
disp(['The firing rate is ' num2str(firing_rate)]); %display the firing rate
disp(['The mean ISI is ' num2str(mean(isi))]); %display the mean ISI
disp(['The cv is ' num2str(cv)]); %display the cv


%The firing rate is the inverse of ISI, thus we can also calculate it by
%using 1/ISI:
%firing_rate=round(1/mean(isi));