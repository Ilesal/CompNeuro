% Code Tested on Matlab R2020b 
% Assignment 1 in Cortical Modelling




%section A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define the hyperparameters, which are fixed during the simulation

tau = 15*10^-3; %characteristic time in seconds
dt = tau/50; %timestep divided by tau to obtain more accurate results
Rm = 10*10^6;  %total membrane resistance in omega (Ω)
V_thresh= -50*10^-3; %threshold in Volt (V)
V_reset = -80*10^-3; %reset membrane potential in Volt (V), which must be lower than the resting potential
E= -70*10^-3; %resting potential in Volt (V) 
Ie =3.5*10^-9; %the external positive current injected into the neuron, in Ampere

%section B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialise the necessary variables that will change with time
 
T= 0:dt:0.3; %time vector (maximum time 300ms) with little increments of dt
V_hat = zeros(size(T)); %inizialise the membrane potential variable V_hat, which will store the estimate of our membrane potential through time
S= zeros(size(T)); %inizialise a spike array discretised with the same size as the time vector T
V_0 = -70*10^-3; %initial membrane potential in Volt (V) equal to the resting potential E
V_hat(1) = V_0; %The first value of V_hat is equal to the initial membrane potential value, ie -70*10^-3


%section C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=2:length(T)  %start from index 2 to compute the next step as a function of the previous step
    if V_hat(t-1)<V_thresh % if the membrane potential value is lower than threshold, go ahead with the differential equation in the next step
        V_hat(t)= V_hat(t-1)+(dt/tau)*(E-V_hat(t-1)+Rm*Ie); %the value of V_hat at index t is the value of V at the previous timestep plus the value of the derivative
    else
        V_hat(t)= V_reset; %when the threshold is crossed (and the neuron fire a spike) reset the membrane potential to -80*10^-3
        S(t)=1;   %store the spike at the particular timestep. If there is a spike store 1, otherwise 0
    end    
end


figure, plot(T,V_hat); %plot the membrane potential and time vector T
xlabel('Time [s]'); %add a x-axis label
ylabel('Membrane Potential [V]'); %add a y-axis label
title('LIF Neuron'); %add the tile

N_spikes=sum(S);% Count the total number of spikes:
% Calculate the firing rate by dividing the total number of spikes with the
% total length of the simulation, ie 300ms:
firing_rate= N_spikes/max(T); %In Hertz
disp(['The firing rate is ' num2str(firing_rate),[' Hz']]); %display the firing rate


%% Second simulation with increased tau

close all
clear all 

tau = 30*10^-3; %characteristic time in seconds
dt = tau/50; %timestep divided by tau to obtain more accurate results
Rm = 10*10^6;  %total membrane resistance in omega (Ω)
V_thresh= -50*10^-3; %threshold in Volt (V)
V_reset = -80*10^-3; %reset membrane potential in Volt (V), which must be lower than the resting potential
E= -70*10^-3; %resting potential in Volt (V) 
Ie =3.5*10^-9; %the external current injected into the neuron


%section B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialise the necessary variables that will change with time

T= 0:dt:0.3; %time vector (maximum time 300ms) with little increments of dt
V_hat = zeros(size(T)); %inizialise the membrane potential variable V_hat, which will store the estimate of our membrane potential through time
S= zeros(size(T)); %inizialise a spike array discretised following the same size as the time vector T
V_0 = -70*10^-3; %initial membrane potential in Volt (V)
V_hat(1) = V_0; %The first value of V_hat is equal to the initial membrane potential value, ie -70*10^-3


%section C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=2:length(T)  %start from index 2 to compute the next step as a function of the previous step
    if V_hat(t-1)<V_thresh % if the membrane potential value is lower than threshold, go ahead with the differential equation in the next step
        V_hat(t)= V_hat(t-1)+(dt/tau)*(E-V_hat(t-1)+Rm*Ie); %the value of V_hat at index t is the value of V at the previous timestep plus the value of the derivative
    else
        V_hat(t)= V_reset; %when the threshold is crossed (and the neuron fire a spike) reset the membrane potential to -80*10^-3
        S(t)=1;   %store the spike at the particular timestep. If there is a spike store 1, otherwise 0
    end    
end


figure, plot(T,V_hat); %plot the membrane potential and time vector T
xlabel('Time [s]'); %add a x-axis label
ylabel('Membrane Potential [V]'); %add a y-axis label
title('LIF Neuron'); %add the tile


N_spikes=sum(S); % Count the total number of spikes:
% Calculate the firing rate by dividing the total number of spikes with the
% total length of the simulation, ie 300ms:
firing_rate= N_spikes/max(T); %In Hertz
disp(['The firing rate is ' num2str(firing_rate),[' Hz']]); %display the firing rate

%% Third simulation with a decreased V_thresh value

close all
clear all 


tau = 15*10^-3; %characteristic time in seconds
dt = tau/50; %timestep divided by tau to obtain more accurate results
Rm = 10*10^6;  %total membrane resistance in omega (Ω)
V_thresh= -60*10^-3; %new threshold value
V_reset = -80*10^-3; %reset membrane potential in Volt (V), which must be lower than the resting potential
E= -70*10^-3; %resting potential in Volt (V) 
Ie =3.5*10^-9; %the external positive current injected into the neuron

%section B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialise the necessary variables that will change with time


T= 0:dt:0.3; %time vector (maximum time 300ms) with little increments of dt
V_hat = zeros(size(T)); %inizialise the membrane potential variable V_hat, which will store the estimate of our membrane potential through time
S= zeros(size(T)); %inizialise a spike array discretised following the same size as the time vector T
V_0 = -70*10^-3; %initial membrane potential in Volt (V)
V_hat(1) = V_0; %The first value of V_hat is equal to the initial membrane potential value, ie -70*10^-3


%section C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=2:length(T)  %start from index 2 to compute the next step as a function of the previous step
    if V_hat(t-1)<V_thresh % if the membrane potential value is lower than threshold, go ahead with the differential equation in the next step
        V_hat(t)= V_hat(t-1)+(dt/tau)*(E-V_hat(t-1)+Rm*Ie); %the value of V_hat at index t is the value of V at the previous timestep plus the value of the derivative
    else
        V_hat(t)= V_reset; %when the threshold is crossed (and the neuron fire a spike) reset the membrane potential to -80*10^-3
        S(t)=1;   %store the spike at the particular timestep. If there is a spike store 1, otherwise 0
    end    
end



figure, plot(T,V_hat); %plot the membrane potential and time vector T
xlabel('Time [s]'); %add a x-axis label
ylabel('Membrane Potential [V]'); %add a y-axis label
title('LIF Neuron'); %add the tile


N_spikes=sum(S); % Count the total number of spikes:
% Calculate the firing rate by dividing the total number of spikes with the
% total length of the simulation, ie 300ms:
firing_rate= N_spikes/max(T); %In Hertz
disp(['The firing rate is ' num2str(firing_rate),[' Hz']]); %display the firing rate

