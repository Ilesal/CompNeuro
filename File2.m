% Code Tested on Matlab R2021a 
% Second Part

%% Question 2.2, 2.3, 2.4
% Implement a noisy integrator and fire neuron
% dV/dt=Is+ε(t) where ε is Gaussian noise

close all
clear all

% Define the parameters

V_thresh=1; % membrane threshold
Rm=1;       % membrane resistance
tau=1;      % time constant
dt=0.01;    % timestep
E=0;        % resting potential
V_reset=-1; % membrane reset
Ie=4;       % current injected
V_0=0;      % initial membrane potential
t_stop=2.5; % max time of simualation, 2.5 unit of time


% Initialise the variables that will change with time

T= 0:dt:t_stop;             % define the time vector
V_hat = zeros(size(T));     % inizialise the membrane potential variable 
V_hat(1) = V_0;             % The first value of V_hat is equal to the initial membrane potential value
S=zeros(size(T));           % inizialise a spike array discretised with the same size as the time vector T
 

          for t=2:length(T) %we want to compute the next step as a function of the previous step
             if V_hat(t-1)<V_thresh %if membrane potential is lower than threshold, we follow the differential equation. As soon as it cross the threshold we have a spike
                epsilon = normrnd (0,0.5*Ie); %noise from a gaussian distribution with mean 0 and std=5
                V_hat(t)= V_hat(t-1)+(Ie+epsilon)*dt; %integrate with euler the differential equation
              else
                  V_hat(t)=V_reset; %neuron fire a spike and we reset to -1
                   S(t)=1; %store the spike fired at the particular timestep

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
title('Noisy Integrate-Fire Neuron'); %add the tile

%Measure the ISI and CV (Question 2.4)

isi=diff(n_spikes); %calculate interspike interval
cv=std(isi)/mean(isi); %calculate the coefficient of variation, to quantify the variability

figure, histogram(isi); %Plot the ISI histogram 
xlim([0 1]) %define x-axis limit
xlabel('Interspike Interval'); %add x-label
ylabel('Frequency'); %add ylabel
text(0.75,3.5,['ISI = ' num2str(mean(isi))],'FontSize',12) %add ISI mean value on the plot
text(0.75,3.3,['CV = ' num2str(cv)],'FontSize',12) %add CV value on the plot
xline((2/Ie),'--r'); %plot with a red vertical dotted line the exact solution ISI=2/Ie


firing_rate= length(n_spikes)/max(T); %quantify the firing rate
disp(['The firing rate is ' num2str(firing_rate)]); %display the firing rate
disp(['The mean ISI is ' num2str(mean(isi))]); %display the mean ISI

    


%% Question 2.5, 2.6, 2.7, 2.8
% In the first trial the IF neuron receives Ie = 4 and in a second trial the neuron receives Ie = 4.2

close all
clear all

% Define the parameters

V_thresh=1; % membrane threshold
Rm=1;       % membrane resistance
tau=1;      % time constant
dt=0.01;    % timestep
E=0;        % resting potential
V_reset=-1; % membrane reset
Ie=[4, 4.2];% current injected in the two trials
V_0=0;      % initial membrane potential
t_stop=2.5; % max time of simualation, 2.5 unit of time


% Initialise the variables that will change with time

T= 0:dt:t_stop;   % define the time vector
n_trials=2;       % number of trials


    for i=1:n_trials  %Question 2.5, run two trials
        I=Ie(i);      % in the first trial I=4 and in the second trial I=4.2
        V_hat= zeros(size(T)); % inizialise the membrane potential variable 
        V_hat(1)=V_0;      % the first value of V_hat is equal to the initial membrane potential value
        S= zeros(size(T));  % inizialise the spike array
    
          for t=2:length(T) %we want to compute the next step as a function of the previous step
             if V_hat(t-1)<V_thresh %if membrane potential is lower than threshold, we follow the differential equation. However, as soon as it cross the threshold we have a spike
                epsilon = normrnd (0,0.5*I); %noise from a gaussian distribution
                V_hat(t)= V_hat(t-1)+(I+epsilon)*dt; %integrate with euler the differential equation
              else
                  V_hat(t)=V_reset; %neuron fire a spike and we reset to -1
                   S(t)=1; %store the spike at the particular time step

                end
            end
            
    V_hat_trials(i,:)=V_hat; %store the V_hat for each trial 
    N_spikes=find(S)*dt; %find the number of spikes
    firing_rates(i)= length(N_spikes)/T(end); 
    isi(i)={diff(find(S)*dt)}; %Question 2.5, store ISI for each trial in two arrays
    mean_isi(i)=mean(isi{i}); %Question 2.7, calculate the mean ISI for each trial
    std_isi(i)=std(isi{i}); % calculate std ISI for each trial
    cv(i)=std_isi/mean_isi; %calculate the coefficient of variation for each trial
%     varian(i)=var(isi{i}); %uncomment to calculate the variance of ISI
%     fano(i)=varian(i)/mean_isi(i); %uncomment to calculate the Fano Factor as variance(ISI)/mean(ISI)
    end
   
    
%Plot the membrane potentials of the two trials
figure, plot(T,V_hat_trials(1,:)), hold on, plot(T,V_hat_trials(2,:))  
ylim([-1.5 1.5]); %define x-axis limit
yline(V_thresh,'r--'); %draw an horizontal line for the membrane threshold value
yline(V_reset,'--c'); %draw an horizontal line for the membrane reset value
legend('Trial 1','Trial 2', 'Vthresh','Vreset')
title('IF Neuron with Noise')

%Question 2.6
%Plot the Interspike interval for the two trials in one plot 
figure()
subplot (2,1,1), histogram(isi{1}) %plot ISI for the first stimulus
xlim([0 1]);%define x-axis limit
ylabel('Frequency') %add y-axis label
title('First trial with Ie=4'); %add title
subplot (2,1,2), histogram(isi{2}) %plot ISI for the second stimulus stimulus
xlim([0 1]); %define x-axis limit
title('Second trial with Ie=4.2');%add title
xlabel('Interspike Interval'); %add x-axis label
ylabel('Frequency'); %add y-axis label



%% Question 2.9 and 2.10 - Run the same simulation for 2000 units of time
 
close all
clear all

% Define the unit-less parameters

V_thresh=1; % membrane threshold
Rm=1;       % membrane resistance
tau=1;      % time constant
dt=0.01;    % timestep
E=0;        % resting potential
V_reset=-1; % membrane reset
Ie=[4, 4.2];% current injected in the two trials
V_0=0;      % initial membrane potential
t_stop=2000; % max time of simualation, 2000 unit of time, question 2.9


% Initialise the variables that will change with time

T= 0:dt:t_stop;   % define the time vector
n_trials=2;       % number of trials


    for i=1:n_trials  %run two trials
        I=Ie(i);      % in the first trial I=4 and in the second trial I=4.2
        V_hat= zeros(size(T)); % inizialise the membrane potential variable 
        V_hat(1)=V_0;      % the first value of V_hat is equal to the initial membrane potential value
        S= zeros(size(T));  % inizialise the spike array
    
          for t=2:length(T) %we want to compute the next step as a function of the previous step
             if V_hat(t-1)<V_thresh %if membrane potential is lower than threshold, we follow the differential equation. However, as soon as it cross the threshold we have a spike
                epsilon = normrnd (0,0.5*I); %noise from a gaussian distribution
                V_hat(t)= V_hat(t-1)+(I+epsilon)*dt; %integrate with euler the differential equation
              else
                  V_hat(t)=V_reset; %neuron fire a spike and we reset to -1
                   S(t)=1; %store the spike at the particular time step

                end
            end
            
    V_hat_trials(i,:)=V_hat; %store the V_hat for each trial 
    N_spikes=find(S)*dt; %find the number of spikes
    firing_rates(i)= length(N_spikes)/T(end); 
    isi(i)={diff(find(S)*dt)}; % store ISI for each trial in two arrays
    mean_isi(i)=mean(isi{i}); % calculate the mean ISI for each trial
    std_isi(i)=std(isi{i}); % calculate std ISI for each trial
    cv(i)=std_isi/mean_isi; %calculate the coefficient of variation for each trial
%     varian(i)=var(isi{i}); %uncomment to calculate the variance of ISI
%     fano(i)=varian(i)/mean_isi(i); %uncomment to calculate the Fano Factor as variance(ISI)/mean(ISI)
    end
   
    
% %Plot the membrane potentials of the two trials
% figure, plot(T,V_hat_trials(1,:)), hold on, plot(T,V_hat_trials(2,:))  
% ylim([-1.5 1.5]); %define x-axis limit
% yline(V_thresh,'r--'); %draw an horizontal line for the membrane threshold value
% yline(V_reset,'--c'); %draw an horizontal line for the membrane reset value
% legend('Trial 1','Trial 2', 'Vthresh','Vreset')
% title('IF Neuron with Noise')


%Plot the Interspike interval for the two trials in one plot - Question 2.6
figure()
subplot (2,1,1), histogram(isi{1}) %plot ISI for the first stimulus
xlim([0 1]);%define x-axis limit
ylabel('Frequency') %add y-axis label
title('First trial with Ie=4'); %add title
subplot (2,1,2), histogram(isi{2}) %plot ISI for the second stimulus stimulus
xlim([0 1]); %define x-axis limit
title('Second trial with Ie=4.2');%add title
xlabel('Interspike Interval'); %add x-axis label
ylabel('Frequency'); %add y-axis label