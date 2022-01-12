% Code Tested on Matlab R2021a 
% Assignment 2 in Cortical Modelling


%% Define the parameters of this simulation

clear all       

dt = 10^-4;    %define the timestep in seconds
max_T= 10 ;    %define the maximum time of simulation
T= 0:dt:max_T; %time vector
alpha = 4;     %define the inhibitory coupling strenght alpha, try here different values for question 2.2


%Initialise the necessary variables that will change with time
r_1= zeros(size(T)); %initialise the r_1 variable to store the rate values during the simulation
r_2= zeros(size(T)); %initialise the r_2 variable
r_1(1) = 1;    %initial rate of neuron 1 derived in the pdf file
r_2(1) = 4;    %initial rate of neuron 2 derived in the pdf file


for t=2:length(T) %start from index 2 to compute the next step as a function of the previous step
r_1(t)= r_1(t-1)+dt*(r_2(t-1)-sqrt(alpha)); %implementation of the approximate solution for r_1 derived in the pdf file
r_2(t)= r_2(t-1)+dt*(-alpha*(r_1(t-1))+alpha); %implementation of the approximate solution for r_1 derived in the pdf file
end

%Implement the exact solutions of the differential equations system
sol_1=sin(sqrt(alpha)*T)+1; %exact solution: r_1(t)=sin(√αt)+1
sol_2=sqrt(alpha)*cos(sqrt(alpha)*T)+sqrt(alpha); %exact solution: r_2(t)=√αcos(√αt)+√α

%Plot the rate of the two neurons
figure, plot(T, r_1, T, r_2) %plot the approximate solutions as a function of time
hold on, plot(T,sol_1, T, sol_2) %plot the exact solutions on the same graph
legend('ApproximateRate1','ApproximateRate2','SolutionRate1','SolutionRate2') %add a legend
ylim([0 5]) %define the y-axis limit
xlabel('Time')
ylabel('Rate')

%Phase Plot Portrait
figure, plot(r_1, r_2)
title('Phase portrait');
xlabel('Rate Neuron 1');
ylabel('Rate Neuron 2');

%Find the Period
[idx,idx]=findpeaks(r_1); % detect the peaks in the signal r_1 or r_2
P=T(idx(2))-T(idx(1)); %use the spacing between peaks to determine periodicity
freq_=1/P; %calculate the frequency of the wave
disp(['The sine wave has a time period of ' num2str(P),[' seconds'], ' and a frequency of '  num2str(freq_),[' Hz']]);


 