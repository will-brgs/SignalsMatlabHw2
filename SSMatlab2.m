%% Signals and Systems Matlab Homework #2
%% Introduction
% * Author:                   Will Burgess
% * Class:                    ESE 351
% * Date:                     Created 1/25/2024, Last Edited 2/06/2024
% * With contributions from:  Mack Larosa, Tasha Igic, Mischa Tranor
%% Variable Initialization
close all 
R = 1e3 ; %Resistance in Ohms
C = 5e-6 ; %Capacitence in Farads
tau = R*C;
sampleFreq = 44.1e3; %Sampling freq in Hz
samplePeriod = 1/sampleFreq;

timeRange = 0:samplePeriod:15*tau;
%% Part 1: Compute discrete-time Convolution for a finite-length Function
%% Part 1 A
oldparam = sympref('HeavisideAtOrigin',1);
figure();
NRange = [3,10,21];
for i = 1:3 % range of 3 due to 3 N values
N = NRange(i);

r = ones(1,N+1);

x = conv(r,r);
val = (-1:1:length(x)-1);

% Subplot results
subplot(3,1,i)
stem(x,'filled') %Use Stem plot for DT
title(['Convolution for N = ', num2str(N)]);
xlabel('n');
ylabel('x[n]');
end
sgtitle('DT Convolution outputs for varying N')
%% Part 1 B
N = 21;
n = -200:1:200;
r = ones(1,N+1);
x = conv(r,r);
%x = x(1:length(n));

MRange = [10,25,50];

figure()
for i = 1:3 %range of 3 due to 3 M values
M = MRange(i);
sum = zeros(1,range(n)+2);
    for k = 1:(range(n)+2)
        if mod(k,M) == 0
            sum(k)=1;
        end   
    end
y = conv(x,sum);

subplot(3,1,i)
stem(linspace(-200,200, length(y)),y,'filled')
title(['Convolution for N = 21 & M = ',num2str(M)]);
xlabel('n');
ylabel('y[n]');
sgtitle('Convolution of x[n] with Infinite Series Impulse Train for Different M Values');
end

%% Part 2 A
R = 1e3; % Resistance in Ohms
C = 5e-6; % Capacitence in Farads
tau = R*C;
sampleFreq = 44.1e3; % Sampling frequency in Hz
samplePeriod = 1/sampleFreq;

LCCDE = zeros(length(sampleTimes),1);
LCCDE(1) = 0;

for n = 2:length(sampleTimes)
LCCDE(n) = inputFunction(n) - inputFunction(n-1) - (LCCDE(n-1) * ((samplePeriod/tau) -1));
end
