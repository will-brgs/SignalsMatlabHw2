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
figure();
NRange = [3,10,21];
for i = 1:3 % range of 3 due to 3 N values
N = NRange(i);
n = 0:1:(N-1);
r = ones(1,N);
x = conv(r,r);

% Subplot results
subplot(3,1,i)
stem(n,x(1:N),'filled') %Use Stem plot for DT
title(['Convolution for N = ',num2str(N)]);
xlabel('n');
ylabel('x[n]');
end
sgtitle('DT Convolution outputs for varying N')
%%
figure()
N = 21;
n = 0:N-1;
r = ones(1,N);
x = conv(r,r);

MRange = [10,25,50];
nRange = -200:1:200;
for i = 1:3 %range of 3 due to 3 M values
M = MRange(i);
sum = zeros(size(nRange));
    for k = -200:200
    sum = sum + dirac(n - k*M);
    end
y = conv(x,sum);

subplot(3,1,i)
stem(nRange,y,'filled')
title(['Convolution for N = 21 & M = ',num2str(M)]);
xlabel('n');
ylabel('y[n]');
sgtitle('Convolution of x[n] with Infinite Series Impulse Train for Different M Values');
end
