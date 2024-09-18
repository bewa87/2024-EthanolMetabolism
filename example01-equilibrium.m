% Example 01: Equilibrium State Of Ethanol Metabolism Model

% Step 01: Definition Of Problem Parameters

% Constant Problem Parameters

a = 0.08;
b = 0.20;
c = 0.05;
d = 0.01;

% Time Vector

T = 100;
h = 1;
t = (0:h:T)';

% Initialize Solution Vectors

A = zeros(length(t),1);
B = zeros(length(t),1);
C = zeros(length(t),1);

% Initial Conditions

A(1) = 0.5;
B(1) = 0;
C(1) = 0;

% Step 02: Computations

A = A(1).*exp(-a*t);

for j = 1:1:(length(t)-1)
  h = t(j+1) - t(j);
  C(j+1) = C(j)/(1 + (h*b)/(1 + h*b) + (h*c)/(d+C(j))) + (h*b*B(j) + h^2*a*b*A(j+1))/((1 + (h*b)/(1 + h*b) + (h*c)/(d+C(j))) * (1 + h*b));
  B(j+1) = (B(j) + h*a*A(j+1) + h*b*C(j+1))/(1 + h*b);
endfor

% Step 03: Plots Of Solution Components

figure(1)

subplot(2,2,1);
plot(t(1:2:end),A(1:2:end),'linewidth',1.25);
title('Amount of ethanol in stomach');
xlabel('t');
ylabel('A(t)');
yticks([0 0.2 0.4 0.6]);

subplot(2,2,2);
plot(t(1:2:end),B(1:2:end),'linewidth',1.25);
title('Amount of ethanol in blood');
xlabel('t');
ylabel('B(t)');
yticks([0 0.04 0.08 0.12]);

subplot(2,2,3);
plot(t(1:2:end),C(1:2:end),'linewidth',1.25);
title('Amount of ethanol in liver');
xlabel('t');
ylabel('C(t)');
yticks([0 0.002 0.004 0.006]);

subplot(2,2,4);
plot(t(1:2:end),A(1:2:end)+B(1:2:end)+C(1:2:end),'linewidth',1.25);
title('Total amount of ethanol in body');
xlabel('t');
ylabel('N(t)');
yticks([0 0.2 0.4 0.6]);
