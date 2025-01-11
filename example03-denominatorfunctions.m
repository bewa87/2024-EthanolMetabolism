% Example 01: Equilibrium State Of Ethanol Metabolism Model

% Step 01: Definition Of Problem Parameters

% Constant Problem Parameters

a = 0.08;
b = 0.20;
c = 0.05;
d = 0.01;

% Time Vector

varphi_1  = @(x) exp(x) - 1;
varphi_2  = @(x) 1 - exp(-x);
varphi_3  = @(x) x;
T         = 100;
h_1_tilde = 2;
h_1       = varphi_1(h_1_tilde);
h_2_tilde = 2;
h_2       = varphi_2(h_2_tilde);
h_3_tilde = 2;
h_3       = varphi_3(h_3_tilde);
t_1       = (0:h_1:T)';
t_2       = (0:h_2:T)';
t_3       = (0:h_3:T)';

% Initialize Solution Vectors

A_1 = zeros(length(t_1),1);
B_1 = zeros(length(t_1),1);
C_1 = zeros(length(t_1),1);
A_2 = zeros(length(t_2),1);
B_2 = zeros(length(t_2),1);
C_2 = zeros(length(t_2),1);
A_3 = zeros(length(t_3),1);
B_3 = zeros(length(t_3),1);
C_3 = zeros(length(t_3),1);

% Initial Conditions

A_1(1) = 0.5;
B_1(1) = 0;
C_1(1) = 0;
A_2(1) = 0.5;
B_2(1) = 0;
C_2(1) = 0;
A_3(1) = 0.5;
B_3(1) = 0;
C_3(1) = 0;

% Step 02: Computations

A_1 = A_1(1).*exp(-a*t_1);
A_2 = A_2(1).*exp(-a*t_2);
A_3 = A_3(1).*exp(-a*t_3);

for j = 1:1:(length(t_1)-1)
  h = t_1(j+1) - t_1(j);
  C_1(j+1) = C_1(j)/(1 + (h*b)/(1 + h*b) + (h*c)/(d+C_1(j))) + (h*b*B_1(j) + h^2*a*b*A_1(j+1))/((1 + (h*b)/(1 + h*b) + (h*c)/(d+C_1(j))) * (1 + h*b));
  B_1(j+1) = (B_1(j) + h*a*A_1(j+1) + h*b*C_1(j+1))/(1 + h*b);
endfor

for j = 1:1:(length(t_2)-1)
  h = t_2(j+1) - t_2(j);
  C_2(j+1) = C_2(j)/(1 + (h*b)/(1 + h*b) + (h*c)/(d+C_2(j))) + (h*b*B_2(j) + h^2*a*b*A_2(j+1))/((1 + (h*b)/(1 + h*b) + (h*c)/(d+C_2(j))) * (1 + h*b));
  B_2(j+1) = (B_2(j) + h*a*A_2(j+1) + h*b*C_2(j+1))/(1 + h*b);
endfor

for j = 1:1:(length(t_3)-1)
  h = t_3(j+1) - t_3(j);
  C_3(j+1) = C_3(j)/(1 + (h*b)/(1 + h*b) + (h*c)/(d+C_3(j))) + (h*b*B_3(j) + h^2*a*b*A_3(j+1))/((1 + (h*b)/(1 + h*b) + (h*c)/(d+C_3(j))) * (1 + h*b));
  B_3(j+1) = (B_3(j) + h*a*A_3(j+1) + h*b*C_3(j+1))/(1 + h*b);
endfor

% Step 03: Plots Of Solution Components

figure(1)

hold on

s1 = subplot(2,2,1);
plot(t_1(1:2:end),A_1(1:2:end),'linewidth',1.25);
set(s1,'title','Amount of ethanol in stomach');
xlabel('t');
ylabel('A(t)');
yticks([0 0.2 0.4 0.6]);

s2 = subplot(2,2,2);
plot(t_1(1:2:end),B_1(1:2:end),'linewidth',1.25);
set(s2,'title','Amount of ethanol in blood');
xlabel('t');
ylabel('B(t)');
yticks([0 0.04 0.08 0.12]);

s3 = subplot(2,2,3);
plot(t_1(1:2:end),C_1(1:2:end),'linewidth',1.25);
set(s3,'title','Amount of ethanol in liver');
xlabel('t');
ylabel('C(t)');
yticks([0 0.002 0.004 0.006]);

s4 = subplot(2,2,4);
plot(t_1(1:2:end),A_1(1:2:end)+B_1(1:2:end)+C_1(1:2:end),'linewidth',1.25);
set(s4,'title','Total amount of ethanol in body');
xlabel('t');
ylabel('N(t)');
yticks([0 0.2 0.4 0.6]);

S = axes('visible','off','title','Denominator function varphi(x) = exp(x) - 1');

hold off

figure(2)

hold on

s1 = subplot(2,2,1);
plot(t_2(1:2:end),A_2(1:2:end),'linewidth',1.25);
set(s1,'title','Amount of ethanol in stomach');
xlabel('t');
ylabel('A(t)');
yticks([0 0.2 0.4 0.6]);

s2 = subplot(2,2,2);
plot(t_2(1:2:end),B_2(1:2:end),'linewidth',1.25);
set(s2,'title','Amount of ethanol in blood');
xlabel('t');
ylabel('B(t)');
yticks([0 0.04 0.08 0.12]);

s3 = subplot(2,2,3);
plot(t_2(1:2:end),C_2(1:2:end),'linewidth',1.25);
set(s3,'title','Amount of ethanol in liver');
xlabel('t');
ylabel('C(t)');
yticks([0 0.002 0.004 0.006]);

s4 = subplot(2,2,4);
plot(t_2(1:2:end),A_2(1:2:end)+B_2(1:2:end)+C_2(1:2:end),'linewidth',1.25);
set(s4,'title','Total amount of ethanol in body');
xlabel('t');
ylabel('N(t)');
yticks([0 0.2 0.4 0.6]);

S = axes('visible','off','title','Denominator function varphi(x) = 1 - exp(-x)');

hold off

figure(3)

hold on

s1 = subplot(2,2,1);
plot(t_3(1:2:end),A_3(1:2:end),'linewidth',1.25);
set(s1,'title','Amount of ethanol in stomach');
xlabel('t');
ylabel('A(t)');
yticks([0 0.2 0.4 0.6]);

s2 = subplot(2,2,2);
plot(t_3(1:2:end),B_3(1:2:end),'linewidth',1.25);
set(s2,'title','Amount of ethanol in blood');
xlabel('t');
ylabel('B(t)');
yticks([0 0.04 0.08 0.12]);

s3 = subplot(2,2,3);
plot(t_3(1:2:end),C_3(1:2:end),'linewidth',1.25);
set(s3,'title','Amount of ethanol in liver');
xlabel('t');
ylabel('C(t)');
yticks([0 0.002 0.004 0.006]);

s4 = subplot(2,2,4);
plot(t_3(1:2:end),A_3(1:2:end)+B_3(1:2:end)+C_3(1:2:end),'linewidth',1.25);
set(s4,'title','Total amount of ethanol in body');
xlabel('t');
ylabel('N(t)');
yticks([0 0.2 0.4 0.6]);

S = axes('visible','off','title','Denominator function varphi(x) = x');

hold off
