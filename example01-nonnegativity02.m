% Example 01: Equilibrium State Of Ethanol Metabolism Model

% Step 01: Definition Of Problem Parameters

% Constant Problem Parameters

a = 0.08;
b = 0.20;
c = 0.05;
d = 0.01;

% Time Vector

T       = 100;
varphi  = @(x) exp(x) - 1;
h_tilde = 0.1;
h       = varphi(h_tilde);
t       = (0:h:T)';

% Initialize Solution Vectors

% Solution Vectors Of NSFDM

A_NSFDM = zeros(length(t),1);
B_NSFDM = zeros(length(t),1);
C_NSFDM = zeros(length(t),1);

% Solution Vectors Of EE

A_EE = zeros(length(t),1);
B_EE = zeros(length(t),1);
C_EE = zeros(length(t),1);

% Solution Vectors Of RK2

A_RK2 = zeros(length(t),1);
B_RK2 = zeros(length(t),1);
C_RK2 = zeros(length(t),1);

% Initial Conditions

A_NSFDM(1) = 0.5;
B_NSFDM(1) = 0;
C_NSFDM(1) = 0;

A_EE(1) = 0.5;
B_EE(1) = 0;
C_EE(1) = 0;

A_RK2(1) = 0.5;
B_RK2(1) = 0;
C_RK2(1) = 0;

% Step 02: Computations

A_NSFDM = A_NSFDM(1).*exp(-a*t);
A_EE    = A_EE(1).*exp(-a*t);
A_RK2   = A_RK2(1).*exp(-a*t);

% Solution Loop For NSFDM

for j = 1:1:(length(t)-1)
  %h = t(j+1) - t(j);
  C_NSFDM(j+1) = C_NSFDM(j)/(1 + (h*b)/(1 + h*b) + (h*c)/(d+C_NSFDM(j))) + (h*b*B_NSFDM(j) + h^2*a*b*A_NSFDM(j+1))/((1 + (h*b)/(1 + h*b) + (h*c)/(d+C_NSFDM(j))) * (1 + h*b));
  B_NSFDM(j+1) = (B_NSFDM(j) + h*a*A_NSFDM(j+1) + h*b*C_NSFDM(j+1))/(1 + h*b);
endfor

% Solution Loop For EE

for j = 1:1:(length(t)-1)
  B_EE(j+1) = B_EE(j) + h*(b*C_EE(j) - b*B_EE(j) + a*A_EE(1)*exp(-a*t(j+1)));
  C_EE(j+1) = C_EE(j) + h*(b*B_EE(j) - b*C_EE(j) - (c*C_EE(j))/(d+C_EE(j)));
endfor

% Solution Loop For RK2

for j = 1:1:(length(t)-1)
  k_1B = b*C_RK2(j) - b*B_RK2(j) + a*A_RK2(1)*exp(-a*t(j+1));
  k_1C = b*B_RK2(j) - b*C_RK2(j) - (c*C_RK2(j))/(d+C_RK2(j));
  k_2B = b*(C_RK2(j)+h*k_1C) - b*(B_RK2(j)+h*k_1B) + a*A_RK2(1)*exp(-a*t(j+1));
  k_2C = b*(B_RK2(j)+h*k_1B) - b*(C_RK2(j)+h*k_1C) - (c*(C_RK2(j)+h*k_1C))/(d+(C_RK2(j)+h*k_1C));
  B_RK2(j+1) = B_RK2(j) + (h/2)*(k_1B + k_2B);
  C_RK2(j+1) = C_RK2(j) + (h/2)*(k_1C + k_2C);
endfor

% Step 03: Plots Of Solution Components

##figure(1)
##plot(t,A_NSFDM)
##title('Amount Of Ethanol A(t) In Stomach')
##xlabel('t')
##ylabel('A(t)')
##
##figure(2)
##plot(t,B_NSFDM)
##title('Amount of Ethanol B(t) In Blood')
##xlabel('t')
##ylabel('B(t)')
##
##figure(3)
##plot(t,C_NSFDM)
##title('Amount Of Ethanol C(t) In Liver')
##xlabel('t')
##ylabel('C(t)')
##
##figure(4)
##plot(t,A_NSFDM+B_NSFDM+C_NSFDM)
##title('Total Amount Of Ethanol N(t) In Body')
##xlabel('t')
##ylabel('N(t)')

##figure(1)
##plot(t,A_RK2)
##title('Amount Of Ethanol A(t) In Stomach')
##xlabel('t')
##ylabel('A(t)')
##
##figure(2)
##plot(t,B_RK2)
##title('Amount of Ethanol B(t) In Blood')
##xlabel('t')
##ylabel('B(t)')
##
##figure(3)
##plot(t,C_RK2)
##title('Amount Of Ethanol C(t) In Liver')
##xlabel('t')
##ylabel('C(t)')
##
##figure(4)
##plot(t,A_RK2+B_RK2+C_RK2)
##title('Total Amount Of Ethanol N(t) In Body')
##xlabel('t')
##ylabel('N(t)')

figure(1)
title('Denominator function exp(h) - 1 for h = 0.1');

subplot(3,1,1);
plot(t(1:25:end),B_NSFDM(1:25:end),'linewidth',1.25);
title('Amount of ethanol in blood (NSFDM)');
xlabel('t');
ylabel('B(t)');
yticks([0 0.04 0.08 0.12]);

subplot(3,1,2);
plot(t(1:25:end),B_EE(1:25:end),'linewidth',1.25);
title('Amount of ethanol in blood (EE)');
xlabel('t');
ylabel('B(t)');
yticks([0 0.04 0.08 0.12]);

subplot(3,1,3);
plot(t(1:25:end),B_RK2(1:25:end),'linewidth',1.25);
title('Amount of ethanol in blood (RK2)');
xlabel('t');
ylabel('B(t)');
yticks([0 0.04 0.08 0.12]);
