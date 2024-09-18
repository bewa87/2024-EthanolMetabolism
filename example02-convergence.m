% Example 02: Convergence

% Step 01: Definition Of Problem Parameters

% Constant Problem Parameters

a = 0.08;
b = 0.20;
c = 0.05;
d = 0.01;

% Initial Conditions

A_zero = 0.5;
B_zero = 0;
C_zero = 0;

% Time Vectors

T     = 100;

h_RK4   = 0.001;
t_RK4   = (0:h_RK4:T)';

h_NSFDM = [4 2 1 0.5 0.25 0.125];
t_h1    = (0:h_NSFDM(1):T)';
t_h2    = (0:h_NSFDM(2):T)';
t_h3    = (0:h_NSFDM(3):T)';
t_h4    = (0:h_NSFDM(4):T)';
t_h5    = (0:h_NSFDM(5):T)';
t_h6    = (0:h_NSFDM(6):T)';

% Initialize Solution Vectors

% Solution Vectors Of RK4

A_RK4 = zeros(length(t_RK4),1);
B_RK4 = zeros(length(t_RK4),1);
C_RK4 = zeros(length(t_RK4),1);
A_RK4 = A_zero.*exp(-a*t_RK4);

% Solution Vectors Of NSFDM1

A_NSFDM1 = zeros(length(t_h1),1);
B_NSFDM1 = zeros(length(t_h1),1);
C_NSFDM1 = zeros(length(t_h1),1);
A_NSFDM1 = A_zero.*exp(-a*t_h1);

% Solution Vectors Of NSFDM2

A_NSFDM2 = zeros(length(t_h2),1);
B_NSFDM2 = zeros(length(t_h2),1);
C_NSFDM2 = zeros(length(t_h2),1);
A_NSFDM2 = A_zero.*exp(-a*t_h2);

% Solution Vectors Of NSFDM3

A_NSFDM3 = zeros(length(t_h3),1);
B_NSFDM3 = zeros(length(t_h3),1);
C_NSFDM3 = zeros(length(t_h3),1);
A_NSFDM3 = A_zero.*exp(-a*t_h3);

% Solution Vectors Of NSFDM4

A_NSFDM4 = zeros(length(t_h4),1);
B_NSFDM4 = zeros(length(t_h4),1);
C_NSFDM4 = zeros(length(t_h4),1);
A_NSFDM4 = A_zero.*exp(-a*t_h4);

% Solution Vectors Of NSFDM5

A_NSFDM5 = zeros(length(t_h5),1);
B_NSFDM5 = zeros(length(t_h5),1);
C_NSFDM5 = zeros(length(t_h5),1);
A_NSFDM5 = A_zero.*exp(-a*t_h5);

% Solution Vectors Of NSFDM6

A_NSFDM6 = zeros(length(t_h6),1);
B_NSFDM6 = zeros(length(t_h6),1);
C_NSFDM6 = zeros(length(t_h6),1);
A_NSFDM6 = A_zero.*exp(-a*t_h6);

% Step 02: Computations

% Loop For RK4

for j = 1:1:(length(t_RK4)-1)
  k_1B   = b*C_RK4(j) - b*B_RK4(j) + a*A_RK4(1)*exp(-a*t_RK4(j+1));
  k_1C   = b*B_RK4(j) - b*C_RK4(j) - (c*C_RK4(j))/(d+C_RK4(j));
  k_2B   = b*(C_RK4(j)+0.5*h_RK4*k_1C) - b*(B_RK4(j)+0.5*h_RK4*k_1B) + a*A_RK4(1)*exp(-a*t_RK4(j+1));
  k_2C   = b*(B_RK4(j)+0.5*h_RK4*k_1B) - b*(C_RK4(j)+0.5*h_RK4*k_1C) - (c*(C_RK4(j)+0.5*h_RK4*k_1C))/(d+(C_RK4(j)+0.5*h_RK4*k_1C));
  k_3B   = b*(C_RK4(j)+0.5*h_RK4*k_2C) - b*(B_RK4(j)+0.5*h_RK4*k_2B) + a*A_RK4(1)*exp(-a*t_RK4(j+1));
  k_3C   = b*(B_RK4(j)+0.5*h_RK4*k_2B) - b*(C_RK4(j)+0.5*h_RK4*k_2C) - (c*(C_RK4(j)+0.5*h_RK4*k_2C))/(d+(C_RK4(j)+0.5*h_RK4*k_2C));
  k_4B   = b*(C_RK4(j)+h_RK4*k_3C) - b*(B_RK4(j)+h_RK4*k_3B) + a*A_RK4(1)*exp(-a*t_RK4(j+1));
  k_4C   = b*(B_RK4(j)+h_RK4*k_3B) - b*(C_RK4(j)+h_RK4*k_3C) - (c*(C_RK4(j)+h_RK4*k_3C))/(d+(C_RK4(j)+h_RK4*k_3C));
  B_RK4(j+1) = B_RK4(j) + (h_RK4/6)*(k_1B+2*k_2B+2*k_3B+k_4B);
  C_RK4(j+1) = C_RK4(j) + (h_RK4/6)*(k_1C+2*k_2C+2*k_3C+k_4C);
endfor

% Loop For NSFDM1

for j = 1:1:(length(t_h1)-1)
  C_NSFDM1(j+1) = C_NSFDM1(j)/(1 + (h_NSFDM(1)*b)/(1 + h_NSFDM(1)*b) + (h_NSFDM(1)*c)/(d+C_NSFDM1(j))) + (h_NSFDM(1)*b*B_NSFDM1(j) + (h_NSFDM(1))^2*a*b*A_NSFDM1(j+1))/((1 + (h_NSFDM(1)*b)/(1 + h_NSFDM(1)*b) + (h_NSFDM(1)*c)/(d+C_NSFDM1(j))) * (1 + h_NSFDM(1)*b));
  B_NSFDM1(j+1) = (B_NSFDM1(j) + h_NSFDM(1)*a*A_NSFDM1(j+1) + h_NSFDM(1)*b*C_NSFDM1(j+1))/(1 + h_NSFDM(1)*b);
endfor

% Loop For NSFDM2

for j = 1:1:(length(t_h2)-1)
  C_NSFDM2(j+1) = C_NSFDM2(j)/(1 + (h_NSFDM(2)*b)/(1 + h_NSFDM(2)*b) + (h_NSFDM(2)*c)/(d+C_NSFDM2(j))) + (h_NSFDM(2)*b*B_NSFDM2(j) + (h_NSFDM(2))^2*a*b*A_NSFDM2(j+1))/((1 + (h_NSFDM(2)*b)/(1 + h_NSFDM(2)*b) + (h_NSFDM(2)*c)/(d+C_NSFDM2(j))) * (1 + h_NSFDM(2)*b));
  B_NSFDM2(j+1) = (B_NSFDM2(j) + h_NSFDM(2)*a*A_NSFDM2(j+1) + h_NSFDM(2)*b*C_NSFDM2(j+1))/(1 + h_NSFDM(2)*b);
endfor

% Loop For NSFDM3

for j = 1:1:(length(t_h3)-1)
  C_NSFDM3(j+1) = C_NSFDM3(j)/(1 + (h_NSFDM(3)*b)/(1 + h_NSFDM(3)*b) + (h_NSFDM(3)*c)/(d+C_NSFDM3(j))) + (h_NSFDM(3)*b*B_NSFDM3(j) + (h_NSFDM(3))^2*a*b*A_NSFDM3(j+1))/((1 + (h_NSFDM(3)*b)/(1 + h_NSFDM(3)*b) + (h_NSFDM(3)*c)/(d+C_NSFDM3(j))) * (1 + h_NSFDM(3)*b));
  B_NSFDM3(j+1) = (B_NSFDM3(j) + h_NSFDM(3)*a*A_NSFDM3(j+1) + h_NSFDM(3)*b*C_NSFDM3(j+1))/(1 + h_NSFDM(3)*b);
endfor

% Loop For NSFDM4

for j = 1:1:(length(t_h4)-1)
  C_NSFDM4(j+1) = C_NSFDM4(j)/(1 + (h_NSFDM(4)*b)/(1 + h_NSFDM(4)*b) + (h_NSFDM(4)*c)/(d+C_NSFDM4(j))) + (h_NSFDM(4)*b*B_NSFDM4(j) + (h_NSFDM(4))^2*a*b*A_NSFDM4(j+1))/((1 + (h_NSFDM(4)*b)/(1 + h_NSFDM(4)*b) + (h_NSFDM(4)*c)/(d+C_NSFDM4(j))) * (1 + h_NSFDM(4)*b));
  B_NSFDM4(j+1) = (B_NSFDM4(j) + h_NSFDM(4)*a*A_NSFDM4(j+1) + h_NSFDM(4)*b*C_NSFDM4(j+1))/(1 + h_NSFDM(4)*b);
endfor

% Loop For NSFDM5

for j = 1:1:(length(t_h5)-1)
  C_NSFDM5(j+1) = C_NSFDM5(j)/(1 + (h_NSFDM(5)*b)/(1 + h_NSFDM(5)*b) + (h_NSFDM(5)*c)/(d+C_NSFDM5(j))) + (h_NSFDM(5)*b*B_NSFDM5(j) + (h_NSFDM(5))^2*a*b*A_NSFDM5(j+1))/((1 + (h_NSFDM(5)*b)/(1 + h_NSFDM(5)*b) + (h_NSFDM(5)*c)/(d+C_NSFDM5(j))) * (1 + h_NSFDM(5)*b));
  B_NSFDM5(j+1) = (B_NSFDM5(j) + h_NSFDM(5)*a*A_NSFDM5(j+1) + h_NSFDM(5)*b*C_NSFDM5(j+1))/(1 + h_NSFDM(5)*b);
endfor

% Loop For NSFDM6

for j = 1:1:(length(t_h6)-1)
  C_NSFDM6(j+1) = C_NSFDM6(j)/(1 + (h_NSFDM(6)*b)/(1 + h_NSFDM(6)*b) + (h_NSFDM(6)*c)/(d+C_NSFDM6(j))) + (h_NSFDM(6)*b*B_NSFDM6(j) + (h_NSFDM(6))^2*a*b*A_NSFDM6(j+1))/((1 + (h_NSFDM(6)*b)/(1 + h_NSFDM(6)*b) + (h_NSFDM(6)*c)/(d+C_NSFDM6(j))) * (1 + h_NSFDM(6)*b));
  B_NSFDM6(j+1) = (B_NSFDM6(j) + h_NSFDM(6)*a*A_NSFDM6(j+1) + h_NSFDM(6)*b*C_NSFDM6(j+1))/(1 + h_NSFDM(6)*b);
endfor

% Step 03: Error Calculation

err1_B = max(abs(B_RK4(1:4000:end)-B_NSFDM1));
err1_C = max(abs(C_RK4(1:4000:end)-C_NSFDM1));
err1   = max(err1_B,err1_C);

err2_B = max(abs(B_RK4(1:2000:end)-B_NSFDM2));
err2_C = max(abs(C_RK4(1:2000:end)-C_NSFDM2));
err2   = max(err2_B,err2_C);

err3_B = max(abs(B_RK4(1:1000:end)-B_NSFDM3));
err3_C = max(abs(C_RK4(1:1000:end)-C_NSFDM3));
err3   = max(err3_B,err3_C);

err4_B = max(abs(B_RK4(1:500:end)-B_NSFDM4));
err4_C = max(abs(C_RK4(1:500:end)-C_NSFDM4));
err4   = max(err4_B,err4_C);

err5_B = max(abs(B_RK4(1:250:end)-B_NSFDM5));
err5_C = max(abs(C_RK4(1:250:end)-C_NSFDM5));
err5   = max(err5_B,err5_C);

err6_B = max(abs(B_RK4(1:125:end)-B_NSFDM6));
err6_C = max(abs(C_RK4(1:125:end)-C_NSFDM6));
err6   = max(err6_B,err6_C);

err    = [err1 err2 err3 err4 err5 err6];

% Step 04: Error Plot

figure(1)
plot(h_NSFDM,err,'marker','+','linewidth',1.25);
title('Convergence plot');
xlabel('h');
ylabel('error');
yticks([0 0.01 0.02 0.03]);
