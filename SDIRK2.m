%% Initialize constants

% Method specific constants
s = 2;
alpha = 1 - sqrt(2)/2;
A = [alpha, 0; 1 - alpha, alpha];
b = [1 - alpha; alpha];
c = [alpha; 1];

% Problem specific constants
t_max = 5;
steps = 1000;
Lambda = [-1, 10; 0, -3];
y_0 = [1; 1];

h = t_max/steps;
t = [(0:h:t_max);(0:h:t_max)];

%% Calculations
% Initialize variables
I = eye(2);
W = zeros(length(y_0), s);
S_i = zeros(length(y_0), 1);
S_n = zeros(length(y_0), 1);
y = zeros(length(y_0), steps);
y(:,1) = y_0;

% Start stepping
for n = 1:steps
    S_n = 0*S_n;
    for i = 1:s
        S_i = 0*S_i;
        for j = 1:i-1
            S_i = S_i + A(i,j)*Lambda*W(:,j);
        end
        W(:,i) = (I - h*A(i,i)*Lambda)^-1 * (y(:,n) + h*S_i);
        S_n = S_n + Lambda*W(:,i)*b(i);
    end
    y(:,n+1) = y(:,n) + h*S_n;
end

% Plot numerical solution
figure
plot(t', y')
title('Numerical solution for SDIRK 2')

% Plot exact solution
% y_exact = @(t) [1/29*(129*exp(-t) - 100*exp(-30*t)); exp(-30*t)];
y_exact = @(t) [6*exp(-t) - 5*exp(-3*t); exp(-3*t)];
figure
plot(t', y_exact(t(1,:))')
title('Exact solution for problem 1')

figure 
plot(t', abs(y'- y_exact(t(1,:))'))
title('Error for SDIRK 2')
        