close all; clear all

% Parameters
S_min = 0;
S_max = 600;
K = 520; %strike price
T = 1; %time of expiry in months
r = 0.0525; %rate of interest (risk-free)
sigma = 0.5313; %volatility

% Grid
Num_time = 500; % time
dt = T/Num_time;
T_grid = (0:Num_time)*dt;

Num_S = 100; % Stock price grid points
ds = (S_max - S_min)/Num_S;
i = 0:Num_S;
S = S_min + (i)*ds;

% Price value vector
p = zeros(Num_S+1, Num_time+1);

% Initial condition for a call option
p(:,1) = max(S - K, 0);

% Boundary conditions
p(1,:) = 0;
p(Num_S+1, :) = S_max - K*exp(-r*(T_grid));

% Discretization
for n = 1:Num_time
    for j = 2:Num_S
        p(j, n+1) = p(j, n) + dt * (r * j * p(j, n) + 0.5 * sigma^2 * j^2 * (p(j+1, n) - 2 * p(j, n) + p(j-1, n)));
    end
end

% Plotting
figure;
plot(S, p(:,1))
xlabel('Stock price')
ylabel('Payoff of option')
title('Solution for Call option at a specific time')