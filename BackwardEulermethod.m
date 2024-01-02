clear all; close all;
%Creating a code to evaluate the value of an option using Backward Euler
%method

%%%%% VARIABLES USED %%%%%

S_min = 350;
S_max = 555;
K = 480; %strike price
T = 1; %time of expiry in months
r = 0.0525; %rate of interest (risk-free)
sigma = 0.1425; %volatility

%%%%% %%%%%
tic
%%% GRID MAKING %%%

Num_time = 5000; %time
dt = T/Num_time;
T_space = T:-dt:0;
T_points = (0:Num_time)*dt; %creating a vector with all the timesteps we consider

Num_S = 100; %Stock price grid points
ds = (S_max - S_min)/Num_S;
S_vec = S_min:ds:S_max; %creating a vector with all the timesteps we consider

i = 0:Num_S;

%%% %%%

p = zeros(Num_S+1, Num_time+1);

%put boundary
p(1,:) = K * exp(-r * (T_space));
p(:,Num_time+1) = max(K - S_vec, 0);
p(Num_S+1, :) = 0;

a = 0.5*(dt*r*i - dt*sigma^2*i.^2); %p(j-1) coeff
b = 1 + dt*r + dt*sigma^2*i.^2; %p(j) coeff
c = 0.5 * (-dt*sigma^2*i.^2 - dt*r*i); %p(j+1) coeff

%A = spdiags([a(2:end)' b(2:end)' c(2:end)'], -1:1, Num_S-1, Num_S-1);
A = diag(a(3:Num_S),-1) + diag(b(2:Num_S)) + diag(c(2:Num_S-1),1);
[L,U] = lu(A);

b = zeros(length(a)-2,1); %creating a vector for the boundary conditions

for j = Num_time:-1:1
    b(1) = a(2) * p(1,j); %we don't consider b(end) here as it is always 0
    p(2:Num_S, j) = A \ (p(2:Num_S,j+1) - b);
end
toc


% good exercise point P(T) > K - S(T) 
exercise_time = NaN;
exercise_strike = NaN;

for j = 2:Num_time
    exercise_point = find( p(j,1) >  (-S_vec(j) + K) );
   
    if ~isempty(exercise_point)
        exercise_time = Num_time*T_points(j);
        
        if exercise_point <= length(S_vec)
            exercise_strike = S_vec(exercise_time);
        else
            exercise_strike = NaN;
        end
        
        break; 
    end
end


figure(1)
%plot(S, p(:,1), '-', S, p(:,Num_time/2), '-', S, p(:,Num_time/4))
S_vec = S_vec(3:end);
p_1 = p(3:end,:);
plot(S_vec, p_1(:,1))
hold on
plot(exercise_strike, K - exercise_strike, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r')
xlabel('Strike price')
ylabel('Payoff of option')
legend('Option Payoff' , 'Exercise Strike')
title('Option Payoff')
text(exercise_strike, K - exercise_strike, sprintf(' %.2f', exercise_strike), 'HorizontalAlignment', 'right');
hold off

figure(2)
mesh(T_points, S_vec, p_1)
title('Solution for In the Money Put option using Backward Euler')
xlabel('time T')
ylabel('Strike price of stock')
zlabel('Payoff price P(t, S)')
colorbar


