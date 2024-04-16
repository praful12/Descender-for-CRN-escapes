% syms t;
% x = sym('x' ,[1 2]);
% 
% f1 = x(2);
% f2 = -x(1);
% 
% f1_func = matlabFunction(f1, 'Vars', {t, x.'}); % Note the change in the 'Vars' argument
% f2_func = matlabFunction(f2, 'Vars', {t, x.'}); % Note the change in the 'Vars' argument
% 
% tspan = [0, 20]; % Replace with your desired time span
% initial_conditions = [0.5; 0.5]; % Use column vector instead of row vector
% 
% [t_sol, X] = ode45(@(t, x) [f1_func(t, x); f2_func(t, x)], tspan, initial_conditions);


syms t;
x = sym('x' ,[1 2]);

f = [x(2);-x(1)];

f_func = matlabFunction(f, 'Vars', {t, x.'}); % Note the change in the 'Vars' argument

tspan = [0, 20]; % Replace with your desired time span
initial_conditions = [0.5; 0.5]; % Use column vector instead of row vector

[t_sol, X] = ode45(@(t, x) f_func(t,x), tspan, initial_conditions);

% Assuming you have already obtained t_sol and X from ode45

% Extract the solution for x1 and x2
x1_sol = X(:, 1);
x2_sol = X(:, 2);

% Plot x1 and x2 as functions of time t
figure;
plot(t_sol, x1_sol, 'b', 'LineWidth', 2); % Blue line for x1
hold on;
plot(t_sol, x2_sol, 'r', 'LineWidth', 2); % Red line for x2
xlabel('Time t');
ylabel('Values');
legend('x1', 'x2');
title('Solution of the System of Differential Equations');
grid on;
hold off;

