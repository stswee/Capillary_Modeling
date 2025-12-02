%% Capillary and density parameters
rho_cap = 1000; % # capillaries/mm2
d_cap = 5*10^-3; % mm
epsilon = rho_cap * (pi * (d_cap/2)^2); % unitless, void fraction
rho_cat = 1.06; % g/cm3, tissue density
rho_bed = (1 - epsilon) * rho_cat; % g/cm3, bed density

%% Fluid parameters
u = 0.114; % cm/s, velocity
d = 0.5; % cm, capillary bed diameter (biopsy is 3-8mm across)
A = pi*(d/2)^2; % cm2, cross-section area
v = u*A; % cm3/s, volumetric flow rate
mu = 4.5; % cP, viscosity (blood viscosity)
mu = mu * 0.0075/1000; % mmHg * s, viscosity (blood viscosity);
rho = 0.994; % g/cm3, fluid density (blood density)

%% Pressure parameters
k1 = 2.5*10^-3; % cm2, linear permeability coefficient (range from 10^-5 to 5*10^-5)
k2 = 5*10^-2; % g/(s^2*mmHg), quadratic permeability coefficient (range from 5*10^-2 to 10^-1)
alpha = ((mu/k1)*u + (rho/k2)*u^2)/(rho_bed * A); % mmHg/g, lumped parameter
P0 = 30; % mmHg, inlet pressure
Pc = 25; % mmHg, critical pressure
Wc = (P0 - Pc)/alpha; % g, critical weight

%% Reaction parameters
W = linspace(0, 10); % g, weight of sample
k = 5*10^-4; % 1/s, effective mass transfer coefficient
ni0 = 1; % mol, entering amount of compound
ni0_tissue = 0.5; % mol, amount of compound in tissue

%% Analytical Solution

% Model 1
n1 = ni0_tissue + (ni0 - ni0_tissue)*exp(-(k*W)/(rho_bed*v));

% Model 2
n2 = (ni0 + ni0_tissue)/2 + exp(-2.*k.*W/(rho_bed*v))*(ni0 - ni0_tissue)/2;

% Model 3
y = @(W, k)((ni0_tissue + (ni0 - ni0_tissue)*exp(-(k*W)/(rho_bed*v))).*(W<=Wc)) + ...
    (ni0_tissue + (ni0 - ni0_tissue).*exp((k*(W-2*Wc))/(rho_bed*v))).*(W > Wc);
n3 = y(W, k);

% Model 4
y = @(W, k)(((ni0 + ni0_tissue)/2 + exp(-2.*k*W/(rho_bed*v))*(ni0 - ni0_tissue)/2).*(W <= Wc) + ...
    ((ni0 + ni0_tissue)/2 + exp(2.*k*(W-2*Wc)/(rho_bed*v))*(ni0 - ni0_tissue)/2).*(W > Wc));
n4 = y(W, k);

% Model 5
% Pressure function
% Time component
A = 10; % mmHg
pulse = 60; % HR
w = 2*pi*60/pulse; % rad/s
B = 30; % mmHg
Ptime = @(t)(A*sin(w*t)+B); % Choose P(t)

% Weight function
% Sizing parameter
delta = Wc; % g
beta = @(w, delta)(exp(-w/delta)); % Choose beta

% Solving
Wspan = linspace(0, 10); % g, weight of sample
Tspan = linspace(0, 3*2*pi/w, 1000); % s, time, set for 3 periods
N = zeros(numel(Wspan), numel(Tspan)); % Preallocate
for i = 1:numel(Tspan)
    t = Tspan(i); % Choose value for t
    % Solve ODE
    [~, N(:,i)] = ode45(@(W,n) fun5(W, n, t, k, ni0_tissue, v, rho_bed, alpha, delta, Pc, Ptime, beta), Wspan, ni0);
    % fun(W, n, t, k, ni_tissue, v, rho_bed, alpha, delta, B, Pc, Ptime, beta)
end
n5 = mean(N, 2);

% Solving (t = 0)
Wspan = linspace(0, 10); % g, weight of sample
Tspan = 0; % s
N = zeros(numel(Wspan), numel(Tspan)); % Preallocate
for i = 1:numel(Tspan)
    t = Tspan(i); % Choose value for t
    % Solve ODE
    [~, N(:,i)] = ode45(@(W,n) fun5(W, n, t, k, ni0_tissue, v, rho_bed, alpha, delta, Pc, Ptime, beta), Wspan, ni0);
    % fun(W, n, t, k, ni_tissue, v, rho_bed, alpha, delta, B, Pc, Ptime, beta)
end
n50 = mean(N, 2);

% Model 6
% Solving
Wspan = linspace(0, 10); % g, weight of sample
Tspan = linspace(0, 3*2*pi/w, 1000); % s, time, set for 3 periods
Mi0 = [ni0; ni0_tissue];
Ni = zeros(numel(Wspan), numel(Tspan)); % Preallocate
Ni_tissue = zeros(numel(Wspan), numel(Tspan)); % Preallocate
for i = 1:numel(Tspan)
    t = Tspan(i); % Choose value for t
    % Solve ODE
    [~, M] = ode45(@(W, n) fun6(W, n, t, k, v, rho_bed, alpha, delta, Pc, Ptime, beta), Wspan, Mi0);
    Ni(:, i) = M(:, 1);
    Ni_tissue(:, i) = M(:, 2);
    % fun(W, n, t, k, ni_tissue, v, rho_bed, alpha, delta, B, Pc, Ptime, beta)
end
n6 = mean(Ni, 2); 

% Solving (t = 0)
Wspan = linspace(0, 10); % g, weight of sample
Tspan = 0; % s
Mi0 = [ni0; ni0_tissue];
Ni = zeros(numel(Wspan), numel(Tspan)); % Preallocate
Ni_tissue = zeros(numel(Wspan), numel(Tspan)); % Preallocate
for i = 1:numel(Tspan)
    t = Tspan(i); % Choose value for t
    % Solve ODE
    [~, M] = ode45(@(W, n) fun6(W, n, t, k, v, rho_bed, alpha, delta, Pc, Ptime, beta), Wspan, Mi0);
    Ni(:, i) = M(:, 1);
    Ni_tissue(:, i) = M(:, 2);
    % fun(W, n, t, k, ni_tissue, v, rho_bed, alpha, delta, B, Pc, Ptime, beta)
end
n60 = mean(Ni, 2); 
%% Plotting
plot(W, n1)
hold on
plot(W, n2)
plot(W, n3)
plot(W, n4)
plot(W, n5)
plot(W, n6)
plot(W, n50)
plot(W, n60)

%% Publication Plots
f = figure(1)
subplot(1, 2, 1);
data_line = plot(W, n1, ...
    W, n2, ...
    W, n3, ...
    W, n4, ...
    W, n5, ...
    W, n6, ...
    'LineWidth', 2)
set(data_line(2), 'LineStyle', '--');
set(data_line(3), 'LineStyle', ':');
set(data_line(4), 'LineStyle', '-.');
set(data_line(5), 'LineStyle', '--');
set(data_line(6), 'LineStyle', ':');
grid on
set(gca, 'FontSize', 20, 'FontName', 'Times')
xlabel('$W (g)$', 'Interpreter', 'latex')
ylabel('$n_i (\frac{mol}{s})$', 'Interpreter', 'latex')

legend({'$Model \; 1$', '$Model \; 2$', '$Model \; 3$', ...
    '$Model \; 4$', '$Model \; 5$', '$Model \; 6$'}, ...
    'Interpreter', 'latex', 'Location', 'Northwest')

subplot(1, 2, 2);
data_line = plot(W, n1, ...
    W, n2, ...
    W, n50, ...
    W, n60, ...
    'LineWidth', 2)
set(data_line(2), 'LineStyle', '--');
set(data_line(3), 'LineStyle', ':');
set(data_line(4), 'LineStyle', '-.');
grid on
set(gca, 'FontSize', 20, 'FontName', 'Times')
xlabel('$W (g)$', 'Interpreter', 'latex')
ylabel('$n_i (\frac{mol}{s})$', 'Interpreter', 'latex')
legend({'$Model \; 1$', '$Model \; 2$', '$Model \; 5 \; (t = 0)$', ...
    '$Model \; 6 \; (t = 0)$'}, ...
    'Interpreter', 'latex', 'Location', 'Northeast')
% Create textbox
annotation(f,'textbox',...
    [0.28 0 0.05 0.045],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FontName','Times');

% Create textbox
annotation(f,'textbox',...
    [0.72 0 0.05 0.045],...
    'String',{'(b)'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FontName','Times');

% Set size
f.Position(3:4) = [1250 835];

% Save
saveas(f, 'ModelAll.png')


%% Pressure function
% Time average pressure
function Ptimeavg = funPtimeavg(t, Ptime)
    % Special case where no time has passed
    if t <= 0
        Ptimeavg = 0; 
    else
        Ptimeavg = integral(Ptime, 0, t)/t;
    end
end

%% ODE (Model 5)
function dndW = fun5(W, n, t, k, ni_tissue, v, rho_bed, alpha, delta, Pc, Ptime, beta)
    % Evaluate pressure
    P = funPtimeavg(100, Ptime) - alpha * W - (funPtimeavg(t, Ptime) - Ptime(t)) * beta(W, delta); 

    % Determine rate based on pressure
    if P > Pc
        ri = -k * (n - ni_tissue)/v; % Rate law
    else
        ri = k * (n - ni_tissue)/v; % Starling law approximation
    end
    
    % Evaluate ODE
    dndW = ri/rho_bed;
end

%% ODE (Model 6)
function dndW = fun6(W, n, t, k, v, rho_bed, alpha, delta, Pc, Ptime, beta)
    ni = n(1); % Bloodstream
    ni_tissue = n(2); % Tissue

    % Evaluate pressure
    P = funPtimeavg(100, Ptime) - alpha * W - (funPtimeavg(t, Ptime) - Ptime(t)) * beta(W, delta); 

    % Determine rate based on pressure
    if P > Pc
        ri = -k * (ni - ni_tissue)/v; % Rate law
    else
        ri = k * (ni - ni_tissue)/v; % Starling law approximation
    end
    
    % Evaluate ODE
    dnidW = ri/rho_bed;
    dnitissuedW = -ri/rho_bed;

    dndW = [dnidW; dnitissuedW];
    
end