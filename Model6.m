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
%W = linspace(0, 10); % g, weight of sample
%k = linspace(0.01, 0.1, 10); % 1/s, effective mass transfer coefficient
k = 5*10^-4; % 1/s, effective mass transfer coefficient (range from 10^-5 to 10^-4)
ni0 = 1; % mol, entering amount of compound
ni_tissue0 = 0.5; % mol, amount of compound in tissue
Mi0 = [ni0; ni_tissue0];

%% Pressure function
% Time component
A = 10; % mmHg
pulse = 60; % HR
w = 2*pi*60/pulse; % rad/s
B = 30; % mmHg
Ptime = @(t)(A*sin(w*t)+B); % Choose P(t)
% Ptimeavg = funPtimeavg(t, Ptime); 

%% Weight function
% Sizing parameter
delta = Wc; % g
beta = @(w, delta)(exp(-w/delta)); % Choose beta

%% Solving
Wspan = linspace(0, 10); % g, weight of sample
Tspan = linspace(0, 3*2*pi/w, 1000); % s, time, set for 3 periods
Ni = zeros(numel(Wspan), numel(Tspan)); % Preallocate
Ni_tissue = zeros(numel(Wspan), numel(Tspan)); % Preallocate
for i = 1:numel(Tspan)
    t = Tspan(i); % Choose value for t
    % Solve ODE
    [~, M] = ode45(@(W, n) fun(W, n, t, k, v, rho_bed, alpha, delta, Pc, Ptime, beta), Wspan, Mi0);
    Ni(:, i) = M(:, 1);
    Ni_tissue(:, i) = M(:, 2);
    % fun(W, n, t, k, ni_tissue, v, rho_bed, alpha, delta, B, Pc, Ptime, beta)
end

%% Plot
% figure(1)
% surf(Wspan, Tspan, Ni')
% grid on
% set(gca, 'FontSize', 12, 'FontName', 'Times')
% xlabel('$W (g)$', 'Interpreter', 'latex')
% ylabel('$t (s)$', 'Interpreter', 'latex')
% zlabel('$n (\frac{mol}{s})$', 'Interpreter', 'latex')
% view(45, 30);
% axis tight
% shading interp

%% Plot
% figure(2)
% surf(Wspan, Tspan, Ni_tissue')
% % Add labels and title
% xlabel('W (g)');
% ylabel('t (s)');
% zlabel('n_i (mol/s)');
% shading('interp')

%% Publication Plot (about a minute to run)
f = figure(1);
tcl = tiledlayout(1,2);

% Bloodstream plot
nexttile()
surf(Wspan, Tspan, Ni')
grid on
set(gca, 'FontSize', 16, 'FontName', 'Times')
xlabel('$W (g)$', 'Interpreter', 'latex')
ylabel('$t (s)$', 'Interpreter', 'latex')
zlabel('$n_i (\frac{mol}{s})$', 'Interpreter', 'latex')
view(45, 30);
axis tight
shading interp

% Tissue plot
nexttile()
surf(Wspan, Tspan, Ni_tissue'u)
grid on
set(gca, 'FontSize', 16, 'FontName', 'Times')
xlabel('$W (g)$', 'Interpreter', 'latex')
ylabel('$t (s)$', 'Interpreter', 'latex')
zlabel('$n_{i,tissue} (\frac{mol}{s})$', 'Interpreter', 'latex')
view(45, 30);
axis tight
shading interp

% Set size
f.Position(3:4) = [1250 625];

% Create textbox
annotation(f,'textbox',...
    [0.25 0 0 0.05],'String',{'(a)'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FontName','Times');

% Create textbox
annotation(f,'textbox',...
    [0.725 0 0 0.05],'String',{'(b)'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',24,...
    'FontName','Times');

% Save
saveas(f, 'Model6.png')

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

%% ODE
function dndW = fun(W, n, t, k, v, rho_bed, alpha, delta, Pc, Ptime, beta)
    ni = n(1); % Bloodstream
    ni_tissue = n(2); % Tissue

    % Evaluate pressure
    P = funPtimeavg(100, Ptime) - alpha * W - (funPtimeavg(100, Ptime) - Ptime(t)) * beta(W, delta); 

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

