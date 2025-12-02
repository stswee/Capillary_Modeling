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
ni_tissue = 0.5; % mol, amount of compound in tissue

%% Pressure function
% Time component
A = 10; % mmHg
pulse = 60; % HR
w = 2*pi*60/pulse; % rad/s
B = 30; % mmHg
Ptime = @(t)(A*sin(w*t)+B); % Choose P(t)

%% Weight function
% Sizing parameter
delta = Wc; % g
beta = @(w, delta)(exp(-w/delta)); % Choose beta

%% Solving
Wspan = linspace(0, 10, 1000); % g, weight of sample
Tspan = linspace(0, 3*2*pi/w, 1000); % s, time, set for 3 periods
N = zeros(numel(Wspan), numel(Tspan)); % Preallocate
for i = 1:numel(Tspan)
    t = Tspan(i); % Choose value for t
    % Solve ODE
    [~, N(:,i)] = ode45(@(W,n) fun(W, n, t, k, ni_tissue, v, rho_bed, alpha, delta, Pc, Ptime, beta), Wspan, ni0);
    % fun(W, n, t, k, ni_tissue, v, rho_bed, alpha, delta, B, Pc, Ptime, beta)
end

%% Plot
surf(Wspan, Tspan, N')
% Add labels and title
xlabel('W (g)');
ylabel('t (s)');
zlabel('n (mol/s)');
shading('interp')

%% Publication Plot
f = figure(1); 
surf(Wspan, Tspan, N')
grid on
set(gca, 'FontSize', 16, 'FontName', 'Times')
xlabel('$W (g)$', 'Interpreter', 'latex')
ylabel('$t (s)$', 'Interpreter', 'latex')
zlabel('$n_i (\frac{mol}{s})$', 'Interpreter', 'latex')
view(45, 30);
axis tight
shading interp

% % Determine where surface intersects initial condition
% ndiff = N' - ni0;
% C = contours(Wspan, Tspan, ndiff, [0 0]);
% mask = C(1,:) > 0;
% C = C(:, mask); % Ignore initial condition
% WL = C(1, 2:end);
% kL = C(2, 2:end);
% nL = interp2(Wspan, Tspan, N', WL, kL, 'cubic');
% mask = WL < 10;
% WL = WL(mask);
% kL = kL(mask);
% nL = nL(mask);
% line(WL, kL, nL, 'Color', 'k', 'LineWidth', 3);

% Set size
f.Position(3:4) = [1250/2 625];

% Save
saveas(f, 'Model5.png')

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
function dndW = fun(W, n, t, k, ni_tissue, v, rho_bed, alpha, delta, Pc, Ptime, beta)
    % Evaluate pressure
    P = funPtimeavg(100, Ptime) - alpha * W - (funPtimeavg(100, Ptime) - Ptime(t)) * beta(W, delta); 

    % Determine rate based on pressure
    if P > Pc
        ri = -k * (n - ni_tissue)/v; % Rate law
    else
        ri = k * (n - ni_tissue)/v; % Starling law approximation
    end
    
    % Evaluate ODE
    dndW = ri/rho_bed;
end


