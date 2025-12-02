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
k = linspace(10^-5, 10^-4); % 1/s, effective mass transfer coefficient
ni0 = 1; % mol, entering amount of compound
ni0_tissue = 0.5; % mol, amount of compound in tissue

%% Analytical Solution
% Create a grid of W and k
[W, k] = meshgrid(W, k);
y = @(W, k)((ni0_tissue + (ni0 - ni0_tissue)*exp(-(k*W)/(rho_bed*v))).*(W<=Wc)) + ...
    (ni0_tissue + (ni0 - ni0_tissue).*exp((k*(W-2*Wc))/(rho_bed*v))).*(W > Wc);

n = y(W, k);
%n = (ni0_tissue + (ni0 - ni0_tissue)*exp(-(k*W)/(rho_bed*v))).*(W<=Wc) + ...
     %(ni0_tissue + (ni0 - ni0_tissue)*exp((k.*(W-2*Wc))/(rho_bed*v))).*(W>Wc);
%n = (ni0_tissue + (ni0 - ni0_tissue)*exp(-(k*W)/(rho_bed*v))).*(W<=Wc)
%n = (ni0_tissue + (ni0 - ni0_tissue)*exp((k.*(W-2*Wc))/(rho_bed*v))).*(W>Wc); 


surf(W, k, n);
% Add labels and title
xlabel('W (g)');
ylabel('k (s^{-1})');
zlabel('n (mol)');
shading('interp')

%% Publication Plot
f = figure(1);
surf(W, k, n);
grid on
set(gca, 'FontSize', 16, 'FontName', 'Times')
xlabel('$W (g)$', 'Interpreter', 'latex')
ylabel('$k_i (s^{-1})$', 'Interpreter', 'latex')
zlabel('$n_i (\frac{mol}{s})$', 'Interpreter', 'latex')
view(45, 30);
axis tight
shading interp
f.Position(3:4) = [1250/2 625];

% % Determine where surface intersects initial condition
% ndiff = n - ni0;
% C = contours(W, k, ndiff, [0 0]);
% mask = C(1,:) > 0; % Ignore initial condition
% C = C(:, mask); % Ignore initial condition
% WL = C(1, 2:end);
% kL = C(2, 2:end);
% nL = interp2(W, k, n, WL, kL);
% line(WL, kL, nL, 'Color', 'k', 'LineWidth', 3);

% Set size
f.Position(3:4) = [1250/2 625];

% Save
saveas(f, 'Model3.png')

