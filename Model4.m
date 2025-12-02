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

% %% Analytical Solution (bloodstream)
% % Create a grid of W and k
% [W, K] = meshgrid(W, k);
% y = @(w, k)(((ni0+ni_tissue)/2 + (ni0 - ni_tissue)/2.*exp(-(k.*w)/(rho_bed*v))).*(w <= Wc) + ...
%     ((ni0+ni_tissue)/2 + (ni0 - ni_tissue)/2.*exp((k.*(w-2*Wc))/(rho_bed*v))).*(w > Wc));
% 
% n = y(W, K);
% figure(1)
% surf(W, K, n);
% % Add labels and title
% xlabel('W (g)');
% ylabel('k (s^{-1})');
% zlabel('n (mol)');
% shading('interp')
% 
% %% Analytical Solution (tissue)
% % Create a grid of W and k
% [W, K] = meshgrid(W, k);
% y = @(w, k)(((ni0+ni_tissue)/2 - (ni0 - ni_tissue)/2.*exp(-(k.*w)/(rho_bed*v))).*(w <= Wc) + ...
%     ((ni0+ni_tissue)/2 - (ni0 - ni_tissue)/2.*exp((k.*(w-2*Wc))/(rho_bed*v))).*(w > Wc));
% 
% n = y(W, K);
% figure(2)
% surf(W, K, n);
% % Add labels and title
% xlabel('W (g)');
% ylabel('k (s^{-1})');
% zlabel('n (mol)');
% shading('interp')

%% Publication Plot
f = figure(1);
tcl = tiledlayout(1, 2);

% Bloodstream plot
nexttile()
[W, K] = meshgrid(W, k);
y = @(W, k)(((ni0 + ni0_tissue)/2 + exp(-2.*k*W/(rho_bed*v))*(ni0 - ni0_tissue)/2).*(W <= Wc) + ...
    ((ni0 + ni0_tissue)/2 + exp(2.*k*(W-2*Wc)/(rho_bed*v))*(ni0 - ni0_tissue)/2).*(W > Wc));

n = y(W, K);
surf(W, K, n);
grid on
set(gca, 'FontSize', 16, 'FontName', 'Times')
xlabel('$W (g)$', 'Interpreter', 'latex')
ylabel('$k (s^{-1})$', 'Interpreter', 'latex')
zlabel('$n_i (\frac{mol}{s})$', 'Interpreter', 'latex')
view(45, 30);
axis tight
shading interp

% % Determine where surface intersects initial condition
% ndiff = n - ni0;
% C = contours(W, K, ndiff, [0 0]);
% mask = C(1,:) > 0; % Ignore initial condition
% C = C(:, mask); % Ignore initial condition
% WL = C(1, 2:end);
% kL = C(2, 2:end);
% nL = interp2(W, K, n, WL, kL);
% line(WL, kL, nL, 'Color', 'k', 'LineWidth', 3);

% Tissue plot
nexttile()
y = @(W, k)(((ni0 + ni0_tissue)/2 - exp(-2.*k*W/(rho_bed*v))*(ni0 - ni0_tissue)/2).*(W <= Wc) + ...
    ((ni0 + ni0_tissue)/2 - exp(2.*k*(W-2*Wc)/(rho_bed*v))*(ni0 - ni0_tissue)/2).*(W > Wc));

n = y(W, K);
surf(W, K, n);
grid on
set(gca, 'FontSize', 16, 'FontName', 'Times')
xlabel('$W (g)$', 'Interpreter', 'latex')
ylabel('$k_i (s^{-1})$', 'Interpreter', 'latex')
zlabel('$n_{i,tissue} (\frac{mol}{s})$', 'Interpreter', 'latex')
view(45, 30);
axis tight
shading interp

% % Determine where surface intersects initial condition
% ndiff = n - ni0_tissue;
% C = contours(W, K, ndiff, [0 0]);
% mask = C(1,:) > 0; % Ignore initial condition
% C = C(:, mask); % Ignore initial condition
% WL = C(1, 2:end);
% kL = C(2, 2:end);
% nL = interp2(W, K, n, WL, kL);
% line(WL, kL, nL, 'Color', 'k', 'LineWidth', 3);

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

% Save figure
saveas(f, 'Model4.png')
