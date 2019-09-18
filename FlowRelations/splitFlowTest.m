% Testing out split flow
% This didn't get anywhere, since the pressure drop for compressible flow
% doesn't tell us how the other thermodynamic states vary--rho is certainly
% not constant (if it were, then we would know T change from pressure drop)

%% Parameters and precompute
Q = 287;
cp = 1000;
cv = cp - Q;
gamma = cp/cv;

machAreaFunction = precomputeMachAreaFunction(gamma);

%% Input
T0 = 300;
p0 = 1000 * 6895; % [psi -> Pa]
M = 0.1;

% Compute state
T(1) = temperatureMachFunction(gamma, M) * T0;
p(1) = pressureMachFunction(gamma, M) * p0;
rho(1) = p / Q / T;
s(1) = cv * log(p(1) / 1) + cp * log(1 / rho(1));
M(1) = M;

%% Section 1: area change
% Area ratio: contraction
Aratio = areaMachFunction(gamma, M);
Aratio = Aratio/2;
M = machAreaFunction(Aratio);

M(2) = M;
T(2) = temperatureMachFunction(gamma, M) * T0;
p(2) = pressureMachFunction(gamma, M) * p0;
rho(2) = p / Q / T;
s(2) = cv * log(p(2) / 1) + cp * log(1 / rho(2));

%% Section 2: pressure loss (head loss)
% Note: due to sudden contraction
% semi-empirical relation; see also
% http://thermopedia.com/content/659/
% Model as incompressible

% Compute u from flow
u = sqrt( 2 * cv * (T0 - T));
rho = p / (Q * T);

coeffLoss = 0.45 + 0.75; % Ballpark on loss coefficients
pDrop = 0.5 * coeffLoss * rho * u^2;

p(3) = p - pDrop;
p0 = p0 - pDrop;

%% Section 3: area change 2
% Area ratio: contraction
Aratio = areaMachFunction(gamma, M);
Aratio = Aratio/2;

M = machAreaFunction(Aratio);
T = temperatureMachFunction(gamma, M) * T0;
p = pressureMachFunction(gamma, M) * p0;

%% Section 4: 