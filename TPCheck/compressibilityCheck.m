% Check whether real gas effects need to be taken into account.
% They probably don't, since most of the time either temperature is high
% enough (e.g. initial gas) or pressure is low (post expansion).

% Dry air criticals
T_c = 132.63; % [K]
P_c_psi = 549.08; % [psi]
P_c = 37.858e5; % [Pa]

% Reduced pressure max and min
P_r_max = 1000/P_c_psi; % Max from reservoir
P_r_min = 95/P_c_psi; % Min Ch16, 001

% Reduced temperature max and min
T_min = 173; % Post expansion
T_max = 320; % Final state
T_r_min = T_min / T_c;
T_r_max = T_max / T_c;

% * At min temperature, maximum pressure, compressibility factor goes down
%   to about 0.7.
% * But min temp goes with min pressure etc., so this is typically not a
%   problem. Compressibility Z = 1 should really be fine for the bulk of
%   the gas.