function [physConst, t0, icAirgun, icBubble] = configAirgun(str, ...
    airgunPressure,airgunLength,airgunPortArea,airgunDepth, ...
    airgunCrossSectionalArea, airgunFiringChamberProfile, ...
    airgunOperatingChamberProfile, bubbleInitialVolume, ...
    shuttleBdryPenaltyStrength)

% Set defaults for last arguments backward compatibility
if nargin == 5
    if strcmpi(str, 'GeneralAirgun')
        error(['Invalid number of input arguments for' ...
               'GeneralAirgun setting.']);
    end
elseif nargin < 5 || nargin ~= 10
    error('Invalid number of input arguments.')
end

switch str
    case 'Bolt1500LL'
        % 'realistic' initial conditions that will be used in paper.        
        % bubble volume is set equal to the airgun volume 
        % temperature and pressure in the bubble is the same as the ambient properties.
        % The airgun is set to have constant values.
        
        t0 = 0;
        
        p0a_imperial = airgunPressure; % air gun pressure [psi]
        p0a = p0a_imperial * 6894.8; % air gun pressure [Pa]
        physConst.p0a = p0a;
        
        physConst.L = airgunLength; % Length of Airgun in meters

        %cross_sectional_area = 24; % [in^2] cross sectional area is calculated by assuming 1m long (39.3701) for 600in3 airgun.
        cross_sectional_area = airgunPortArea;
        airgunLengthImperial = physConst.L/0.0254; % length in inches
        V_imperial = airgunLengthImperial * cross_sectional_area;
        V = V_imperial * 1.63871e-5; % air gun volume [m^3]
        
        A_imperial = airgunPortArea; % air gun port area [in^2]
        A = A_imperial * 6.4516e-4; % air gun port area [m^2]
        physConst.A = A;
        
        
       
        physConst.rho_inf = 1e3; % density [kg/m^3]
        physConst.pa = 1e5; % atmospheric pressure [Pa]
        depth = airgunDepth; % depth [m]
        g = 9.8; % gravitational acceleration [m/s^2]
        physConst.p_inf = physConst.pa + physConst.rho_inf*g*depth; % ambient pressure at depth [Pa]
        
        physConst.c_v = 718; % heat capacity of air at constant volume [J/kgK]
        physConst.c_inf = 1482; % speed of sound in water [m/s]
        physConst.Q = 287.06; % specific gas constant for dry air [J/kgK]
        physConst.R_G = physConst.Q; % DUPLICATE
        physConst.Tinf = 288; % temperature assumed constant throughout the system [K]
        physConst.gamma = 1.4; % ratio of heat capacities for dry air
        
        physConst.AirgunCutoffTime = 0.010; % time when air gun stops firing [s]
        
                
        % Air gun
        p = physConst.p0a;
        T = physConst.Tinf;
        Q = physConst.Q;
        c_v = physConst.c_v;
        
        rho = p/(Q*T);
        e = c_v*rho*T;
        
        icAirgun.rho0 = @(x)0*x + rho;
        icAirgun.rv0  = @(x)0*x;
        icAirgun.e0   = @(x)0*x + e;
        icAirgun.p0   = @(x)0*x + p;
        
        % Bubble
        p = physConst.p_inf;
        T = physConst.Tinf;
        Q = physConst.Q;
        c_v = physConst.c_v;
        
        %V_airgun_const = 250 * 1.63871e-5;
        V_airgun_const = 600 * 1.63871e-5;
        icBubble.R = (3/(4*pi) * V_airgun_const)^(1/3);
        %icBubble.R = (3/(4*pi) * V)^(1/3);
        icBubble.Rdot = 0;
        Tvary = T;
        icBubble.m = p*V_airgun_const / (Q*Tvary); 
        icBubble.E = c_v * icBubble.m * Tvary;
        %icBubble.m = p*V_airgun_const / (Q*T); 
        %icBubble.E = c_v * icBubble.m * T;
    
    case 'GeneralAirgun' % Based on the TPS
        % Modified version to accommodate additional parameters relevant to
        % the airgun
        
        % Initial time [s]
        t0 = 0;
        % Airgun pressure conversion
        p0a_imperial = airgunPressure; % Air gun pressure [psi]
        p0a = p0a_imperial * 6894.8;   % Air gun pressure [Pa]
        physConst.p0a = p0a;
        % Airgun length [m]
        physConst.L = airgunLength;

%         cross_sectional_area = 24; % [in^2] cross sectional area is calculated by assuming 1m long (39.3701) for 600in3 airgun.
% %         cross_sectional_area = airgunPortArea;
%         airgunLengthImperial = physConst.L/0.0254; % length in inches
%         V_imperial = airgunLengthImperial * cross_sectional_area;
%         V = V_imperial * 1.63871e-5; % air gun volume [m^3]

        % Cross-sectional area [in^2]
        crossSectionalAreaImperial = airgunCrossSectionalArea; % 12.5 in^2 in [Watson 2019]
        % Compute volume
        airgunLengthImperial = airgunLength/0.0254; % [in]
        V_imperial = crossSectionalAreaImperial * airgunLengthImperial; % [cui]
        physConst.airgunVolume = V_imperial * 1.63871e-5; % Air gun volume [m^3]
        
        %% Water constants
        physConst.rho_inf = 1e3; % density [kg/m^3]
        physConst.pa = 1e5;      % atmospheric pressure [Pa]
        physConst.Tinf = 288;    % temperature (constant) throughout [K]
        depth = airgunDepth;     % depth [m]
        g = 9.8;                 % gravitational acceleration [m/s^2]
        % Compute ambient pressure at depth [Pa]
        physConst.p_inf = physConst.pa + physConst.rho_inf*g*depth;
        
        %% (Air calorically perfect gas) constants
        physConst.c_v = 718;    % heat capacity of air at constant volume [J/kgK]
        physConst.c_inf = 1482; % speed of sound in water [m/s]
        physConst.Q = 287.06;   % specific gas constant for dry air [J/kgK]
        physConst.gamma = 1.4;  % ratio of heat capacities for dry air
        
        %% Fixed open-time model: unused
        % Time when air gun stops firing [s]
        % Warning: suggested to be 0.010 s in paper, but specified as 0.04
        % in airgun1d/HEAD
        % physConst.AirgunCutoffTime = 0.04;
        
        %% Compute airgun initial conditions
        p = physConst.p0a;
        T = physConst.Tinf;
        Q = physConst.Q;
        c_v = physConst.c_v;
        rho = p/(Q*T);
        e = rho*c_v*T;
        % Set uniform, quiescent airgun initial conditions
        icAirgun.rho0 = @(x)0*x + rho;
        icAirgun.rv0  = @(x)0*x;
        icAirgun.e0   = @(x)0*x + e;
        icAirgun.p0   = @(x)0*x + p;
        
        %% Set bubble initial conditions and gas air parameters
        p = physConst.p_inf;
        T = physConst.Tinf;
        Q = physConst.Q;
        c_v = physConst.c_v;
        
        % Convert to metric the initial bubble volume argument
        % Note: LeightonWatson/airgun1d/HEAD uses a fixed presribed volume
        % 600, off by just a little from the chamber volume [cui]
        bubbleInitialVolumeMetric = bubbleInitialVolume * 1.63871e-5; 
        icBubble.R = (3/(4*pi) * bubbleInitialVolumeMetric)^(1/3);
        icBubble.Rdot = 0;
        icBubble.m = p*bubbleInitialVolumeMetric / (Q*T); 
        icBubble.E = icBubble.m * c_v * T;
        
        % Convert cross sectional area from [in^2] to [m^2]
        physConst.crossSectionalArea = airgunCrossSectionalArea*0.00064516;
        % Set left pressure = right pressure
        physConst.p_R0 = p0a; % [Pa]
        
        % Set fixed shuttle assembly mass [kg]
        physConst.shuttleAssemblyMass = 63 * .454;
        
        % Convert port area from [in^2] to [m^2]
        physConst.APortTotal = airgunPortArea*0.00064516;
        % Ludvig uses: 0.0464 [m^2] ~ 72 in^2
        
        % Set fixed operating chamber length [m]
        physConst.operatingChamberLength = 3.009 *0.0254;
        
        % Set firing chamber and operating chamber profiles A(x)
        physConst.airgunFiringChamberProfile = ...
            airgunFiringChamberProfile;
        physConst.airgunOperatingChamberProfile = ...
            airgunOperatingChamberProfile;
        
        % Ludvig paramters for airgun area:
%         A_L = 0.0586; % [m^2]
%         A_R = 0.0506; % [m^2]
        % Artificial shuttle left and right areas
        physConst.shuttle_area_left = pi/4 * (11.2 * 0.0254)^2; % [m^2]
        % physConst.shuttle_area_left = physConst.crossSectionalArea; % [m^2]
        
        % Reasonable guess to shuttle right area: ratio using Ludvig's
        % code, although the areas are different [deprecate]
%         physConst.shuttle_area_right = 0.0506/0.0586* ...
%             physConst.crossSectionalArea; % [m^2]

        % Front of shuttle: projected area of cushioned side of piston
        physConst.shuttle_area_right = pi/4 * ( ...
            (11.1 * 0.0254)^2); % [m^2]        
        % Rear of shuttle: projected area of piston, minus the shaft area
        physConst.shuttle_area_right_rear = pi/4 * ( ...
            (11.1 * 0.0254)^2 - (2.1 * 0.0254)^2); % [m^2]
        
        % Lead-in length where shuttle can move without exposing the air
        physConst.portLead = 0.35 * 0.0254; % [m]
        
        physConst.flangeDepth = 3 * 0.0254; % [m]
        % Approximate the flage ID to be equal to the chamber
        physConst.plugVolumeFn = @(xi) ...
            physConst.crossSectionalArea * (xi + physConst.flangeDepth);

        % Set shuttle penalty parameter
        % Note: ~1e11 makes sense based on linear elasticity
        physConst.shuttleBdryPenaltyStrength = shuttleBdryPenaltyStrength;
        
    otherwise
        error();
    end
end