function [physConst, t0, icAirgun, icBubble] = configAirgun(str,airgunPressure,airgunLength,airgunPortArea,airgunDepth)

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

%         cross_sectional_area = 24; % [in^2] cross sectional area is calculated by assuming 1m long (39.3701) for 600in3 airgun.
% %         cross_sectional_area = airgunPortArea;
%         airgunLengthImperial = physConst.L/0.0254; % length in inches
%         V_imperial = airgunLengthImperial * cross_sectional_area;
%         V = V_imperial * 1.63871e-5; % air gun volume [m^3]

        % Rewriting airgun geometry section
        % Port area input is currently the cross-sectional area rather than
        % the max area of the circumferential port
        crossSectionalAreaImperial = airgunPortArea; % Was manually 12.5; % [in^2] (case in Watson long paper) 
        airgunLengthImperial = airgunLength/0.0254; % [in] driven by airgunLength
        V_imperial = crossSectionalAreaImperial * airgunLengthImperial; % [cui]
        V = V_imperial * 1.63871e-5; % air gun volume [m^3]
        
        A_imperial = airgunPortArea; % air gun port area [in^2]
        A = A_imperial * 6.4516e-4; % air gun port area [m^2]
        physConst.A = A;
        
        %% Water constants
        physConst.rho_inf = 1e3; % density [kg/m^3]
        physConst.pa = 1e5; % atmospheric pressure [Pa]
        physConst.Tinf = 288; % temperature assumed constant throughout the system [K]
        depth = airgunDepth; % depth [m]
        g = 9.8; % gravitational acceleration [m/s^2]
        physConst.p_inf = physConst.pa + physConst.rho_inf*g*depth; % ambient pressure at depth [Pa]
        
        %% (Air calorically perfect gas) constants
        physConst.c_v = 718; % heat capacity of air at constant volume [J/kgK]
        physConst.c_inf = 1482; % speed of sound in water [m/s]
        physConst.Q = 287.06; % specific gas constant for dry air [J/kgK]
        physConst.R_G = physConst.Q; % DUPLICATE
        physConst.gamma = 1.4; % ratio of heat capacities for dry air
        
        %% Fixed open-time model
        % Time when air gun stops firing [s]
        physConst.AirgunCutoffTime = 0.04; % Should be 0.010, but specified as 0.04 in airgun1d/HEAD
        
        % Air gun
        p = physConst.p0a;
        T = physConst.Tinf;
        Q = physConst.Q;
        c_v = physConst.c_v;
        
        rho = p/(Q*T);
        e = c_v*rho*T;
        
        %% Airgun initial conditions
        icAirgun.rho0 = @(x)0*x + rho;
        icAirgun.rv0  = @(x)0*x;
        icAirgun.e0   = @(x)0*x + e;
        icAirgun.p0   = @(x)0*x + p;
        
        %% Bubble initial conditions and gas air parameters
        p = physConst.p_inf;
        T = physConst.Tinf;
        Q = physConst.Q;
        c_v = physConst.c_v;
        
        % Make the widely-used V_bubble = V_airgun assumption
        %V_airgun_const = 250 * 1.63871e-5;
%         V_airgun_const = 600 * 1.63871e-5;

        % HACK: Overwrites current V with fixed prescribed volume as in
        % airgun1d/HEAD
        V = 600 * 1.63871e-5; 
        VBubbleInitial = V; % Equal to airgun volume [m^3]
        
        icBubble.R = (3/(4*pi) * VBubbleInitial)^(1/3);
        icBubble.Rdot = 0;
        icBubble.m = p*VBubbleInitial / (Q*T); 
        icBubble.E = c_v * icBubble.m * T;
        
        % Additions by Fred
        % Convert cross_sectional_area from [in^2] to [m^2]
        crossSectionalArea = crossSectionalAreaImperial*0.00064516;
        physConst.cross_sectional_area = crossSectionalArea;
        physConst.operating_chamber_length = 3 *2.54/100; % [m]
        physConst.p_R0 = p0a; % [Pa] // left pressure = right pressure
        physConst.shuttleAssemblyMass = 63 * .454; % 63/.454; % [kg]
        physConst.APortTotal = 0.0464; % [m^2] // Credits Ludvig
        
        % Ludvig paramters for airgun area:
%         A_L = 0.0586; % [m^2]
%         A_R = 0.0506; % [m^2]
        % Artificial shuttle left and right areas
        physConst.shuttle_area_left = crossSectionalArea; % [m^2]
        % Reasonable guess to shuttle right area: ratio using Ludvig's
        % code, although the raw areas are different
        % Note: in practice the area on the right is higher, but when the
        % valve is clicked an additional pressure channel makes the
        % effective area on the left higher
        physConst.shuttle_area_right = 0.0506/0.0586*crossSectionalArea; % [m^2]
        % Override; HACK: using arbitrary number for fun
%         physConst.shuttle_area_right = 0.8*crossSectionalArea; % [m^2]

        % Local numerical parameters
        physConst.shuttleBdryPenaltyStrength = 1e11; % 1e7; % [ ... ]

    otherwise
        error();
    end
end