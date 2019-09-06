classdef DiscrAirgunShuttle < DiscrAirgun
    properties
        shuttle0           % Initial shuttle state [pos; vel]
        opChamberPressure0 % Initial pressure in op chamber
        port0              % Initial port state [rho; T]
    end
    
    methods
        function obj = DiscrAirgunShuttle(m,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth)
            % Meta flags
            REVERT_MODEL = false; % Skips port area, reverting to LW model
            DEBUG = true; % Debug flag
            
            % Call parent constructor
            obj = obj@DiscrAirgun(m,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth);
            % Alias important objects
            physConst = obj.physConst;
            schm = obj.schm;
            % Replace description
            obj.description = ...
                'Airgun augmented with port-region CV and shuttle';
            
            % Define initial shuttle state: position; velocity
            obj.shuttle0 = [1e-3; 0]; % NOTE: position non-zero or else singular--for CV treatment TODO: check sensitivity
            % Set initial port region state: [p; rho; T]
            % from the initial airgun state
            rho0 = obj.q0(end-2); % Density (rho)
            rv0 = obj.q0(end-1);  % Rho * v
            e0 = obj.q0(end);   % Volumetric stagnation energy e
            T0 = (e0-0.5*rv0^2/rho0)/physConst.c_v/ rho0;
            %             p0 = physConst.p0a; % Consistent!
            obj.port0 = [rho0; 0; T0];
%             obj.port0 = [0; 0; 0]; % Reassigning port to be rho dot, 0
            
            % Set boundary conditions: L: wall; R: pressure when subsonic
            closure_l = schm.boundary_condition('l', 'wall');
            closure_r_out_sub = schm.boundary_condition('r', 'outflow');
            closure_r_out_sub_vel = schm.boundary_condition('r', 'outflow_vel');
            closure_r_out_sub_rhovel = schm.boundary_condition('r', 'outflow_rhovel');
%             % Pressure extraction operator
%             L = @(~,u,~)(physConst.gamma -1)*[0 -1/2*u 1];
            
            %% Redefine RHS to include evolution of shuttle and port-region
            function [dq, dBubble, dShuttle,miscStates] = ...
                    RHS(q,t,bubble,shuttle,miscStates)
                %% Unpacking data
                % Compute primitive variables at right of PDE domain
                q_R = schm.e_R'*q;
                p_R = schm.p(q_R);
                rho_R = q_R(1);
                u_R =  q_R(2)/q_R(1);
                T_R = p_R / rho_R / physConst.Q;
                c_R = sqrt(physConst.gamma * physConst.Q * T_R);
                M_R = u_R / c_R;                % Mach number
                
                % Extract shuttle variables
                posShuttle = shuttle(1);
                velShuttle = shuttle(2);
                % Compute primitive bubble variables
                pBubble = bubblePressure(bubble, physConst);
                TBubble = bubble(4) / physConst.c_v / bubble(3);
                rhoBubble = pBubble / physConst.Q / TBubble;
                
                % Port state
                TPort = T_R; % With temperature always taken from upstream

                % Initialize local flags for sonic/subsonic this timestep
                isSonicFlags = [false, false]; % Internal and port
                
                % Geometry
                % Approximate the total port length as the full travel of
                % the shuttle: the % of the travel is thus the % of the
                % full port area that is exposed
                APortExposed = physConst.APortTotal * ...
                    (posShuttle / physConst.operating_chamber_length);

                %% Catch special flow states
                flowState = schm.flowStateR(q);
                if flowState == scheme.Euler1d.SUBSONIC_INFLOW
                    warning(['Inflow @ t = ' num2str(t) ...
                        '; u|x=0 = ' num2str(u_R)]);
                elseif flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW
                    isSonicFlags(1) = true;
                end
                
                %% Determine port flow state
                % Compute subsonic mass flow rate
                pPort = pBubble; % Match downstream
                rhoPort = pPort / physConst.Q / TPort;
                massFlowPort = physConst.cross_sectional_area * ...
                    (- rhoPort * velShuttle ...
                    + rho_R * u_R);
                % IDEA: replace u_R with near-boundary average and
                % smoothing
                uPortGuess = massFlowPort / rhoPort / APortExposed;
                c_port = sqrt(physConst.gamma * physConst.Q * ...
                    TPort);
                
%                 % Messwithit
% %                 rhoPort = pBubble / physConst.Q / TBubble;
%                 uPortGuess = massFlowPort / rhoPort / APortExposed;
%                 c_port = sqrt(physConst.gamma * physConst.Q * ...
%                     TBubble);
                % Enforce port sonic -> subsonic transition once only
                if t < 0 % 1e-6 % "subsonic" only
                    isSonicFlags(2) = false;
                else    
                    if miscStates(1) == 0 % Accelerating flow
                        if uPortGuess >= c_port
                            % True sonic regime
                            miscStates(1) = 1;
                        end
                        isSonicFlags(2) = false;
                    elseif miscStates(1) == 1 % Stabilization-needed state
    %                     % Can transition to subsonic only
                        if uPortGuess < 0.99*c_port && uPortGuess > 0.0*c_port
                            % True sonic regime
                            miscStates(1) = 2;  %miscStates(1) = 1;
                        else
                            isSonicFlags(2) = true;
                        end
%                         isSonicFlags(2) = true;
                    end

                    if miscStates(1) == 2 % Subsonic flow
                        if uPortGuess >= 1.01*c_port
                            miscStates(1) = 1;
                            isSonicFlags(2) = true;
                        else
                            isSonicFlags(2) = false;
                        end
                    end
                end
                
                % Sonic port treatment: compute choked mass flow rate
                if isSonicFlags(2)
                    % pPort =  TODO: treat pressure properly; it's missing
                    % Recompute port thermodynamics with upstream
                    % properties:
                    
                    pPort = p_R;
                    
                    rhoPort = pPort / physConst.Q / TPort;
                    massFlowPort = c_port * rhoPort * APortExposed;
                end                
                
                %% Build step report (to console)
                % Just building a string
                str = sprintf("TIME %.3e", t);
                if isSonicFlags(1)
                    str = str + " | INT    SONIC";
                else
                    str = str + " | INT SUBSONIC";
                end
                if isSonicFlags(2)
                    str = str + " | PORT    SONIC";
                else
                    str = str + " | PORT SUBSONIC";
                end
                % Data
                str = str + " | SHUTPOS: " + sprintf("%.4e", shuttle(1));
                str = str + " | SHUTVEL: " + sprintf("%.4e", shuttle(2));
                str = str + " | UPORT*: " + sprintf("%.4e", uPortGuess);
                str = str + " | MdotPORT: " + sprintf("%.4e", massFlowPort);
                fprintf(str + "\n");
                
                %% Assert step
                % @post
                % massFlowPort, isSonicFlags[] are valid
                assert(all(shuttle) >= 0);
%                 assert(massFlowPort >= -1e-4)
                
                % DELETE:
%                 if ~REVERT_MODEL
%                     % Check boundary condition
%                     % Cases: sonic correction to boundary condition
%                     % Assemble port region's state evolution vector
%                 else
%                     %% Short-circuit port region dynamics: set == airgun
%                     dPort = [0; 0];
%                     TPort = T_R;
%                     pPort = p_R;
%                     rhoPort = rho_R;
%                     velocityPort = u_R;
%                     ePort = rhoPort * physConst.c_v * TPort;
%                     APortExposed = physConst.APortTotal;
%                     pBubble = bubblePressure(bubble, physConst);
%                 end
                
                %% Shuttle (and operating chamber) dynamics
                % Compute shuttle state evolution // initial ~80g accel
                dShuttle = shuttleEvolve(shuttle, p_R, physConst);
                
                %% Airgun PDE evolution
                % Compute airgun state evolution with left BC
                dq = schm.D(q) + closure_l(q);
                if ~isSonicFlags(1)
                    if isSonicFlags(2) % Sonic at port: momentum BC
                        rhovel_a = massFlowPort / ...
                            physConst.cross_sectional_area + ...
                            rhoPort * velShuttle;
                        dq = dq + closure_r_out_sub_rhovel(q, rhovel_a);
                    else % Subsonic at port: pressure matching
                        dq = dq + closure_r_out_sub(q, pPort);
                    end
                end
                
                %% Bubble evolution
                % Compute port specific stagnation energy
                velocityPort = massFlowPort / rhoPort / APortExposed;
                
                ePort = rhoPort * physConst.c_v * TPort;
                dBubble = bubbleRHS(bubble, ...
                    rhoPort, ...
                    velocityPort, ...
                    ePort, ...
                    pPort, ...
                    physConst.A, ...APortExposed, ...
                    physConst);
                
                plot(q(1:3:end)); hold on;
                drawnow; pause(0.1);
            end
            obj.RHS = @RHS;
        end
    end
end