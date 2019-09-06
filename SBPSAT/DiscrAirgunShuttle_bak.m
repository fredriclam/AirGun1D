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
            REVERT_MODEL = false; % Skips port area, reverting to Watson 2019
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
            obj.shuttle0 = [0; 0];
            % Set initial port region state: [p; rho; T]
            % from the initial airgun state
            rho0 = obj.q0(end-2); % Density (rho)
            rv0 = obj.q0(end-1);  % Rho * v
            e0 = obj.q0(end);   % Volumetric stagnation energy e
            T0 = (e0-0.5*rv0^2/rho0)/physConst.c_v/ rho0;
            %             p0 = physConst.p0a; % Consistent!
            obj.port0 = [rho0; T0];
            obj.port0 = [0; 0]; % Reassigning port to be rho dot, 0
            
            % Set boundary conditions: L: wall; R: pressure when subsonic
            closure_l = schm.boundary_condition('l', 'wall');
            closure_r_out_sub = schm.boundary_condition('r', 'outflow');
            % Pressure extraction operator
            L = @(~,u,~)(physConst.gamma -1)*[0 -1/2*u 1];
            
            %% Redefine RHS to include evolution of shuttle and port-region
            function [dq, dBubble, dShuttle, dPort] = ...
                    RHS(q,t,bubble,shuttle,port)
                %% Shuttle (and operating chamber) dynamics
                % Get pressure at x = 0 (right)
                q_R = schm.e_R'*q;
                % Compute primitive variables at right of PDE domain
                p_R = schm.p(q_R);
                rho_R = q_R(1);
                u_R =  q_R(2)/q_R(1);
                T_R = p_R / rho_R / physConst.Q;
                % Compute shuttle state evolution // initial ~80g accel
                dShuttle = shuttleEvolve(shuttle, p_R, physConst);
                assert(all(shuttle) >= 0);
                
                if ~REVERT_MODEL
                    %% Port region conservation
                    % Flags for coupling regions
                    useVelocityBdry = false;
                    % Unpack and alias states
                    rhoDotPort = port(1);
                    rhoLastPort = port(2);
                    posShuttle = shuttle(1);
                    velShuttle = shuttle(2);
                    pBubble = bubblePressure(bubble, physConst);
                    TBubble = bubble(4) / physConst.c_v / bubble(3);
                    rhoBubble = pBubble / physConst.Q / TBubble;
                    rhoPort = rhoBubble;
                    % Approximating the total port length as the full travel of
                    % the shuttle: the % of the travel is thus the % of the
                    % full port area that is exposed
                    APortExposed = physConst.APortTotal * ...
                        (posShuttle / physConst.operating_chamber_length);
                    
                    % Compute mass flow rate to bubble (intermediate result)
                    massFlowRateToBubble = (rho_R * u_R * ...
                        physConst.cross_sectional_area + ...
                        - (rho * velShuttle + rhoDotPort * velPos) * ...
                        physConst.cross_sectional_area);
                    if APortExposed <= 0
                        massFlowRateToBubble = 0;
                        velocityPort = 0;
                    else
                        velocityPort = massFlowRateToBubble / ...
                            (rhoPort * APortExposed);
                    end
                    % Compute sound speed at port
                    cPort = sqrt(physConst.gamma * physConst.Q * TPort);
                    
                    % Cases: sonic or not sonic at port
                    if APortExposed >= 0 && velocityPort < cPort % :subsonic
                        % Pressure matching with the bubble
                        pPort = pBubble;
                        % Compute new port region density from p, T
                        rhoPort = pPort / physConst.Q / TPort;
                        % Unused density change: assume instant response
                        if DEBUG
                            
                        end
                    elseif APortExposed < 0
                        % Clamp port area to be non-negative
                        APortExposed = max(APortExposed, 0);
                        warning('Port is very closed')
                    else % Fix outflow to sonic if it wants to be supersonic
                        % Choke flow velocity
                        velocityPort = cPort;
                        % Recompute mass flow rate (choked)
                        massFlowRateToBubble = velocityPort * ...
                            (rhoPort * APortExposed);
                        % Compute new port region density from mass balance
                        % %                     if velShuttle > 0
                        % Full mass balance:
                        mdot = ...
                            rho_R * u_R * physConst.cross_sectional_area...
                            - massFlowRateToBubble;
                        assert(mdot >= 0);
                        %                         drho = (mdot / physConst.cross_sectional_area ...
                        %                             - rhoPort * velShuttle) / posShuttle;
                        
                        % Velocity-match PDE region
                        useVelocityBdry = true;
                        uBdry = mdot / rho_R / ...
                            physConst.cross_sectional_area;
                        rhoPort = rho_R;
                        pPort = p_R;
                        
                        % Wrong (linearized about xi = 0):
                        %                         rhoPort = 1/ velShuttle * ...
                        %                             (rho_R * u_R - ...
                        %                             massFlowRateBubble/...
                        %                                 physConst.cross_sectional_area);
                        % %                     else
                        % %                         drho = 0;
                        % %                     end
                        % Compute (discontinuous pressure) from rho, T
                        %                     pPort = rhoPort * physConst.Q * TPort;
                        % CHECK: does this ^ work?
                    end
                    % Compute temperature evolution
                    if posShuttle ~= 0
                        dT = -physConst.gamma * TPort * ...
                            velShuttle / posShuttle + ...
                            physConst.gamma / rhoPort / posShuttle * ...
                            (rho_R * u_R * T_R - ...
                            massFlowRateToBubble / ...
                            physConst.cross_sectional_area / ...
                            TBubble);
                    else
                        dT = 0;
                    end
                    
                    % Assemble port region's state evolution vector
                    TPort = TBubble;
                    pPort = pBubble;
                    drho = 0;
                    dT = 0;
                    dPort = [drho; dT];
                    assert(all(~isnan(dPort)));
                else
                    %% Short-circuit port region dynamics: set == airgun
                    dPort = [0; 0];
                    TPort = T_R;
                    pPort = p_R;
                    rhoPort = rho_R;
                    velocityPort = u_R;
                    ePort = rhoPort * physConst.c_v * TPort;
                    APortExposed = physConst.APortTotal;
                    pBubble = bubblePressure(bubble, physConst);
                end
                
                
                if DEBUG
                    disp([num2str(pPort) '    ' num2str(bubble(1))]);
                end
                
                %% Airgun
                % Get flow state at right end of PDE domain
                flowState = schm.flowStateR(q);
                % If the flow is inward
                if flowState == scheme.Euler1d.SUBSONIC_INFLOW
                    warning('Subsonic inflow condition.');
                end
                % Compute airgun state evolution with left BC
                dq = schm.D(q) + closure_l(q);
                
                %% Apply conditions on airgun and calculate q_hat
                if flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW
                    % No bc required
                else % flowState == scheme.Euler1d.SUBSONIC_OUTFLOW
                    % Pressure matching: airgun and port region
                    %                 dq = dq + closure_r_out_sub(q,pPort); % Rather than pBubble, use pPort
                    dq = dq + closure_r_out_sub(q,pBubble); % Rather than pBubble, use pPort
                end
                
                %% Bubble evolution
                % Compute port specific stagnation energy
                ePort = rhoPort * physConst.c_v * TPort;
                dBubble = bubbleRHS(bubble, ...
                    rhoPort, ...
                    velocityPort, ...
                    ePort, ...
                    pPort, ...
                    physConst.A, ...APortExposed, ...
                    physConst);
            end
            obj.RHS = @RHS;
        end
    end
end