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
            obj.shuttle0 = [0.001; 0]; % NOTE: position non-zero or else singular TODO: check sensitivity
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
%             % Pressure extraction operator
%             L = @(~,u,~)(physConst.gamma -1)*[0 -1/2*u 1];
            
            %% Redefine RHS to include evolution of shuttle and port-region
            function [dq, dBubble, dShuttle, dPort] = ...
                    RHS(q,t,bubble,shuttle,port)
                %% Unpacking data
                % Compute primitive variables at right of PDE domain
                q_R = schm.e_R'*q;
                p_R = schm.p(q_R);
                rho_R = q_R(1);
                u_R =  q_R(2)/q_R(1);
                T_R = p_R / rho_R / physConst.Q;
                c_R = sqrt(physConst.gamma * physConst.Q * T_R);
                % Extract shuttle variables
                posShuttle = shuttle(1);
                velShuttle = shuttle(2);
                % Compute primitive bubble variables
                pBubble = bubblePressure(bubble, physConst);
                TBubble = bubble(4) / physConst.c_v / bubble(3);
                rhoBubble = pBubble / physConst.Q / TBubble;
                
                % Port state
                TPort = T_R; % Temperature always taken from upstream
                % REVISION 0.0
                
%                 % Extract port variables % Old2019-08-09
%                 rhoPort = port(1);
%                 rhoDotPort = port(2);
%                 TPort = port(3);

                % Initialize local flags for sonic/subsonic this timestep
                isSonicFlags = [false, false]; % Internal and port

                %% Catch special flow states
                flowState = schm.flowStateR(q);
                if flowState == scheme.Euler1d.SUBSONIC_INFLOW
                    warning(['Inflow @ t = ' num2str(t) ...
                        '; u|x=0 = ' num2str(u_R)]);
                elseif flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW
                    isSonicFlags(1) = true;
                end
                
                %% Determine port flow state
                % Treat as subsonic:
                pPort = pBubble; % Try bubble pressure
                rhoPort = pPort / physConst.Q / TPort;
                massFlowPort = A_cs * (rhoPort * velShuttle - ...
                    rho_R * u_R);
                
                % Compute port average velocity and check against sonic
                uPortGuess = massFlowPort / rhoPort / A_port;
                c_port = sqrt(physConst.gamma * physConst.Q * ...
                        TPort);
                if uPortGuess >= c_port
                    disp('Port supersonic detected...')
                    isSonicFlags(2) = true;
                else
                    % Explicitly set sonic for port to false
                    isSonicFlags(2) = false;
                end
                
                % Sonic port treatment: compute choked mass flow rate
                if isSonicFlags(2)
                    % pPort =  TODO: treat pressure properly; it's missing
                    % Recompute port thermodynamics with upstream
                    % properties:
                    pPort = p_R;
                    rhoPort = pPort / physConst.Q / TPort;
                    massFlowPort = c_port * rhoPort * A_port;
                end
                
                % @post
                % massFlowPort, isSonicFlags[] are valid
                
                %% Build step report (to console)
                
                %% Shuttle (and operating chamber) dynamics
                % Compute shuttle state evolution // initial ~80g accel
                dShuttle = shuttleEvolve(shuttle, p_R, physConst);
                % Assert no negative position, velocity
                assert(all(shuttle) >= 0);
                
                %% 
                
                if ~REVERT_MODEL
                    %% Port region conservation
                    % Flags for coupling regions
                    useVelocityBdry = false;
                    
                    

%                     % Compute rho-dot of the bubble (same rho as output)
%                     RBubble = bubble(1);
%                     dRBubble = bubble(2);
%                     dmBubble = dBubble(3);
%                     dVBubble = (4*pi*RBubble^2) * dRBubble;
%                     VBubble = 4/3*pi*RBubble^3;
%                     rhoDotBubble = (dmBubble - rhoBubble * ...
%                         dVBubble) / VBubble;
                    
                    % Approximate the total port length as the full travel of
                    % the shuttle: the % of the travel is thus the % of the
                    % full port area that is exposed
                    APortExposed = physConst.APortTotal * ...
                        (posShuttle / physConst.operating_chamber_length);
                    
                    % Compute mass flow rate to bubble (intermediate result)
                    massFlowRateToBubble = (rho_R * u_R * ...
                        physConst.cross_sectional_area + ...
                        - (rhoPort * velShuttle + rhoDotPort * posShuttle) * ...
                        physConst.cross_sectional_area);
                    if massFlowRateToBubble < 0
                       disp('!') 
                    end
                    velocityPort = massFlowRateToBubble / ...
                            (rhoPort * APortExposed);
                        
                    if APortExposed <= 0
                        massFlowRateToBubble = 0;
                        velocityPort = 0;
                        % Recompute compression effects
                        useVelocityBdry = true;
                        uBdry = 1 / rho_R * (rhoPort * velShuttle + ...
                            rhoDotPort * posShuttle) + ...
                            massFlowRateToBubble / physConst.cross_sectional_area;
                    end

                    % Compute sound speed at port
                    c_port = sqrt(physConst.gamma * physConst.Q * TPort);
                    
                    % Cases: sonic correction to boundary condition
                    if velocityPort < c_port % :subsonic
                        % Pressure matching with the bubble
                        pPort = pBubble;
                    else % Fix outflow to sonic if it wants to be supersonic
                        % Choke flow velocity
                        velocityPort = c_port;
                        % Recompute mass flow rate (choked)
                        massFlowRateToBubble = velocityPort * ...
                            (rhoPort * APortExposed);
                        useVelocityBdry = true;
                        uBdry = 1 / rho_R ...
                            / physConst.cross_sectional_area * ...
                            (massFlowRateToBubble + ...
                            (rhoPort * velShuttle + ...
                            rhoDotPort * posShuttle) * ...
                            physConst.cross_sectional_area);
                        % Sonic in interior
                        if uBdry > c_R
                            warning('doubly sonic...logically improbable')
                            uBdry = c_R;
                        end
%                         % Full mass balance:
%                         mdot = ...
%                             rho_R * u_R * physConst.cross_sectional_area...
%                             - massFlowRateToBubble;
%                         assert(mdot >= 0);
%                         %                         drho = (mdot / physConst.cross_sectional_area ...
%                         %                             - rhoPort * velShuttle) / posShuttle;
%                         
%                         % Velocity-match PDE region
%                         useVelocityBdry = true;
%                         uBdry = mdot / rho_R / ...
%                             physConst.cross_sectional_area;
%                         rhoPort = rho_R;
%                         pPort = p_R;
%                         
%                         % Wrong (linearized about xi = 0):
%                         %                         rhoPort = 1/ velShuttle * ...
%                         %                             (rho_R * u_R - ...
%                         %                             massFlowRateBubble/...
%                         %                                 physConst.cross_sectional_area);
%                         % %                     else
%                         % %                         drho = 0;
%                         % %                     end
%                         % Compute (discontinuous pressure) from rho, T
%                         %                     pPort = rhoPort * physConst.Q * TPort;
%                         % CHECK: does this ^ work?
                    end
                    % Compute temperature evolution
%                     if posShuttle ~= 0
%                         dT = -physConst.gamma * TPort * ...
%                             velShuttle / posShuttle + ...
%                             physConst.gamma / rhoPort / posShuttle * ...
%                             (rho_R * u_R * T_R - ...
%                             massFlowRateToBubble / ...
%                             physConst.cross_sectional_area / ...
%                             TBubble);
%                     else
%                         dT = 0;
%                     end
%                     assert(abs(dT) < 1e9);
% %                     
                    % Assemble port region's state evolution vector
                    dT = 0;
                    dPort = [0; 0; dT];
                    TPort = T_R;
                    pPort = pBubble;
%                     dPort = [drho; dT];
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
                
                
                % Compute airgun state evolution with left BC
                dq = schm.D(q) + closure_l(q);
                
                %% Apply conditions on airgun and calculate q_hat
                if flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW
                    % No bc required
                else % flowState == scheme.Euler1d.SUBSONIC_OUTFLOW
                    % Pressure matching: airgun and port region
                    %                 dq = dq + closure_r_out_sub(q,pPort); % Rather than pBubble, use pPort
                    if useVelocityBdry
                        dq = dq + closure_r_out_sub_vel(q,uBdry);
                    else
                        % Classic pressure bdry: pPort == pBubble
                        dq = dq + closure_r_out_sub(q,pPort); % Rather than pBubble, use pPort
                    end
                end
                
%                 if DEBUG
%                     disp([num2str(pPort) '    ' num2str(bubble(1))]);
%                 end
                
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