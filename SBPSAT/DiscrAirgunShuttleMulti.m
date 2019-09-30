classdef DiscrAirgunShuttleMulti < DiscrAirgun
    properties
        shuttle0            % Initial shuttle state [pos; vel]
        plug0               % Initial plug control volume state [mass; en]
        machAreaFunction    % Precomputed M(A/A*) function interpolant
    end
    
    methods
        function obj = DiscrAirgunShuttleMulti(m,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth, airgunCrossSectionalArea, REVERT_MODEL, ...
                airgunFiringChamberProfile, ...
                airgunOperatingChamberProfile, bubbleInitialVolume, ...
                shuttleBdryPenaltyStrength)
            if nargin == 8 && ~REVERT_MODEL
                error('Incorrect # of arguments for unreverted model.')
            end
            
            % Call parent constructor
            obj = obj@DiscrAirgun(m,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth);
            if ~REVERT_MODEL
                % Override the configuration with updated, backward-compatible
                % setting
                [physConst, t0, icAirgun, icBubble] = ...
                    configAirgun('GeneralAirgun', ...
                        airgunPressure, ...
                        airgunLength, ...
                        airgunPortArea, ...
                        airgunDepth, ...
                        airgunCrossSectionalArea, ...
                        airgunFiringChamberProfile, ...
                        airgunOperatingChamberProfile, ...
                        bubbleInitialVolume, ...
                        shuttleBdryPenaltyStrength);
                obj.physConst = physConst;
            end
            
            % Alias commonly used objects
            physConst = obj.physConst;
            schm = obj.schm;
            gamma_ = physConst.gamma;
            Q_ = physConst.Q;
            % Replace this object's description
            obj.description = ...
                'Airgun augmented with port-region CV and shuttle';
            
            % Set initial port region state: [p; rho; T]
            % from the initial airgun state
            rho0 = obj.q0(end-2); % Density (rho)
            rv0 = obj.q0(end-1);  % Rho * v
            e0 = obj.q0(end);     % Volumetric stagnation energy e
            T0 = (e0-0.5*rv0^2/rho0)/physConst.c_v/ rho0; % Temperature
            p0 = physConst.p0a;   % Initial pressure
            
            % Define initial shuttle state: position; velocity
            % NOTE: position non-zero or else singular--for CV treatment
            % TODO: check sensitivity
            % TODO: incorporate initial volume encompassed by shuttle but
            % not relevant to the port area calculation
            obj.shuttle0 = [1e-3; % [m]
                            0];   % [m/s]

            if ~REVERT_MODEL
                obj.plug0 = [rho0*physConst.plugVolumeFn(obj.shuttle0(1)); % [kg]
                    rho0*physConst.plugVolumeFn(obj.shuttle0(1))* ...
                    physConst.c_v*T0];                      % [J]
            else
                obj.plug0 = [0; 0];
            end
                       
                     
            %% Create boundary condition operators
            closure_l = schm.boundary_condition('l', 'wall');
            % Outflow with pressure (for unchoked-everywhere flow)
            closure_r_out_sub = schm.boundary_condition('r', 'outflow');
            % Outflow with velocity (for choked flow)
            closure_r_out_sub_vel = ...
                schm.boundary_condition('r', 'outflow_vel');
            % Outflow with rho*velocity (alternative to the latter; unused)
            closure_r_out_sub_rhovel = ...
                schm.boundary_condition('r', 'outflow_rhovel');
            % Wall on the right
            closure_r_closed = schm.boundary_condition('r', 'wall');
            % Precompute mach area function M(A/A*) and, for subsonic,
            % the mach pressure function M(p/p0)
            if ~REVERT_MODEL
                obj.machAreaFunction = precomputeMachAreaFunction(gamma_);
            end
            
            %% Redefine RHS to include evolution of shuttle and port-region
            function [dq, dBubble, dShuttle, dPlug, monitor] = ...
                    RHS(q, t, bubble, shuttle, plug)
                %% Precompute
                % Compute primitive variables at right of PDE domain
                q_R = schm.e_R'*q;
                p_R = schm.p(q_R);
                rho_R = q_R(1);
                u_R =  q_R(2)/q_R(1);
                en_R = q_R(3);
                T_R = p_R / rho_R / Q_;
                c_R = sqrt(physConst.gamma * Q_ * T_R);
                M_R = u_R / c_R;
                % Compute bubble variables
                pBubble = bubblePressure(bubble, physConst);
                TBubble = bubble(4) / physConst.c_v / bubble(3);
                rhoBubble = pBubble / Q_ / TBubble;
                % Extract shuttle variables
                posShuttle = shuttle(1);
                velShuttle = shuttle(2);
                % Compute geometry
                if REVERT_MODEL
                    % Fixed outlet area
%                     APortExposed = physConst.APortTotal;
                    % Fix outlet area to equal the cross-sectional area
                    APortExposed = physConst.A;
                else
                    % Approximate the total port length as the full travel of
                    % the shuttle: the % of the travel is thus the % of the
                    % full port area that is exposed
                    APortExposed = physConst.APortTotal * ...
                        (posShuttle / physConst.operatingChamberLength);
%                     APortExposed = 0;
                end

                % Initialize local flags for sonic/subsonic this timestep
                isSonicFlags = [false, false]; % Internal and port
                
                % Capture flow state in PDE domain
                flowState = schm.flowStateR(q);
                if flowState == scheme.Euler1d.SUBSONIC_INFLOW
                    if ~REVERT_MODEL || t < physConst.AirgunCutoffTime
                        warning(['Inflow @ t = ' num2str(t) ...
                            '; u|x=0 = ' num2str(u_R)]);
                    end
                elseif flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW
                    isSonicFlags(1) = true;
                end

                if REVERT_MODEL
                    if t <= physConst.AirgunCutoffTime
                        massFlowPort = rho_R * u_R * APortExposed;
                        velocityPort = u_R;
                    else
                        massFlowPort = 0;
                        velocityPort = 0;
                    end
                    rhoPort = rho_R;
                    pPort = p_R;
                    TPort = T_R;
                else
                    % Plug flow:
                    % Split velocity at constant pressure to the plug
                    uUpstream = u_R - shuttle(2);
                    
                    %% Case split
                    % Case 1: subsonic everywhere
                    % Conditions:
                    %   - the internal pressure is insufficient to expand
                    %     to the characteristic outflow pressure
                    % Should be an edge case indicating insufficient
                    % internal pressure.
                    % Effects:
                    %   - bubble pressure determines flow state upstream
                    
                    % Case 2: internally sonic
                    % Conditions:
                    %   - the resulting flow to the port must be possible
                    %     (i.e., must be converging in absence of friction)
                    % Effects:
                    %   - the port flow is determined
                    
                    % Case 3: firing chamber shock formation
                    % Conditions:
                    %   - the flow was internally sonic, but the resulting
                    %     steady flow to the port was not possible
                    % e.g. if the flow is choked inside and the port is
                    % closed too rapidly
                    % Effects:
                    %   - the port becomes sonic
                    %   - shock formation must come through Euler domain BC
                    
                    % Case 4: sonic at port
                    % Conditions:
                    %   - All else fails
                    % Should be the most commonly encountered case.
                    
                    % Compute critical back pressure at characteristic
                    % bubble velocity
                    
%                     pCrit = pressureMachFunction(bubble(2) / ...
%                             sqrt(gamma_ * Q_ * TBubble));
%                     ABubble = 4 * pi * bubble(1).^2; (gamma_ * Q_);
%                     ASonicAs = APortExposed;
                    % Stagnation pressure required for choked flow at neck,
                    % estimated using characteristic bubble velocities
                    %                                                       % TODO: Check using p0, rather than p_L
                    pCrit = pBubble / pressureMachFunction(gamma_, ...
                            bubble(2) / sqrt(gamma_ * Q_ * TBubble));
                    if p_R < pCrit
                        caseNum = 1;
                        assert(~isSonicFlags(1));
%                         assert(uUpstream >= 0);                             % TODO: find better workaround
                        [vel_a, pPort, TPort, rhoPort, cPort, ...
                            massFlowPort, pUpstream] = ...
                            resolveSubsonic(obj, APortExposed, ...
                            uUpstream/c_R, p_R, T_R, ...
                            APortExposed, ... % Set the jet diameter equal to the port
                            bubble(2)/sqrt(gamma_ * Q_ * TBubble), ... 
                            ... Set jet mach number equal to bubble exp rate
                            pBubble, TBubble);
                        
                        % Airgun velocity is what it is--BC moved to
                        % pressure
                        vel_a = u_R;
                    elseif isSonicFlags(1)
                        caseNum = 2;
                        [velPort, pPort, TPort, rhoPort, cPort, massFlowPort] ...
                        = resolveSonicChamber(obj, APortExposed, uUpstream/c_R, p_R, T_R);
                    elseif false
                        caseNum = 3;
                        error('Shock code doesn''t exist yet');
                    else
                        caseNum = 4;
                        [velUp, pPort, TPort, rhoPort, cPort, massFlowPort]...
                        = resolveSonicPort(obj, APortExposed, uUpstream/c_R, p_R, T_R);
                        % Add constraint velocity of plug flow
                        vel_a = velUp + shuttle(2);
                        isSonicFlags(2) = true;
                    end 
                end

                %% Console output
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
%                 str = str + " | UPORT*: " + sprintf("%.4e", uPortGuess);
                str = str + ...
                    " | MdotPORT: " + sprintf("%.4e", massFlowPort);
%                 str = str + ...
%                     " | APORT/A*: " + sprintf("%.4e", APortOnASonic);
%                 fprintf(str + "\n");
                
                %% State monitoring
                if ~REVERT_MODEL
                    monitor(1,1) = M_R;
                    if APortExposed > 0
                        uPort = massFlowPort/(rhoPort * APortExposed);
                    else
                        uPort = 0;
                    end
                    monitor(2,1) = uPort / sqrt(gamma_ * Q_ * TPort);
                    monitor(3,1) = caseNum;
                    monitor(4,1) = u_R;
                    monitor(5,1) = uPort;
                    monitor(6,1) = shuttle(2);
                else
                    % Do nothing
                end

                %% Shuttle (and operating chamber) dynamics
                % Compute shuttle state evolution
                % Early data: should be initial ~80g accel
                if ~REVERT_MODEL
                    volPlug = physConst.plugVolumeFn(shuttle(1));
                    assert(volPlug > 0);
                    assert(plug(1) > 0);
                    assert(plug(2) > 0);
                    rhoPlug = plug(1) /  volPlug;
                    TPlug = plug(2) / plug(1) / physConst.c_v;
                    pPlug = rhoPlug * Q_ * TPlug;                           % TODO: output and complete
                    
%                     dShuttle = shuttleEvolve(shuttle, p_R, physConst);
                    dShuttle = shuttleEvolve(shuttle, pPlug, physConst);
                else
                    dShuttle = 0*shuttle;
                end
                
                % Check for overshoot at weakly enforced sonic BC
%                 assert(~REVERT_MODEL || u_R < 500)
                
                %% Airgun PDE evolution
                % Compute airgun state evolution with left BC
                dq = schm.D(q) + closure_l(q);
                if REVERT_MODEL
                    if velocityPort == 0 && t > 0.0
                        dq = dq + closure_r_closed(q);
                        assert(~REVERT_MODEL || t >= 0.010)
                    elseif ~isSonicFlags(1) % Subsonic (in chamber)
                        dq = dq + closure_r_out_sub(q, pBubble);
                    end
                else
                    if ~isSonicFlags(1) % Subsonic in chamber
                        if vel_a == 0
                            dq = dq + closure_r_closed(q);
                        elseif isSonicFlags(2) % Sonic at port: momentum BC
                            dq = dq + closure_r_out_sub_vel(q, vel_a);
                        else % Subsonic at port: pressure matching
%                             warning('Unchecked subsonic port BC.')
                            dq = dq + closure_r_out_sub(q, pUpstream);
                        end
                    else
                        warning('Sonic-in-chamber case not treated.')
                    end
                end
                
                %% Bubble evolution
                % Compute port velocity
                if APortExposed > 0
                    velocityPort = massFlowPort / rhoPort / APortExposed;
                else
                    velocityPort = 0;
                end
                % Compute port specific stagnation energy
                if REVERT_MODEL
                    ePort = en_R;
                else
                    % Compute total energy per volume
                    ePort = rhoPort * physConst.c_v * TPort + ...
                        0.5 * rhoPort * velocityPort^2;
                end
                % Compute bubble differential
                dBubble = bubbleRHS(bubble, ...
                    rhoPort, ...
                    velocityPort, ...
                    ePort, ...
                    pPort, ...
                    APortExposed, ...physConst.A, ...APortExposed, ...
                    physConst);
                
                %% Plug evolution
                if REVERT_MODEL
                    dPlug = [0;0];
                else
                    dPlug = [nan; nan];
                    dPlug(1) = rho_R * u_R * physConst.crossSectionalArea;
                    dPlug(2) = dPlug(1) * physConst.c_v * T_R + ...
                        0.5 * dPlug(1) * u_R.^2;
                end
            end
            obj.RHS = @RHS;
        end
        
        
        
        
        % Subsonic port, subsonic firing chamber
        % Bubble parameters matter; information can propagate upstream.
        % Note that at the bubble interface (assumed spherical) the
        % velocity is equal to the bubble expansion velocity \dot{R}.
        % 
        % We assume the expansion is isentropic, and furthermore that
        % shocks are absent (although the subsonic resolution may be
        % invoked right after the float becomes unchoked, and in reality
        % there could still be reflecting shocks present).
        %
        % Computes the velocity boundary condition on the PDE domain
        % Input:
        %   obj
        %   APortExposed   Area of port at current opening state
        %   M_R            Upstream mach number
        %   p_R            Upstream pressure
        %   T_R            Upstream temperature
        %   ADownstream
        %   MDownstream
        %   pDownstream
        %   TDownstream
        % Output:
        %   vel_a          Upstream velocity
        %   pPort          Port pressure
        %   TPort          Port temperature
        %   rhoPort        Port density
        %   cPort          Port sound speed
        %   massFlowPort   Port mass flow rate
        %   pUpstream      Upstream pressure
        % Pre:
        %   APortExposed positive or negative (positive part taken).
        %   M_R > 0
        function [vel_a, pPort, TPort, rhoPort, cPort, massFlowPort, pUpstream] ...
                 = resolveSubsonic(obj, APortExposed, ...
                M_R, p_R, T_R, ...
                ADownstream, MDownstream, pDownstream, TDownstream)
            %% Pre
%             assert(M_R >= 0);
            % Workaround:
            M_R = max([0, M_R]);
            %% Take positive part of the exposed area
            APortExposed = max([0, APortExposed]);
            %% Compute upstream pressure from downstream information
            gamma_ = obj.physConst.gamma;
            Q_ = obj.physConst.Q;
            
            ARatio = areaMachFunction(gamma_, MDownstream);
            ASonic = ADownstream / ARatio;
%             assert(APortExposed / ASonic <= 1);
            
            pTotal = pDownstream / ...
                pressureMachFunction(gamma_, MDownstream);
            TTotal = TDownstream / ...                                      % TODO: compare to upstream T
                temperatureMachFunction(gamma_, MDownstream);
            MUpstream = obj.machAreaFunction(...
                obj.physConst.crossSectionalArea / ASonic);
            pUpstream = pTotal * pressureMachFunction(gamma_, MUpstream);
            
            % Enforced density
            rho_R = pUpstream / (Q_ * T_R);
            u_R = M_R * sqrt(gamma_ * Q_ * T_R);
                                                                            % TODO: Check p_R - pUpstream violation
            massFlowPort = rho_R * u_R * obj.physConst.crossSectionalArea;  % Compare value to port mass flow rate
            
            %% Compute port properties
            MPort = obj.machAreaFunction(...
                APortExposed / ASonic);
            pPort = pTotal * pressureMachFunction(gamma_, MPort);
            TPort = TTotal * temperatureMachFunction(gamma_, MPort);
            rhoPort = pPort / (Q_ * TPort);
            cPort = sqrt(gamma_ * Q_ * TPort);
            
            % Compute chamber exit velocity (just unwinding the Mach #)
            vel_a = M_R * sqrt(gamma_ * Q_ * T_R);
        end        
        
        % Sonic port, subsonic firing chamber
        % Computes the velocity boundary condition on the PDE domain
        % Input:
        %   obj
        %   APortExposed   Area of port at current opening state
        %   M_R            Upstream mach number
        %   p_R            Upstream pressure
        %   T_R            Upstream temperature
        % Output:
        %   vel_a          Upstream velocity
        %   pPort          Port pressure
        %   TPort          Port temperature
        %   rhoPort        Port density
        %   cPort          Port sound speed
        %   massFlowPort   Port mass flow rate
        % Pre:
        %   APortExposed positive or negative (positive part taken).
        function [vel_a, pPort, TPort, rhoPort, cPort, massFlowPort]...
                = resolveSonicPort(obj, APortExposed, M_R, p_R, T_R)
            %% Compute velocity boundary condition
            gamma_ = obj.physConst.gamma;
            Q_ = obj.physConst.Q;
            % Set sonic area to be the positive part of the
            % current port area
            A_sonic = max([0,APortExposed]);            
            % Compute chamber outlet velocity as dependent
            if APortExposed <= 0
                vel_a = 0;
            else
                % Compute upstream mach number (set M_a as
                % boundary condition)
                M_a = obj.machAreaFunction(...
                    obj.physConst.crossSectionalArea / A_sonic);
                vel_a = M_a * sqrt(gamma_ * Q_ * T_R);
            end            
            %% Compute choked properties at port
                                                                            % TODO: check significance of these ratios
            pPort = p_R * pressureMachFunction(gamma_, 1) ...
                / pressureMachFunction(gamma_, M_R);
            TPort = T_R * temperatureMachFunction(gamma_, 1)...
                / temperatureMachFunction(gamma_, M_R);
            rhoPort = pPort / Q_ / TPort;
            cPort = sqrt(gamma_ * Q_ * ...
                TPort);
            massFlowPort = rhoPort * cPort * APortExposed;
            
            %% Sanity check
            assert(TPort > 0 && ...
                   pPort > 0 && ...
                   rhoPort > 0 && ...
                   cPort > 0 && ...
                   isreal(cPort))
        end
        
        
        % Subsonic port, sonic firing chamber (subsonic inflow)
        % Computes the velocity boundary condition on the PDE domain
        % Input:
        %   obj
        %   APortExposed   Area of port at current opening state
        %   M_L            Upstream mach number
        %   p_L            Upstream pressure
        %   T_L            Upstream temperature
        % Output:
        %   velPort        Downstream velocity
        %   pPort          Port pressure
        %   TPort          Port temperature
        %   rhoPort        Port density
        %   cPort          Port sound speed
        %   massFlowPort   Port mass flow rate
        % Pre:
        %   APortExposed >= 0
        %   M_L <= 1
        %   p_L > 0
        %   T_L > 0
        function [velPort, pPort, TPort, rhoPort, cPort, massFlowPort] ...
                = resolveSonicChamber(obj, APortExposed, M_L, p_L, T_L)
            %% Precheck
            assert(APortExposed >= 0 && ...
                   M_L <= 1 && ...
                   p_L > 0 && ...
                   T_L > 0);
            %% Compute port conditions using forward relations
            ASonic = obj.physConst.crossSectionalArea / ...
                areaMachFunction(M_L);
            portAreaRatio = APortExposed/ASonic;
            % Check for shock formation due to limited port size
            if portAreaRatio < 1
                error('Shock formation. Redirect.');
            end
            gamma_ = obj.physConst.gamma;
            Q_ = obj.physConst.Q;
            MPort = obj.machAreaFunction(portAreaRatio);
            pPort = p_L * pressureMachFunction(gamma_, MPort) / ...
                pressureMachFunction(gamma_, M_L);
            TPort = T_L * temperatureMachFunction(gamma_, MPort) / ...
                temperatureMachFunction(gamma_, M_L);
            rhoPort = pPort / Q_ / TPort;
            
            cPort = sqrt(gamma_ * Q_ * TPort);
            velPort = cPort * MPort;            
            massFlowPort = rhoPort * velPort * APortExposed;
        end
    end    
end

