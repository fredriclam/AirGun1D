classdef DiscrAirgunShuttleMulti < DiscrAirgun
    properties
        shuttle0            % Initial shuttle state [pos; vel]
        plug0               % Initial plug control volume state [mass; en]
        machAreaFunction    % Precomputed M(A/A*) function interpolant
    end
    
    methods
        function obj = DiscrAirgunShuttleMulti(nx,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth, airgunCrossSectionalArea, REVERT_MODEL, ...
                airgunFiringChamberProfile, ...
                airgunOperatingChamberProfile, bubbleInitialVolume, ...
                shuttleBdryPenaltyStrength)
            if nargin == 8 && ~REVERT_MODEL
                error('Incorrect # of arguments for unreverted model.')
            end
            
            DEBUG = false;
            
            % Call parent constructor
            obj = obj@DiscrAirgun(nx,order,airgunPressure,...
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
            obj.shuttle0 = [1e-6; % [m]
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
            function [dq, dBubble, dShuttle, dPlug, p_RTarget, u_RTarget, monitor] = ...
                    RHS(q, t, bubble, shuttle, plug, p_RTarget, u_RTarget)
                %% Set flag for bypassing plug modeling
                BYPASS_PLUG =  true;
                
                %% Precompute
                % Compute primitive variables at right of PDE domain
                if REVERT_MODEL
                    q_R = schm.e_R'*q;
                else
%                     q_R = q(end-5:end-3); % Go one node in for weak enforcement effects
%                     if sign(q(end-1)) * sign(q_R(2)) < 0
%                         q_R(2) = 0;
%                     end
                    q_R = schm.e_R'*q;
                end

                % Replace p_R with target pressure
                if REVERT_MODEL
                    p_R = schm.p(q_R);
                    u_R =  q_R(2)/q_R(1);
                else
                    % p_RTarget empty means no pressure enforcement
                    if ~isempty(p_RTarget)
                        p_R = p_RTarget;
                    else
                        p_R = schm.p(q_R);
                        if p_R < 0
                            p_R = schm.p(q(end-8:end-6));
                        end
                    end
                    
                    if ~isempty(u_RTarget)
%                         u_R = u_RTarget; % TODO: USE
                        u_R = q_R(2)/q_R(1);
                    else
                        u_R = q_R(2)/q_R(1);
                    end
                    
                                            
                end
                % Default to empty target pressure, velocity constraint
                p_RTarget = [];
                u_RTarget = [];
                
                rho_R = q_R(1);
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

                    if p_R < 0 || rho_R < 0 || T_R < 0
                        % Re-extract from next point in (hopefully some
                        % interlacing property)
                        q_R = q(end-8:end-6);
                        rho_R = q_R(1);
                        en_R = q_R(3);
                        T_R = p_R / rho_R / Q_;
                        c_R = sqrt(physConst.gamma * Q_ * T_R);
                        M_R = u_R / c_R;
                    end
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
%                     isSonicFlags(1) = true;
                end
                
                % reversion: don't trust weak boundary enforcement
                if REVERT_MODEL
                    if M_R >= 1
                        isSonicFlags(1) = true;
                    end
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
                    if BYPASS_PLUG
                        uUpstream = u_R;
                    else
                        uUpstream = u_R - shuttle(2);
                    end
%                     
                    
                    
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
                    charMBubble = bubble(2) / sqrt(gamma_ * Q_ * TBubble);
                    p0bubble = pBubble / pressureMachFunction(gamma_, ...
                            charMBubble);
                    % Compute the critical pressure at which the flow
                    % doesn't go supersonic from downstream
                    pCrit = p0bubble*pressureMachFunction(gamma_,1);
                    
                    % Compute flow throat pressure
                    delayedPMeas = false;
                    if M_R < 1e-5
                        % No-flow case
                            pThroat = p_R;
                    elseif APortExposed < physConst.crossSectionalArea
%                         try
                            % Take port pressure through isentropic expansion
                            sonicArea = physConst.crossSectionalArea / ...
                                areaMachFunction(gamma_, M_R);
                            % HACK: port smaller than sonic area
                            % i.e., round to sonic area
                            % TODO: need to generate a shock
                            if APortExposed / sonicArea < 1
                                APortExposed = sonicArea;
                            end
                            MPort = obj.machAreaFunction( ...
                                APortExposed / sonicArea);
                            pPort = p_R / pressureMachFunction(gamma_, M_R)...
                                * pressureMachFunction(gamma_, MPort);
                            pThroat = pPort;
%                         catch
%                             warning('zing')
%                             delayedPMeas = true;
%                             pThroat = Inf;
%                         end
                    else
                        % Take upstream pressure (chocked at A_cs)
                        pThroat = p_R;
                    end
                    
                    % Case overlap: can't have low pressure sonic expanding
                    % into higher pressure. In practice, the equations
                    % should be self-correcting
                    if((pThroat < pCrit && isSonicFlags(1)))
                        warning('Sonic-pressure contradiction');
                    end
                    
                    if pThroat > pCrit && APortExposed > physConst.crossSectionalArea % isSonicFlags(1)
                        % Resolve subsonic expansion or shock
%                         try
                            % Sonic expanding into subsonic                 TODO: could it be supersonic?                                         
                            caseNum = 3;
                            [velPort, pPort, TPort, rhoPort, cPort, massFlowPort] ...
                            = resolveSonicChamber(obj, APortExposed, uUpstream/c_R, p_R, T_R);
%                         catch e
%                             % Sonic with contraction
%                             caseNum = 4;
%                             % Pass error as warning and resolve
%                             warning(e.message);t
%                             error("Unresolved port-shock case in top level.");
%                         end
                        isSonicFlags(1) = true;
                    elseif pThroat < pCrit
                        caseNum = 1;
                        if DEBUG
                            assert(uUpstream >= 0);
                        end
                        [vel_a, pPort, TPort, rhoPort, cPort, ...
                            massFlowPort, pUpstream] = ...
                            resolveSubsonic(obj, APortExposed, ...
                            uUpstream/c_R, p_R, T_R, ...
                            charMBubble, pBubble);
                        % Airgun velocity is what it is--BC moved to
                        % pressure; should be discarded
                        vel_a = NaN;
                        % Save the target pressure
                        p_RTarget = pUpstream;
                    else
%                         try
                            caseNum = 2;
                            [velUp, pPort, TPort, rhoPort, cPort, massFlowPort]...
                            = resolveSonicPort(obj, APortExposed, uUpstream/c_R, p_R, T_R);
%                             assert(pPort >= pCrit);                         % HACK
                            % Add constraint velocity of plug flow
                            if BYPASS_PLUG
                                vel_a = velUp;
                            else
                                vel_a = velUp + shuttle(2);
                            end
                            u_RTarget = vel_a;
%                             p_RTarget = pPort;
                            isSonicFlags(2) = true;
%                         catch
%                             % Subsonic hack
%                             caseNum = 1;
% %                             assert(pPort >= pBubble);
% %                             assert(physConst.crossSectionalArea > APortExposed)
% %                             assert(~isSonicFlags(1));
%                             %                         assert(uUpstream >= 0); % TODO: find better workaround
%                             [vel_a, pPort, TPort, rhoPort, cPort, ...
%                                 massFlowPort, pUpstream] = ...
%                                 resolveSubsonic(obj, APortExposed, ...
%                                 uUpstream/c_R, p_R, T_R, ...
%                                 APortExposed, ... % Set the jet diameter equal to the port
%                                 bubble(2)/sqrt(gamma_ * Q_ * TBubble), ...
%                                 ... Set jet mach number equal to bubble exp rate
%                                 pBubble, TBubble);
%                             % Airgun velocity is what it is--BC moved to
%                             % pressure
%                             vel_a = u_R;
%                             isSonicFlags(2) = false;
%                             % End hack
%                         end
                    end 
                    
                    if delayedPMeas
                        pThroat = pPort;
                    end
                    
                    
                    
                    % Core block [V1]
%                     if false %~isSonicFlags(1) && p < pBubble%p_R < pCrit
%                         caseNum = 1;
%                         assert(~isSonicFlags(1));
% %                         assert(uUpstream >= 0);                             % TODO: find better workaround
%                         [vel_a, pPort, TPort, rhoPort, cPort, ...
%                             massFlowPort, pUpstream] = ...
%                             resolveSubsonic(obj, APortExposed, ...
%                             uUpstream/c_R, p_R, T_R, ...
%                             APortExposed, ... % Set the jet diameter equal to the port
%                             bubble(2)/sqrt(gamma_ * Q_ * TBubble), ... 
%                             ... Set jet mach number equal to bubble exp rate
%                             pBubble, TBubble);
%                         
%                         % Airgun velocity is what it is--BC moved to
%                         % pressure
%                         vel_a = u_R;
%                     elseif isSonicFlags(1)
%                         caseNum = 2;
%                         [velPort, pPort, TPort, rhoPort, cPort, massFlowPort] ...
%                         = resolveSonicChamber(obj, APortExposed, uUpstream/c_R, p_R, T_R);
%                     elseif false
%                         caseNum = 3;
%                         error('Shock code doesn''t exist yet');
%                     else
%                         try
%                             caseNum = 4;
%                             [velUp, pPort, TPort, rhoPort, cPort, massFlowPort]...
%                             = resolveSonicPort(obj, APortExposed, uUpstream/c_R, p_R, T_R);
%                             assert(pPort >= pCrit);                         % HACK
%                             % Add constraint velocity of plug flow
%                             vel_a = velUp + shuttle(2);
%                             isSonicFlags(2) = true;
%                         catch
%                             % Subsonic hack
%                             caseNum = 1;
% %                             assert(pPort >= pBubble);
% %                             assert(physConst.crossSectionalArea > APortExposed)
% %                             assert(~isSonicFlags(1));
%                             %                         assert(uUpstream >= 0); % TODO: find better workaround
%                             [vel_a, pPort, TPort, rhoPort, cPort, ...
%                                 massFlowPort, pUpstream] = ...
%                                 resolveSubsonic(obj, APortExposed, ...
%                                 uUpstream/c_R, p_R, T_R, ...
%                                 APortExposed, ... % Set the jet diameter equal to the port
%                                 bubble(2)/sqrt(gamma_ * Q_ * TBubble), ...
%                                 ... Set jet mach number equal to bubble exp rate
%                                 pBubble, TBubble);
%                             % Airgun velocity is what it is--BC moved to
%                             % pressure
%                             vel_a = u_R;
%                             isSonicFlags(2) = false;
%                             % End hack
%                         end
%                     end 
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

                %% State regularization
                if M_R < 0 && ~REVERT_MODEL
                    M_R = 0;
                    u_R = 0;
                    warning('Negative mach number')
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
                    % HACK: override pPlug with p_R
                    pPlug = p_R;
                    dShuttle = shuttleEvolve(shuttle, p_R, physConst);
                else
                    dShuttle = 0*shuttle;
                end
                
                % Check for overshoot at weakly enforced sonic BC
%                 assert(~REVERT_MODEL || u_R < 500)

                %% State monitoring
                if ~REVERT_MODEL
                    assert(all(isreal([M_R, u_R, c_R, cPort, ...
                        pThroat, massFlowPort, rhoPort, pPort])));
                    assert(all([M_R, u_R, c_R, cPort, ...
                        pThroat, massFlowPort, rhoPort, pPort, TPort] >= 0));
                    
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
                    monitor(7,1) = pThroat/pCrit;
                    monitor(8,1) = (physConst.APortTotal * ...
                        (posShuttle / physConst.operatingChamberLength))...
                        /physConst.crossSectionalArea;
                    monitor(9,1) = pPlug;
                    monitor(10,1) = p_R;
                    monitor(11,1) = pPort;
                    monitor(12,1) = pCrit;
                    monitor(13,1) = t;
                    monitor(14,1) = dShuttle(2) * physConst.shuttleAssemblyMass * shuttle(2);
                else
                    % Do nothing
                end

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
                            if schm.flowStateR(q) ~= scheme.Euler1d.SUPERSONIC_OUTFLOW
                                dq = dq + closure_r_out_sub_vel(q, vel_a); % TODO: CHECK IF IT'S REAL
                            end
%                             disp(vel_a)
                        else % Subsonic at port: pressure matching
%                             warning('Unchecked subsonic port BC.')
                            dq = dq + closure_r_out_sub(q, pUpstream);
                        end
                    else
%                         warning('Sonic-in-chamber case not treated.')
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
                MUpstream, pUpstream, TUpstream, ...
                MDownstream, pDownstream)
            %% Pre
%             assert(MUpstream >= 0);
            % Workaround: capping the mach number due to overloss to the
            % plug flow
            MUpstream = max([0, MUpstream]);

            %% Take positive part of the exposed area
            APortExposed = max([0, APortExposed]);

            %% Compute upstream pressure from downstream information [V1]
            gamma_ = obj.physConst.gamma;
            Q_ = obj.physConst.Q;
            
            % Compute the stagnation pressure from output
            pTotal = pDownstream / ...
                pressureMachFunction(gamma_, MDownstream);
            % Compute upstream pressure
            pUpstream = pTotal * pressureMachFunction(gamma_, MUpstream);
            % Compute stagnation energy from upstream
            TTotal = TUpstream / temperatureMachFunction(gamma_, MUpstream);             % Q: Which one to use?
            
            % Compute sonic area from upstream
            ASonic = obj.physConst.crossSectionalArea / ...
                areaMachFunction(gamma_, MUpstream);
            % Use relative area for M Port
            MPort = obj.machAreaFunction(APortExposed/ASonic);
            % Compute the rest of the state at port
            TPort = TTotal * temperatureMachFunction(gamma_, MPort);
            cPort = sqrt(gamma_ * Q_ * TPort);
            pPort = pTotal * pressureMachFunction(gamma_, MUpstream);
            rhoPort = pPort / (Q_ * TPort);
            massFlowPort = MPort * cPort * rhoPort * APortExposed;
            
            % Compute chamber exit velocity (just unwinding the Mach #)
            vel_a = MUpstream * sqrt(gamma_ * Q_ * TUpstream);
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
            
            % Numerical tolerancing
            if obj.physConst.crossSectionalArea / A_sonic > 1 ...
                && obj.physConst.crossSectionalArea / A_sonic < 1.2
                A_sonic = obj.physConst.crossSectionalArea;
            end
            
            % Compute chamber outlet velocity as dependent
            if APortExposed <= 0
                vel_a = 0;
            else
                % Compute upstream mach number (set M_a as
                % boundary condition)
                areaRatio = obj.physConst.crossSectionalArea / A_sonic;
                % Fix the area ratio
                areaRatio = max(areaRatio, 1);
                M_a = obj.machAreaFunction(areaRatio);
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
            if M_L > 1
%                 warning('Shuttle velocity acting up!')
                M_L = min(1, M_L);
            end
            assert(APortExposed >= 0 && ...
                   M_L <= 1 && ...
                   p_L > 0 && ...
                   T_L > 0);
            %% Define aliases and constants
            gamma_ = obj.physConst.gamma;
            Q_ = obj.physConst.Q;
            %% Compute port conditions using forward relations
            ASonic = obj.physConst.crossSectionalArea / ...
                areaMachFunction(gamma_, M_L);
            portAreaRatio = APortExposed/ASonic;
            % Check for shock formation due to limited port size
            if portAreaRatio < 1
                portAreaRatio = 1;
%                 warning('Not really shock formation. Redirect.');
            end
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

