classdef DiscrAirgunShuttleMulti < DiscrAirgun
    properties
        shuttle0           % Initial shuttle state [pos; vel]
        opChamberPressure0 % Initial pressure in op chamber
    end
    
    methods
        function obj = DiscrAirgunShuttleMulti(m,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth, REVERT_MODEL)
            % Meta flags
% %             REVERT_MODEL = true; % Skips port area, reverting to LW model
            
            % Call parent constructor
            obj = obj@DiscrAirgun(m,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth);
            % Alias commonly used objects
            physConst = obj.physConst;
            schm = obj.schm;
            gamma_ = physConst.gamma;
            % Replace this object's description
            obj.description = ...
                'Airgun augmented with port-region CV and shuttle';
            
            % Define initial shuttle state: position; velocity
            % NOTE: position non-zero or else singular--for CV treatment
            % TODO: check sensitivity
            % TODO: incorporate initial volume encompassed by shuttle but
            % not relevant to the port area calculation
            obj.shuttle0 = [1e-6;
                            0];

            % Set initial port region state: [p; rho; T]
            % from the initial airgun state
            rho0 = obj.q0(end-2); % Density (rho)
            rv0 = obj.q0(end-1);  % Rho * v
            e0 = obj.q0(end);     % Volumetric stagnation energy e
            T0 = (e0-0.5*rv0^2/rho0)/physConst.c_v/ rho0; % Temperature
            p0 = physConst.p0a;   % Initial pressure
            
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
            
            % Isentropic steady-flow thermodynamic ratios:
            % Ratio p/p0 (p at M to stagnation pressure)
            pressure_ratio = @(M) (1 + 0.5*(gamma_-1)* M .^2 ).^ ...
                (-gamma_/(gamma_-1));
            % Ratio T/T0 (T at M to stagnation temperature)
            temperature_ratio = @(M) 1./ (1 + 0.5*(gamma_-1)* M .^2 );
            % A/A* as function of Mach number
            area_ratio = @(M) ...
                ((gamma_+1)/2)^(-(gamma_+1)/2/(gamma_-1)) * ...
                (1 + (gamma_-1)/2 * M.^2 ).^ ...
                ((gamma_+1)/2/(gamma_-1)) ...
                ./ M;
            
            %% Redefine RHS to include evolution of shuttle and port-region
            function [dq, dBubble, dShuttle,miscStates] = ...
                    RHS(q,t,bubble,shuttle,miscStates)
                %% Precompute
                % Compute primitive variables at right of PDE domain
                q_R = schm.e_R'*q;
                p_R = schm.p(q_R);
                rho_R = q_R(1);
                u_R =  q_R(2)/q_R(1);
                en_R = q_R(3);
                T_R = p_R / rho_R / physConst.Q;
                c_R = sqrt(physConst.gamma * physConst.Q * T_R);
                M_R = u_R / c_R;
                % Compute bubble variables
                pBubble = bubblePressure(bubble, physConst);
                TBubble = bubble(4) / physConst.c_v / bubble(3);
                rhoBubble = pBubble / physConst.Q / TBubble;
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
                        (posShuttle / physConst.operating_chamber_length);
                end

                % Initialize local flags for sonic/subsonic this timestep
                isSonicFlags = [false, false]; % Internal and port
                % Capture flow state in PDE domain
                flowState = schm.flowStateR(q);
                if flowState == scheme.Euler1d.SUBSONIC_INFLOW
                    warning(['Inflow @ t = ' num2str(t) ...
                        '; u|x=0 = ' num2str(u_R)]);
                elseif flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW
                    isSonicFlags(1) = true;
                end
                
                % Compute ratio of chamber cross-sectional area to sonic
                % area of the flow for control
                A_on_Asonic_ratio = ...
                    ((gamma_+1)/2)^(-(gamma_+1)/2/(gamma_-1)) * ...
                    (1 + (gamma_-1)/2*M_R^2 )^((gamma_+1)/2/(gamma_-1)) ...
                    ./ M_R;
                % Compute sonic area of flow
                A_sonic = physConst.cross_sectional_area / ...
                    A_on_Asonic_ratio;
                % Compute sonic area ratio at the port
                APortOnASonic = APortExposed/A_sonic;
                
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
                    % Sonic port
                    A_sonic = APortExposed;
                    % Use positive part of APortExposed only
                    if APortExposed <= 0
                        APortExposed = 0;
                        vel_a = 0;
                    else
                        % Compute upstream mach number (set M_a as
                        % boundary condition)
                        M_a = fzero( @(M) ...
                            ((gamma_+1)/2)^(-(gamma_+1)/2/(gamma_-1)) * ...
                            (1 + (gamma_-1)/2 * M^2 )^ ...
                            ((gamma_+1)/2/(gamma_-1)) ./ M - ...
                            physConst.cross_sectional_area/A_sonic, ...
                            [1e-10,1-1e-10]);
                        % TODO: try using just velocity BC as a substitute
                        vel_a = M_a * c_R;
                    end
                    
                    % Proposed model: if bubble pressure > port pressure, the
                    % flow is not consistent; switch to subsonic
                    % TODO: check code for subsonic case
                    if p_R < pBubble
                        error('p_R < pBubble. Should do subsonic flow');
                        
                        % Inverse problem for finding subsonic mach number at port
                        obj_func = @(M_port) APortExposed / ...
                            physConst.cross_sectional_area- ...
                            area_ratio(M_port) / area_ratio(M_R);
                        % Special case: M is zero
                        if ~isfinite(area_ratio(M_R))
                            M_port = 0;
                        else
                            % Solve for subsonic mach root
                            M_port = fzero(obj_func, [1e-10,1 - 1e-10]);
                            assert(M_port >= 0 && M_port <= 1);
                        end
                        % Compute thermodynamic properties using forward
                        % relations
                        pPort = pressure_ratio(M_port) / ...
                            pressure_ratio(M_R) * p_R;
                        % Compute pressure upstream through isentropic
                        % flow relation
                        p_a = pressure_ratio(M_R) / ...
                            pressure_ratio(M_port) * pBubble;
                        % Compute temperature downstream through isentropic
                        % flow relation
                        TPort = temperature_ratio(M_port) / ...
                            temperature_ratio(M_R) * T_R;
                        rhoPort = pPort / physConst.Q / TPort;
                        % Averaged density over the control volume
                        % (substitute)
                        rhoPortAveraged = mean([rho_R, rhoPort]);
                        % Suppose the mass flow rate is determined by inflow to
                        % control volume minus volumetric dilation; use
                        % averaged density between two endpoints of isentropic
                        % flow as approximation to the bulk-averaged density
                        massFlowPort = physConst.cross_sectional_area * ...
                            (- rhoPort * velShuttle ...
                            ...(- rhoPortAveraged * velShuttle ...
                            + rho_R * u_R);                    
                        % Explicit subsonic
                        isSonicFlags(2) = false;
                    else                    
                        % Flag port as sonic
                        isSonicFlags(2) = true;
                    end
                    
                    % Compute choked properties at port
                    % TODO: check significance of these ratios
                    pPort = pressure_ratio(1) / pressure_ratio(M_R) * p_R;
                    TPort = temperature_ratio(1) / temperature_ratio(M_R) * T_R;
                    rhoPort = pPort / physConst.Q / TPort;
                    cPort = sqrt(physConst.gamma * physConst.Q * ...
                        TPort);
                    massFlowPort = rhoPort * cPort * APortExposed;
                    
                    % Sanity check
                    assert(TPort > 0 && pPort > 0 && rhoPort >0 && ...
                        cPort > 0 && isreal(cPort))
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
                str = str + ...
                    " | APORT/A*: " + sprintf("%.4e", APortOnASonic);
%                 fprintf(str + "\n");
                
                %% Shuttle (and operating chamber) dynamics
                % Compute shuttle state evolution
                % Early data: should be initial ~80g accel
                dShuttle = shuttleEvolve(shuttle, p_R, physConst);
                
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
                            warning('Unchecked subsonic port BC.')
                            dq = dq + closure_r_out_sub(q, p_a);
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
            end
            obj.RHS = @RHS;
        end
        
        
        % Compute port dynamics based on mass conservation on control
        % volume (0-dimensional model with average rho or T dependent on
        % upstream or downstream parameters).
        function massConservationBasedPhysics()
            error('Don''t use me yet.')
            if false % Mass conservation choking
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
               
                % Compute rho * velocity in airgun
                rhovel_a = massFlowPort / ...
                    physConst.cross_sectional_area + ...
                    rhoPort * velShuttle;
            end
        end
    end
end

