classdef DiscrAirgunShuttle < DiscrAirgun
    properties
        shuttle0           % Initial shuttle state [pos; vel]
        opChamberPressure0 % Initial pressure in op chamber
        port0              % Initial port state [rho; T]
%         logger             % Event message string
    end
    
    methods
        function obj = DiscrAirgunShuttle(m,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth)
            % Call parent constructor
            obj = obj@DiscrAirgun(m,order,airgunPressure,...
                airgunLength,airgunPortArea,...
                airgunDepth);
            
            %% Replacing parent constructor
            default_arg('m',100)
            default_arg('order',4)

            % Extract physConst object from parent construction
            physConst = obj.physConst;
%             [physConst, t0, icAirgun, icBubble] = configAirgun('Bolt1500LL',airgunPressure,airgunLength,airgunPortArea,airgunDepth);
%             A  = physConst.A;
% Compute number of cells based on input m (cells per m)
            xlim = {-physConst.L, 0};
            m = ceil((xlim{2}-xlim{1})*m+1);
            schm = scheme.Euler1d(m,xlim,order,[],[],true);
            
%             fprintf('Starting at time t0 = %f\n',t0)

            
%   
%             closure_l = schm.boundary_condition('l', 'wall');
%             closure_r_out_sub = schm.boundary_condition('r', 'outflow');  %% should set pressure
%             L = @(~,u,~)(physConst.gamma -1)*[0 -1/2*u 1];
%             
%             
%             H = schm.H;
% 
%             obj.t0 = t0;
% 
%             rho0 = icAirgun.rho0(schm.u);
%             rv0 = icAirgun.rv0(schm.u);
%             e0 = icAirgun.e0(schm.u);
% 
%             obj.q0 = reshape([rho0'; rv0'; e0'],3*m,1);
% 
% 
%             obj.bubble0 = [
%                 icBubble.R;
%                 icBubble.Rdot;
%                 icBubble.m;
%                 icBubble.E;
%             ];
% 
%             obj.order = order;
%             obj.RHS = @RHS;
%             obj.H = H;
%             obj.schm = schm;
%             obj.h = schm.h;
%             obj.physConst = physConst;
            
            
            %% 
            
            
            
            % Replace description
            obj.description = ...
                'Airgun augmented with port-region CV and shuttle';
%             % Initialize logger (character array for events)
%             obj.logger = '';
            
            
            % Define initial shuttle state: position; velocity
            obj.shuttle0 = [0; 0];
            % Set initial port region state: [p; rho; T]
            % from the initial airgun state
            rho0 = obj.q0(end-2); % Density rho
            rv0 = obj.q0(end-1);  % Rho * v
            e0 = obj.q0(end);   % Volumetric stagnation energy e
            T0 = (e0-0.5*rv0^2/rho0)/physConst.c_v/ rho0;
            % p0 = physConst.p0a; % Consistent!
            obj.port0 = [rho0; T0];
            
            schm = scheme.Euler1d(m,xlim,order,[],[],true);
            closure_l = schm.boundary_condition('l', 'wall');
            closure_r_out_sub = schm.boundary_condition('r', 'outflow');  %% should set pressure
            L = @(~,u,~)(physConst.gamma -1)*[0 -1/2*u 1];
            
            %% Redefine RHS to include evolution of shuttle and port-region
            function [dq, dBubble, dShuttle, dPort] = ...
                    RHS(q,t,bubble,shuttle,port)

                %% Shuttle (op chamber)
                % Get pressure at x = 0 (right) 
                q_R = schm.e_R'*q;
                % Compute primitive variables at right of PDE domain
                p_R = schm.p(q_R);
                rho_R = q_R(1);
                u_R =  q_R(2)/q_R(1);
                T_R = p_R / rho_R / physConst.Q;
                % Compute shuttle state evolution // initial ~80g accel
                dShuttle = shuttleEvolve(shuttle, p_R, physConst);
                
                %% Port region
                % Unpack and alias states
                rhoPort = port(1);
                TPort = port(2);
                posShuttle = shuttle(1);
                velShuttle = shuttle(2);
                pBubble = bubblePressure(bubble, physConst);
                massBubble = bubble(3);
                energyBubble = bubble(4);
                TBubble = bubble(4) / physConst.c_v / bubble(3);
                % Approximating the port length as the full travel of the
                % shuttle: the % of the travel is thus the % of the full
                % port area
                APortExposed = physConst.APortTotal * ...
                    (posShuttle / physConst.operating_chamber_length);
                
                % Compute mass flow rate to bubble (intermediate result)
                massFlowRateBubble = rhoPort * (u_R - velShuttle) * ...
                    physConst.cross_sectional_area;
                % Compute resulting port velocity
                velocityPort = massFlowRateBubble / ...
                    (rhoPort * APortExposed);
                % Compute sound speed at port
                cPort = sqrt(physConst.gamma * physConst.Q * TPort);
                
                % Cases: sonic or not sonic at port
                if APortExposed >= 0 && velocityPort < cPort % :subsonic
                    % Pressure matching
                    pPort = pBubble;
                    % Compute new port region density from p, T
                    rhoPort = pPort / physConst.Q / TPort;
                else % Fix outflow to sonic if it wants to be supersonic
                    % Clamp port area to be non-negative
                    APortExposed = max(APortExposed, 0);
                    % Choke flow velocity
                    velocityPort = cPort;
                    % Recompute mass flow rate (choked)
                    massFlowRateBubble = velocityPort * ...
                        (rhoPort * APortExposed);
                    % Compute new port region density from mass balance
                    if velShuttle > 0
                        rhoPort = 1/ velShuttle * ...
                            (rho_R * u_R - ...
                            massFlowRateBubble/...
                                physConst.cross_sectional_area);
                    else
                    end
                    % Compute (discontinuous pressure) from rho, T
                    pPort = rhoPort * physConst.Q * TPort;
                    % TODO: does this ^ work?
                end
                % Unused density change: assume instant response
                drho = 0;
                % Compute temperature evolution
                if posShuttle ~= 0
                    dT = -physConst.gamma * TPort * ...
                        velShuttle / posShuttle + ...
                        physConst.gamma / rhoPort / posShuttle * ...
                        (rho_R * u_R * T_R - ...
                        massFlowRateBubble / ...
                        physConst.cross_sectional_area / ...
                        TBubble);
                else
                    dT = 0;
                end
                % Assemble port region's state evolution vector
                dPort = [drho; dT];
                
                %% Airgun
                % Get flow state at right end of PDE domain
                flowState = schm.flowStateR(q);
                % If the flow is inward
                if flowState == scheme.Euler1d.SUBSONIC_INFLOW
%                     dq = q.*0;
%                     dBubble = bubbleRHS(bubble, 0, 0, 0, 0, 0, physConst);
                    warning('Subsonic inflow condition.');
%                     return
                end
                % Compute airgun state evolution with left BC
                dq = schm.D(q) + closure_l(q);
                
                %% Apply conditions on airgun and calculate q_hat
%                 p_b = bubblePressure(bubble, physConst);
                if flowState == scheme.Euler1d.SUBSONIC_OUTFLOW
                    % Pressure matching: airgun and port region
                    dq = dq + closure_r_out_sub(q,pPort); % Rather than pBubble, use pPort
                    
%                     % Set
%                     pIn = [3];
%                     pOut = [1 2];
%                     permutation = [pIn pOut];
%                     invPermutation(permutation) = 1:3
%                     
%                     % q_r = schm.e_R'*q;
%                     % w = inv(schm.T(q_r))*q_r;
%                     
%                     % Ltemp = L(q_r(1), q_r(2)/q_r(1), q_r(3));
%                     % g = p_b;
%                     % T = schm.T(q_r);
%                     
%                     % wHat = zeros(3,1);
%                     % wHat(pIn) = inv(Ltemp*T(:,pIn))*(g - Ltemp*T(:,pOut)*w(pOut)); % g := p_b
%                     % wHat(pOut) = w(pOut);
%                     
%                     % qHat = T*wHat;
                elseif flowState == scheme.Euler1d.SUPERSONIC_OUTFLOW
                    % No bc required
%                     qHat = schm.e_R'*q;
                else
                    % Will not run due to earlier if
                    q_r = schm.e_R'*q;
                    printExpr('q_r(1)');
                    printExpr('q_r(2)');
                    printExpr('q_r(3)');
                    
                    % q
%                     c = schm.c(q_r)
%                     v = q_r(2)/q_r(1)
                    error('Undefined behaviour. flowState = %d, t = %f', flowState, t)
                end
                % Extract conserved variables at the right boundary
%                 qHat = schm.e_R'*q;
%                 rho_a = qHat(1);
%                 v_a = qHat(2)/qHat(1);
%                 e_a = qHat(3);
%                 p_a = schm.p(qHat);
                %dBubble = bubbleRHS(bubble, rho_a, v_a, e_a, p_a, A, physConst);
                
%                 % turbulent energy dissipation (uncomment to use)
%                 C = 0; %coefficient of turbulent energy dissipation
%                 gun_area = 12.5;
%                 port_area = 12;
%                 deltaP = C*rho_a*v_a^2*(gun_area/port_area)^2;
%                 % Bubble evolution
%                 dBubble = bubbleRHS(bubble, rho_a, v_a, e_a, p_a - deltaP, A, physConst);

                %% Bubble evolution
                % Compute port specific stagnation energy
                ePort = rhoPort * physConst.c_v * TPort;
                
                dBubble = bubbleRHS(bubble, ...
                    rhoPort, ...
                    velocityPort, ...
                    ePort, ...
                    pPort, ...
                    APortExposed, ...
                    physConst);
            end
            obj.RHS = @RHS;
        end
        
%         % Append line to logger
%         function obj = log(msg)
%             obj.logger = [obj.logger '\n' msg];
%         end
    end
end