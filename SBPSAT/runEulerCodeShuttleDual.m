function [sol, q1, bubble1, shuttle1, plug1, ...
    q2, bubble2, shuttle2, plug2, monitorStates] = ...
    runEulerCodeShuttleDual(nx,airgunPressure,...
    airgunLength,airgunCrossSecArea,airgunPortArea,airgunDepth, ...
    airgunFiringChamberProfile, ...
    airgunOperatingChamberProfile, bubbleInitialVolume, ...
    shuttleBdryPenaltyStrength)

    % Initialize airgun models (reverted and new model)
    orderSBP = 3;
    dRevert = DiscrAirgunShuttleMulti(nx,orderSBP,airgunPressure,airgunLength, ...
        airgunPortArea,airgunDepth,0,true);
    dNew = DiscrAirgunShuttleMulti(nx,orderSBP,airgunPressure,airgunLength, ...
        airgunCrossSecArea,airgunDepth,airgunPortArea,false, ...
        airgunFiringChamberProfile, ...
        airgunOperatingChamberProfile, bubbleInitialVolume, ...
        shuttleBdryPenaltyStrength);
    
    % Grab initial states for each subsystem (same for both systems)
    q0 = dRevert.q0;
    bubble0 = dRevert.bubble0;
    shuttle0 = dNew.shuttle0;
    plug0 = dNew.plug0;
    % Grab ODE's RHS (different)
    RHSRevert = dRevert.RHS;
    RHSShuttle = dNew.RHS;
    
    lastPlot = 0;
    monitorStates = [];
    p_RTarget = [];
    u_RTarget = [];
    function dy = odefun(t,y)
        % Manually partition y into each subsystem
        ind1 = 1;
        ind2 = length(q0);
        q1 = y(ind1:ind2);
        
        ind1 = ind2+1;
        ind2 = ind1+length(bubble0)-1;
        bubble1 = y(ind1:ind2);
        
        ind1 = ind2+1;
        ind2 = ind1+length(shuttle0)-1;
        shuttle1 = y(ind1:ind2);
        
        ind1 = ind2+1;
        ind2 = ind1+length(plug0)-1;
        plug1 = y(ind1:ind2);
        
        % Update reverted model
        [dq, dBubble, dShuttle, dPlug, ~] = RHSRevert(q1,t,bubble1,shuttle1,plug1);
        dyRevert = [dq; dBubble; dShuttle; dPlug];

        % Extract concatenated data
        ind1 = ind2+1;
        ind2 = ind1+length(q0)-1;
        q2 = y(ind1:ind2);
        
        ind1 = ind2+1;
        ind2 = ind1+length(bubble0)-1;
        bubble2 = y(ind1:ind2);
        
        ind1 = ind2+1;
        ind2 = ind1+length(shuttle0)-1;
        shuttle2 = y(ind1:ind2);
        
        ind1 = ind2+1;
        ind2 = ind1+length(plug0)-1;
        plug2 = y(ind1:ind2);
        
        % Update new model
        [dq, dBubble, dShuttle, dPlug, p_RTarget, u_RTarget, monitor] = ...
            RHSShuttle(q2,t,bubble2,shuttle2,plug2,p_RTarget,u_RTarget);
        dyShuttle = [dq; dBubble; dShuttle; dPlug];
        
        % Concatenate dy
        dy = [dyRevert; dyShuttle];
        
        % DEBUG: state monitoring
        monitorStates = [monitorStates, monitor];

        % DEBUG: Console output
%         disp(bubble(1));
%         disp(['t = ' num2str(t)]);
% plot(q(2:3:length(q)) ./ q(1:3:length(q)))

        %% Piggyback plotting every millisecond
        if t - lastPlot >= 1e-3
            figure(99); clf;
            subplot(3,1,1);
            v1 = q1(2:3:end) ./ q1(1:3:end);
            plot(dRevert.schm.x(2:3:end), v1);
            hold on
            v2 = q2(2:3:end) ./ q2(1:3:end);
            plot(dNew.schm.x(2:3:end), v2);
            
            ylabel('u [m/s]')
            title(num2str(t))
            
            subplot(3,1,2);
            T1 = (q1(3:3:end) - 0.5 * q1(2:3:end).^2 ./ q1(1:3:end)) ./ q1(1:3:end) ./dRevert.physConst.c_v;
            plot(dRevert.schm.x(2:3:end), T1);
            hold on
            T2 = (q2(3:3:end) - 0.5 * q2(2:3:end).^2 ./ q2(1:3:end)) ./ q2(1:3:end) ./dNew.physConst.c_v;
            plot(dNew.schm.x(2:3:end), T2);
            ylabel('T [K]')
            
            subplot(3,1,3);
            p1 = (dRevert.physConst.gamma - 1) * (q1(3:3:end) - 0.5 * q1(2:3:end).^2 ./ q1(1:3:end));
            plot(dRevert.schm.x(2:3:end), p1);
            ylabel('p [Pa]')
%             equilibriumPressure = ...
%                 dRevert.physConst.shuttle_area_right / ...
%                 dRevert.physConst.cross_sectional_area * ...
%                 dRevert.physConst.p_R0;
            hold on
%             plot(dRevert.schm.x(2:3:end), equilibriumPressure*ones(size(dRevert.schm.x(2:3:end))));
            
            p2 = (dNew.physConst.gamma - 1) * (q2(3:3:end) - 0.5 * q2(2:3:end).^2 ./ q2(1:3:end));
            plot(dNew.schm.x(2:3:end), p2);
            
            drawnow;
            
            figure(98);
            m_rear = shuttle2(3);
            E_rear = shuttle2(4);
            m_front = shuttle2(5);
            E_front = shuttle2(6);
            
            % compute pressure
            % Compute mass flow rates based on choked or OFF # TODO: implement
% sub-choked flow rates
V_rear = dNew.physConst.shuttle_area_right * shuttle2(1) + 1e-10; % A quick approximation
% Constrained density
rho_rear = m_rear / V_rear;
% Constrained temperature
T_rear = E_rear / (m_rear * dNew.physConst.c_v);
% Pressure
p_rear = rho_rear * dNew.physConst.Q * T_rear;

V_front_max = dNew.physConst.shuttle_area_right * dNew.physConst.operatingChamberLength;

% Compute mass flow rates based on choked or OFF # TODO: implement
% sub-choked flow rates
V_front = V_front_max - V_rear; % A quick approximation
% Constrained density
rho_front = m_front / V_front;
% Constrained temperature
T_front = E_front / (m_front * dNew.physConst.c_v);
% Pressure
p_front = rho_front * dNew.physConst.Q * T_front;
            
            subplot(2,2,1)
            bar([m_rear, m_front])
            title 'm'
            subplot(2,2,2)
            bar([E_rear, E_front])
            title 'E'
            subplot(2,2,3)
            bar([p_rear, p_front])
            title 'p'
            subplot(2,2,4)
            bar([T_rear, T_front]) 
            title 'T'
            drawnow;
            
            lastPlot = t;
        end
    end

    % HACK: Override bubble volume
    bubble0revert = bubble0;
    y0 = [q0; bubble0; shuttle0; plug0; q0; bubble0revert; shuttle0; plug0];
    
    tspan = [0; 0.010]; % Simulation tmin to tmax (used for set-test)
    tspan = [0; 0.030]; % Simulation tmin to tmax (used for set-test)
%     tspan = [0; 0.100]; % Simulation tmin to tmax (used for set-test)
%     tspan = [0; 0.400]; % Simulation tmin to tmax (used for set-test)

%     tspan = [0; 0.250]; % Simulation tmin to tmax (specify here)
%     tspan = [0; 2]; % Simulation tmin to tmax (specify here)
%     options = odeset('RelTol',1e-6);
    options = odeset('RelTol',1e-3);%, 'MaxStep', 1e-3);
    
    sol = ode45(@odefun, tspan, y0,options);

% %     %k = 0.000003/2
% %     %sol.y = LeightonsRK4(@odefun,tspan(2),k,y0);
% %     %sol.x = 0:k:tspan(2);
% %         
% %     %sol.q = q;

    ind1 = 1;
    ind2 = length(q0);
    q1 = sol.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(bubble0)-1;
    bubble1 = sol.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(shuttle0)-1;
    shuttle1 = sol.y(ind1:ind2, :);
    
    ind1 = ind2+1;
    ind2 = ind1+length(plug0)-1;
    plug1 = sol.y(ind1:ind2, :);
    
    ind1 = ind2+1;
    ind2 = ind1+length(q0)-1;
    q2 = sol.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(bubble0)-1;
    bubble2 = sol.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(shuttle0)-1;
    shuttle2 = sol.y(ind1:ind2, :);
    
    ind1 = ind2+1;
    ind2 = ind1+length(plug0)-1;
    plug2 = sol.y(ind1:ind2, :);
    
    
    
%     assert(ind2 == length(y0));
%     shuttle = [];
%     port = [];
end
