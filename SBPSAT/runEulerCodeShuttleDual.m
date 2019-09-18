function [sol, q1, bubble1, shuttle1, ...
    q2, bubble2, shuttle2] = ...
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
    shuttle0 = dRevert.shuttle0;
    % Grab ODE's RHS (different)
    RHSRevert = dRevert.RHS;
    RHSShuttle = dNew.RHS;
    
    % Manual hack -- unused channel
    miscStates = [0];
    
    lastPlot = 0;
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
        
        % Update reverted model
        [dq, dBubble, dShuttle, miscStates] = RHSRevert(q1,t,bubble1,shuttle1,miscStates);
        dyRevert = [dq; dBubble; dShuttle];

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
        
        % Update new model
        [dq, dBubble, dShuttle, miscStates] = RHSShuttle(q2,t,bubble2,shuttle2,miscStates);
        dyShuttle = [dq; dBubble; dShuttle];
        
        % Concatenate dy
        dy = [dyRevert; dyShuttle];

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
            lastPlot = t;
        end
    end

    y0 = [q0; bubble0; shuttle0; q0; bubble0; shuttle0];
    
    tspan = [0; 0.080]; % Simulation tmin to tmax (specify here)
%     tspan = [0; 2]; % Simulation tmin to tmax (specify here)
%     options = odeset('RelTol',1e-6);
    options = odeset('RelTol',1e-6);
    
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
    ind2 = ind1+length(q0)-1;
    q2 = sol.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(bubble0)-1;
    bubble2 = sol.y(ind1:ind2, :);

    ind1 = ind2+1;
    ind2 = ind1+length(shuttle0)-1;
    shuttle2 = sol.y(ind1:ind2, :);
%     shuttle = [];
%     port = [];
end
