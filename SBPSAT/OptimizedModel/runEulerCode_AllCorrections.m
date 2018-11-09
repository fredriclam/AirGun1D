function sol = runEulerCode_AllCorrections(nx,airgunPressure,airgunLength,airgunPortArea,airgunDepth, dissp_boolean, beta, Mfac, Vinit)
%sol = runEulerCode(nx)
    d = DiscrAirgun_AllCorrections(nx,3,airgunPressure,airgunLength,airgunPortArea,airgunDepth, dissp_boolean, beta, Mfac, Vinit);
    
    q0 = d.q0;
    t0 = d.t0;
    bubble0 = d.bubble0;
    RHS = d.RHS;
    
    N = length(bubble0);
    
    function dy = odefun(t,y)
        bubble = y(1:N);
        q = y(N+1:end);

        [dq, dBubble] = RHS(q,t,bubble);
        dy = [dBubble; dq];
    end

    y0 = [bubble0; q0];
    tspan = [0; 2];
    options = odeset('RelTol',1e-6);
    
    sol = ode45(@odefun, tspan, y0,options);

    %k = 0.000003/2
    %sol.y = LeightonsRK4(@odefun,tspan(2),k,y0);
    %sol.x = 0:k:tspan(2);
        
    %sol.q = q;
end