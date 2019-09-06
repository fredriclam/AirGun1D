% Derived from makeFig02_SimulationCWT

clear all; clc;
nx = 200;  % Number of grid points per 1 m of air gun length
% Useless: tmin = 0;  % Simulation time min
% Useless: tmax = 10; % Simulation time max

r = 7.5;         % Distance from source to receiver [m]
c_inf = 1482;   % Speed of sound in water [m/s]
rho_inf = 1000; % Density of water [kg/m^3]

% in2_m2 = 0.00064516; % conversion from in^2 to m^2
% psi_pa = 6894.76; % conversion from psi to pa
% 
% gamma = 1.4; % ratio of heat capacities
% Q = 287.06; % specific gas constant for dry air [J/kgK]
% T_inf = 288; % temperature assumed constant throughout the system [K]

aP = 2000; % air gun pressure for slope/rise time plot
aL = 0.6; % air gun length [m]
aA = 16; % air gun port area [in^2] % cross-sectional area = port area
aD = 7.5; % air gun depth [m]

% run solve
[sol, q, bubble, shuttle, ...
    q2, bubble2, shuttle2] = runEulerCodeShuttle(nx, aP, aL, aA, aD);

t = sol.x; % time

%% Post processing

[pPres, R] = computePressure(bubble, sol, t, rho_inf, c_inf, r, aD);
[pPres2, R2] = computePressure(bubble2, sol, t, rho_inf, c_inf, r, aD);


%% %%% plotting %%%    

%% Chamber properti

%% Bubble radius
figure(1); clf;
tAxis = 1000*sol.x;
plot(tAxis, R);
xlabel('Time [ms]');
ylabel('Bubble radius [m]');
hold on
plot(tAxis, R2);

%% Chamber density
figure(3); clf;
rho = q(1:3:end,:);
xAxis = linspace(0,1.2,size(rho, 1))';
tAxis = 1000*sol.x;
[xx, tt] = meshgrid(xAxis, tAxis);
contourf(xx, tt, rho', 'LineStyle', 'none');
xlabel('Position [m]')
ylabel('Time [ms]')
caxis([0, 160]);
colorbar;

windowPos = get(gcf,'position');
set(gcf,'position',[windowPos(1), windowPos(2), ...
    560, 280]);
colormap bone;

title('rho [kg / m^3]')

%% Chamber particle velocity
figure(4); clf;
u = q(2:3:end,:) ./ q(1:3:end,:);
xAxis = linspace(0,1.2,size(rho, 1))';
tAxis = 1000*sol.x;
[xx, tt] = meshgrid(xAxis, tAxis);
contourf(xx, tt, rho', 'LineStyle', 'none');
xlabel('Position [m]')
ylabel('Time [ms]')
caxis([0, 350]);
colorbar;
windowPos = get(gcf,'position');
set(gcf,'position',[windowPos(1), windowPos(2), ...
    560, 280]);
colormap bone;

title('u [m/s]')

%% Boundary velocity record
% Compute first time the shuttle closes
closingIndex = find(shuttle(1,:) <= 0, 1);
figure(5); clf;
subplot(3,1,1);
u_R = u(end,:);
plot(tAxis, u_R);
title('Boundary velocity [m]')
hold on
yL = ylim;
plot(tAxis(closingIndex)*[1,1], [-100 200], 'k')
ylim(yL)

subplot(3,1,2);
plot(tAxis, shuttle(1,:));
title('Virtual shuttle pos [m]')

hold on
yL = ylim;
plot(tAxis(closingIndex)*[1,1], [-100e-3 200e-3], 'k')
ylim(yL)

subplot(3,1,3);
plot(tAxis, shuttle(2,:));
title('Virtual shuttle velocity [m]')

hold on
yL = ylim;
plot(tAxis(closingIndex)*[1,1], [-100 100], 'k')
ylim(yL)

%%
figHand1 = figure(2); clf;
figPos = get(figHand1,'Position');
set(figHand1,'Position',[figPos(1) figPos(2) 600 900]);

%% bubble volume
V = 4/3*pi*R.^3;
subplot(4,1,1);
plot(t*1000,V);
% xmax = 500;
% xlim([0 xmax]);
% hold on;
xlabel('Time, t (ms)');
ylabel('Volume (m^3)');
% ylim([0 1.5])

%% acoustic pressure
subplot(4,1,2);
h = plot((tInterp-r/c_inf)*1000, pPres*1e-5*r);
% xlim([0 xmax]);
% set(h.Parent,'YTick',[-4 0 4]);
ylabel('\Delta p (bar m)');
xlabel('Time, t-r/c_\infty (ms)');
% ylim([-4 6])

%% CWT %%
if false
    Fs = 1/dt;
    [wt,f] = cwt(pPres,Fs);
    [m,n] = size(wt);
    subplot(4,1,[3 4]);
    h = pcolor(repmat(tInterp-r/c_inf,m,1)*1000, repmat(f,1,n), abs(wt));
    shading interp
    colormap parula
    ylim([5 220]);
    currentXLim = xlim;
    xlim([0 currentXLim(2)]);
    xlabel('Time, t-r/c_\infty (ms)');
    ylabel('Frequency (Hz)');
end

% handaxes1 = axes('Position',[0.47 0.26 0.4 0.2]);
% set(handaxes1,'YColor',[0 0 0]);
% h = plot((tInterp-r/c_inf)*1000, pPres*1e-5*r);
% xlim([0 2])
% ylim([0 6])
% set(h.Parent,'YColor',[1 1 1]);
% set(h.Parent,'XColor',[1 1 1]);
% xlabel('Time (ms)');
% ylabel('\Delta p (bar m)');
% hold on;
% [xx,idx] = max(pPres*1e-5*r);
% % g = hline(xx);
% % set(g,'Color','k');
% % set(g,'LineStyle','--');
% % g = vline((tInterp(idx)-r/c_inf)*1000);
% % set(g,'Color','k');
% % set(g,'LineStyle','--');
% 
% h = text(0.06, 5, 'peak pressure');
% set(h,'FontSize',24);
% h = text(1.1, 1.5, 'rise');
% set(h,'FontSize',24);
% h = text(1.1, 0.7, 'time');
% set(h,'FontSize',24);

disp('Driver OK')

%% Pressure pulse
function [pPres, R] = computePressure(bubble, sol, t, rho_inf, c_inf, r, aD)
    R = bubble(1,:); % bubble radius [m]
    V = 4/3*pi*R.^3; % bubble volume [m^3]
    U = bubble(2,:); % bubble wall velocity [m/s]
    mass = bubble(3,:); % bubble mass [kg]
    [~,solDY] = deval(sol, t); % Numerical differentiation
    % A = solDY(end-2,:); % acceleration
    A = solDY(end-4,:); % acceleration :: For 4 added states (2 shuttle, 3 port)
    % warning('Please change above line 43, it`s not right: accessing {q, bubble} data. Consider exporting from solution method');
    [tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation
    [tGhost, pGhost] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r + 2*aD); %ghost

    dt = 1e-6;
    tInterp = min(tDir):dt:max(tDir); % Interpolate direct, ghost waves to same uniform grid
    pDirInterp = pchip(tDir, pDir, tInterp);
    pGhostInterp = pchip(tGhost, pGhost, tInterp);
    pPres = pDirInterp - pGhostInterp; % High-to-low acoustic impedance: phase reversal
end