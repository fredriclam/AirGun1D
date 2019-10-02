% Launch script for testing models with diff from [LW2019]

clear;
clc;
addpath .\FlowRelations

%% Set parameters
nx = 100;       % Number of grid points per 1 m of air gun length

% % Reference test parameters
% r = 75;         % Distance from source to receiver [m]
% c_inf = 1482;   % Speed of sound in water [m/s]
% rho_inf = 1000; % Density of water [kg/m^3]
% airgunPressure = 2000;     % [psi]
% airgunLength = 0.6;        % [m]
% airgunPortArea = 72;       % [in^2]
% airgunCrossSecArea = 16;   % [in^2]
% airgunDepth = 7.5;         % [m]
% 
% % Compression factor as function of shuttle position
% airgunFiringChamberProfile = @(a) error('Placeholder; not implemented');
% accelerationLength = 2*0.0254; % [m]
% airCushionLength = 1*0.0254; % [m]
% % Compression factor as function of shuttle position
% airgunOperatingChamberProfile = @(xi) (xi - accelerationLength < 0) * 1 ...
%     + (xi - accelerationLength > 0) * ...
%     (airCushionLength / (airCushionLength - (xi - accelerationLength)));
% bubbleInitialVolume = 600; % [cui]
% shuttleBdryPenaltyStrength = 1e11; % [N/m]

% % New set for data matching -- ballpark numbers
r = 3;                                            % Distance from source to receiver [m]
c_inf = 1482;                                     % Speed of sound in water [m/s]
rho_inf = 1000;                                   % Density of water [kg/m^3]
airgunPressure = 1000;                            % [psi]
airgunLength = 20800 / (pi*10.020^2/4) * 0.0254;  % [m] For 20800 cui
airgunPortArea = 2.5 * 0.50 * pi * 11;            % [in^2] % OD; covers 0.5 of cyl
airgunCrossSecArea = pi*10.020^2/4;               % [in^2]
airgunDepth = 7.5;                                % [m]
bubbleInitialVolume = 600; % [cui]
shuttleBdryPenaltyStrength = 1e11; % [N/m]
% Compression factor as function of shuttle position
airgunFiringChamberProfile = @(a) error('Placeholder; not implemented');
accelerationLength = (3.009-0.542)*0.0254; % [m]
airCushionLength = 0.542*0.0254;           % [m]
% Compression factor as function of shuttle position
airgunOperatingChamberProfile = @(xi) (xi - accelerationLength < 0) * 1 ...
    + (xi - accelerationLength > 0) * ...
    (airCushionLength / (airCushionLength - (xi - accelerationLength)));

% Run solve for both models
[sol, q, bubble, shuttle, plug, ...
    q2, bubble2, shuttle2, plug2, monitorStates] = runEulerCodeShuttleDual(nx, ...
    airgunPressure, airgunLength, airgunCrossSecArea, ...
    airgunPortArea, airgunDepth, airgunFiringChamberProfile, ...
    airgunOperatingChamberProfile, bubbleInitialVolume, ...
    shuttleBdryPenaltyStrength);

t = sol.x; % time

%% Post processing
gamma = 1.4;

[~,solDY] = deval(sol, t); % Numerical differentiation (second output arg)
% Cut into the two halves
DY1 = solDY(1:size(solDY,1)/2,:);
DY2 = solDY(size(solDY,1)/2+1:end,:);
[pPres, R, tInterp] = computePressure(bubble, DY1, t, rho_inf, c_inf, r, airgunDepth, size(q,1));
[pPres2, R2, tInterp2] = computePressure(bubble2, DY2, t, rho_inf, c_inf, r, airgunDepth, size(q,1));

%% Plot 1: density
figure(1001); clf;
tmax = 30;

subplot(1,2,1);
rho = q(1:3:end,:);
xAxis = linspace(0,1.2,size(rho, 1))';
tAxis = 1000*sol.x;
[xx, tt] = meshgrid(xAxis, tAxis);
contourf(xx, tt, rho', 'LineStyle', 'none');

subplot(1,2,2);
rho2 = q2(1:3:end,:);
xAxis = linspace(0,1.2,size(rho2, 1))';
tAxis = 1000*sol.x;
[xx, tt] = meshgrid(xAxis, tAxis);
contourf(xx, tt, rho2', 'LineStyle', 'none', 'levelstep', 2);

hold on;
plot(shuttle2(1,:)+1, tAxis, 'g');

for i = 1:2
    subplot(1,2,i);
    title('rho [kg / m^3]');
    ylim([0, tmax]);
    caxis([0, 160]);
    colorbar;
    xlabel('Position [m]');
    ylabel('Time [ms]');
end

windowPos = get(gcf,'position');
set(gcf,'position',[windowPos(1), windowPos(2), ...
    960, 350]);
colormap bone;


%% Plot 2: pressure
figure(1002); clf;
tmax = 30;

subplot(1,3,1);
p = (gamma-1)*( q(3:3:end,:)-1/2*q(2:3:end,:).^2./q(1:3:end,:) );
xAxis = linspace(0,1.2,size(p, 1))';
tAxis = 1000*sol.x;
[xx, tt] = meshgrid(xAxis, tAxis);
contourf(xx, tt, p', 'LineStyle', 'none');

subplot(1,3,2);
p2 = (gamma-1)*( q2(3:3:end,:)-1/2*q2(2:3:end,:).^2./q2(1:3:end,:) );
xAxis = linspace(0,1.2,size(p2, 1))';
tAxis = 1000*sol.x;
[xx, tt] = meshgrid(xAxis, tAxis);
contourf(xx, tt, p2', 'LineStyle', 'none', 'LevelStep', 0.1e6);

hold on;
plot(shuttle2(1,:)+1, tAxis, 'g');

for i = 1:2
    subplot(1,3,i);
    title('p [Pa]');
    ylim([0, tmax]);
    caxis([0, 14e6]);
    colorbar;
    xlabel('Position [m]');
    ylabel('Time [ms]');
end

windowPos = get(gcf,'position');
set(gcf,'position',[windowPos(1), windowPos(2), ...
    1080, 350]);
colormap bone;

subplot(1,3,3);
plot(shuttle2(1,:), tAxis, 'k', 'LineWidth', 1);
xlabel('Shuttle position [m]');
ylabel('Time [ms]');


%% Boundary velocity record and shuttle dynamics
% Compute first time the shuttle closes
figure(1003); clf;

u2 = q2(2:3:end,:) ./ q2(1:3:end,:);
u_R2 = u2(end,:);

subplot(2,1,1);
plot(shuttle2(1,:), tAxis, 'LineWidth', 1);
title('Shuttle pos [m]')
xlabel('xi [m]')
ylabel('t [ms]')
hold on
plot(0.0762 *[1,1], [0, tAxis(end)], '-');
plot(accelerationLength*[1,1], [0, tAxis(end)], '--');

subplot(2,1,2);
plot(shuttle2(2,:), tAxis, 'LineWidth', 1);
title('Shuttle velocity [m/s]')
hold on
plot(u_R2, tAxis, 'k', 'LineWidth', 1);
xlabel('xi-dot or u [m/s]')
ylabel('t [ms]')
legend({'Shuttle','u_R'});

%% Bubble characteristics
figure(1004); clf;
figPos = get(gcf,'Position');
set(gcf,'Position',[figPos(1) 0 600 900]);

% bubble volume
V = 4/3*pi*R.^3;
subplot(3,1,1);
plot(t*1000,R, 'LineWidth', 1);
xlabel('t [ms]');
ylabel('Radius [m]');
hold on
plot(t*1000,R2,'k','LineWidth', 1);

subplot(3,1,2);
plot(t*1000,V, 'LineWidth', 1);
xlabel('t [ms]');
ylabel('Volume [m^3]');

hold on;
V2 = 4/3*pi*R2.^3;
plot(t*1000, V2, 'k', 'LineWidth', 1);
legend({'Instant open 10 ms','Shuttle'},'location','best')

% acoustic pressure
subplot(3,1,3);
h = plot((tInterp-r/c_inf)*1000, pPres*1e-5*r, 'LineWidth', 1);
ylabel('\Delta p (bar m)');
xlabel('Time, t-r/c_\infty (ms)');
hold on;
plot((tInterp2-r/c_inf)*1000, pPres2*1e-5*r, 'k', 'LineWidth', 1);



%% Down-bore temperature msmt comparison
figure(1005); clf;
cv = 718;
T = (q(3:3:end,:) - ...
    0.5 * q(2:3:end,:).^2 ./ q(1:3:end,:)) ./ q(1:3:end,:) / cv;
T2 = (q2(3:3:end,:) - ...
    0.5 * q2(2:3:end,:).^2 ./ q2(1:3:end,:)) ./ q2(1:3:end,:) / cv;

T_L = T(30,:);
T_L2 = T2(30,:);

subplot(1,3,1);
plot(tAxis, T_L);

subplot(1,3,2);
plot(tAxis, T_L2);

for i = 1:2
    subplot(1,3,i);
    title('Closed-end temperature');
%     xlim([0, tmax]);
%     caxis([0, 14e6]);
    xlabel('Time [ms]');
    ylabel('T [K]');
    xlim([0, 80]);
    ylim([120, 300]);
end

figPos = get(gcf, 'position');
set(gcf, 'position', [figPos(1:2), 1300, 420])

%% Attempt to load data
% Should not work on remote devices
subplot(1,3,3);

try
    load ../../LinuxShare/HiTest_Data/HiTestData_v1.mat;
    T_exper_K = 5/9*(HiTestData(24).iNetCh1Data - 32) + 273.15;
    plot(1000*(HiTestData(24).iNetTimeAxisT(3710:3790) - ...
        HiTestData(24).iNetTimeAxisT(3710)), ...
        T_exper_K(3710:3790)); % [K] vs [s]
    xlim([0, 80]);
    ylim([120, 300]);
catch
    disp('Failed to locate data. Ignoring exp. data plot')
end
clear HiTestData;
title('Closed-end temperature data');

%% Phase plot
figure(1006); clf;

subplot(1,2,1);
subsonicStates =       monitorStates(:,monitorStates(3,:)==1);
chamberLimitedStates = monitorStates(:,monitorStates(3,:)==2);
shockStates =          monitorStates(:,monitorStates(3,:)==3);
portLimitedStates =    monitorStates(:,monitorStates(3,:)==4);

plot(subsonicStates(1,:), subsonicStates(2,:), 'b.'); hold on
plot(chamberLimitedStates(1,:), chamberLimitedStates(2,:), 'g.');
plot(shockStates(1,:), shockStates(2,:), 'k.');
plot(portLimitedStates(1,:), portLimitedStates(2,:), 'r.');

xlabel('$M_\mathrm{a}$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$M_\mathrm{port}$', 'Interpreter', 'latex', 'FontSize', 14)
yL = ylim;
ylim([0, yL(2)]);
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');

subplot(1,2,2);
plot(subsonicStates(8,:), subsonicStates(7,:), 'b.'); hold on
plot(chamberLimitedStates(8,:), chamberLimitedStates(7,:), 'g.');
plot(shockStates(8,:), shockStates(7,:), 'k.');
plot(portLimitedStates(8,:), portLimitedStates(7,:), 'r.');
plot(monitorStates(8,:),monitorStates(7,:),'-');

xlabel('$A_\mathrm{p} / A_\mathrm{cs}$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$p_\mathrm{t} / p^* $', 'Interpreter', 'latex', 'FontSize', 14)
yL = ylim;
ylim([0, yL(2)]);
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');
set(gca,'YScale','log')
ylim([1, 1e3])

% %% Plot bubble radius
% figure(1); clf;
% tAxis = 1000*sol.x;
% plot(tAxis, R);
% xlabel('Time [ms]');
% ylabel('Bubble radius [m]');
% hold on
% plot(tAxis, R2);

% %% Chamber density contours
% figure(101); clf;
% rho = q(1:3:end,:);
% xAxis = linspace(0,1.2,size(rho, 1))';
% tAxis = 1000*sol.x;
% [xx, tt] = meshgrid(xAxis, tAxis);
% contourf(xx, tt, rho', 'LineStyle', 'none');
% xlabel('Position [m]')
% ylabel('Time [ms]')
% caxis([0, 160]);
% colorbar;
% 
% windowPos = get(gcf,'position');
% set(gcf,'position',[windowPos(1), windowPos(2), ...
%     560, 280]);
% colormap bone;
% 
% title('rho [kg / m^3]')
% ylim([0, 10])
% 
% %% Chamber particle velocity contours
% figure(102); clf;
% u = q(2:3:end,:) ./ q(1:3:end,:);
% xAxis = linspace(0,1.2,size(rho, 1))';
% tAxis = 1000*sol.x;
% [xx, tt] = meshgrid(xAxis, tAxis);
% contourf(xx, tt, u', 'LineStyle', 'none', 'LevelStep', 20);
% xlabel('Position [m]')
% ylabel('Time [ms]')
% caxis([0, 350]);
% colorbar;
% windowPos = get(gcf,'position');
% set(gcf,'position',[windowPos(1), windowPos(2), ...
%     560, 280]);
% colormap bone;
% 
% title('u [m/s]')
% ylim([0, 10])

% %% Boundary velocity record and shuttle dynamics
% % Compute first time the shuttle closes
% figure(2); clf;
% subplot(3,1,1);
% u_R = u(end,:);
% 
% u2 = q2(2:3:end,:) ./ q2(1:3:end,:);
% u_R2 = u2(end,:);
% 
% plot(tAxis, u_R, 'LineWidth', 1);
% title('Boundary velocity [m]')
% hold on
% yL = ylim;
% ylim(yL)
% 
% % Superimpose shuttle-limited data
% closingIndex = find(shuttle2(1,:) <= 0, 1);
% if isempty(closingIndex)
%     closingIndex = 1;
% end
% hold on
% plot(tAxis, u_R2, 'k', 'LineWidth', 1);
% yL = ylim;
% plot(tAxis(closingIndex)*[1,1], [-1000 10000], 'k--')
% ylim(yL);
% legend({'Instant open 10 ms','Shuttle'});
% 
% subplot(3,1,2);
% plot(tAxis, shuttle2(1,:), 'LineWidth', 1);
% title('Virtual shuttle pos [m]')
% hold on
% yL = ylim;
% plot(tAxis(closingIndex)*[1,1], [-1e4 1e4], 'k--')
% ylim(yL);
% 
% subplot(3,1,3);
% plot(tAxis, shuttle2(2,:), 'LineWidth', 1);
% title('Virtual shuttle velocity [m]')
% hold on
% yL = ylim;
% plot(tAxis(closingIndex)*[1,1], [-1e4 1e4], 'k--')
% ylim(yL)

% %%
% figHand1 = figure(3); clf;
% figPos = get(figHand1,'Position');
% set(figHand1,'Position',[figPos(1) figPos(2) 600 900]);
% 
% %% bubble volume
% V = 4/3*pi*R.^3;
% subplot(2,1,1);
% plot(t*1000,V, 'LineWidth', 1);
% % xmax = 500;
% % xlim([0 xmax]);
% % hold on;
% xlabel('Time, t (ms)');
% ylabel('Volume (m^3)');
% % ylim([0 1.5])
% 
% hold on;
% V2 = 4/3*pi*R2.^3;
% plot(t*1000, V2, 'k', 'LineWidth', 1);
% legend({'Instant open 10 ms','Shuttle'},'location','best')
% 
% %% acoustic pressure
% subplot(2,1,2);
% h = plot((tInterp-r/c_inf)*1000, pPres*1e-5*r, 'LineWidth', 1);
% % xlim([0 xmax]);
% % set(h.Parent,'YTick',[-4 0 4]);
% ylabel('\Delta p (bar m)');
% xlabel('Time, t-r/c_\infty (ms)');
% % ylim([-4 6])
% 
% hold on;
% plot((tInterp2-r/c_inf)*1000, pPres2*1e-5*r, 'k', 'LineWidth', 1);

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

%% Temperature tracking
if false
    c_v = 718;
    T =((q(3:3:end,:) - 0.5 * rho .* u.^2)/c_v);
    T2 =((q2(3:3:end,:) - 0.5 * q2(1:3:end,:) .* u2.^2)/c_v);
end


%% Console report
disp('Launch script finished.')

folderNameRoot = "testLaunch";
index = 1;
while exist(folderNameRoot + sprintf('%02d',index),'dir')
    index = index + 1;
end
folderName = folderNameRoot + sprintf('%02d',index);
mkdir(folderName);
cd(folderName);
for i = 1001:1006
    figure(i);
    savefig(num2str(i));
end

save('runningData','monitorStates')

cd ..

disp("Figures saved to "+ folderName);

%% Define pressure pulse at a distance
function [pPres, R, tInterp] = computePressure(bubble, DY, t, rho_inf, c_inf, r, airgunDepth, qLength)
    R = bubble(1,:); % bubble radius [m]
    V = 4/3*pi*R.^3; % bubble volume [m^3]
    U = bubble(2,:); % bubble wall velocity [m/s]
    mass = bubble(3,:); % bubble mass [kg]
    
    % A = solDY(end-2,:); % acceleration
    A = DY(qLength+2,:); % Seek acceleration entry (located after q, and R)
    % warning('Please change above line 43, it`s not right: accessing {q, bubble} data. Consider exporting from solution method');
    [tDir, pDir] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r); % direct pressure perturbation
    [tGhost, pGhost] = pressure_eqn(t', R', U', A', rho_inf, c_inf, r + 2*airgunDepth); %ghost

    dt = 1e-6;
    tInterp = min(tDir):dt:max(tDir); % Interpolate direct, ghost waves to same uniform grid
    pDirInterp = pchip(tDir, pDir, tInterp);
    pGhostInterp = pchip(tGhost, pGhost, tInterp);
    pPres = pDirInterp - pGhostInterp; % High-to-low acoustic impedance: phase reversal
end