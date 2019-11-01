function [dz, p_rear, p_front] = shuttleEvolve(z, p_L, physConst, opChamber)
% F = ma on shuttle, with isentropic compression of closed op. chamber
%   Input:
%     z = [pos; vel; m_L; E_L; m_R; E_R]
%     p_L = pressure on piston left flange
%     physConst = struct of physical constnats
%     opChamber = object of type OperatingChamber
%   Output: dz = [dpos; dvel; ...] (time derivatives of state z)

% Compute pressure based on (effective) length of operating chamber
% hard-coded version:
% p_R = physConst.p_R0 * ...
%     (physConst.operatingChamberLength/...
%       (physConst.operatingChamberLength-z(1)))^physConst.gamma;

% Unpack state:
m_rear = z(3);
E_rear = z(4);
m_front = z(5);
E_front = z(6);

% % Fraction of op chamber length for pressure ramp up
% rampUpFraction = 0.10;
% % rampUpFraction = 5/6;
% 
% % TODO: Not complete!
% p_R_front = physConst.p_R0 * ...
%     (physConst.airgunOperatingChamberProfile(z(1)))^physConst.gamma;
% p_R_rear = physConst.p_R0 * ...
%     min([1,z(1)/(rampUpFraction*physConst.operatingChamberLength)]);
% assert(p_R_front >= 0);

% Compute elastic penalty for constraint z(1) >= 0
penaltyForce = physConst.shuttleBdryPenaltyStrength * ...
    (z(1) < 0) * abs(z(1));
A_L = physConst.shuttle_area_left;
A_R = physConst.shuttle_area_right;

% Compute mass flow rates
% Geometrically constrained density
rho_rear = m_rear / opChamber.rearVolume(z(1));
% Energy-constrained temperature
T_rear = E_rear / (m_rear * physConst.c_v);
% Ideal gas pressure
p_rear = rho_rear * physConst.Q * T_rear;
% Geometrically constrained density
rho_front = m_front / opChamber.frontVolume(z(1));
% Energy-constrained temperature
T_front = E_front / (m_front * physConst.c_v);
% Ideal gas pressure
p_front = rho_front * physConst.Q * T_front;

% Define local aliases
g = physConst.gamma;
c_p = physConst.c_v + physConst.Q;

% Compute direction of flow and pressure ratio
flowSignL2R = sign(p_rear - p_front);
pRatio = min([p_rear/p_front, p_front/p_rear]);
% Compute mach number, capping to choked
M = min([1, machPressureFunction(g, pRatio)]);
% Compute effective stagnation properties
if flowSignL2R >= 0
    pMax = p_rear;
    rhoMax = rho_rear;
    TMax = T_rear;
else
    pMax = p_front;
    rhoMax = rho_front;
    TMax = T_front;
end

efficiencyFactor = 1.0;

flowL2R = efficiencyFactor * ...
    flowSignL2R * opChamber.gapArea(z(1)) * M * ...
    sqrt(g * pMax * rhoMax) * ...
    (1 + (g-1)/2 * M^2)^(0.5*(-g-1)/(g-1));

dz = [...(z(1) >= 0) * ... % This line disables motion when to the left
      z(2);
      1/physConst.shuttleAssemblyMass ...
      * (p_L*A_L - p_front*A_R ...
      + p_rear*physConst.shuttle_area_right_rear ...
      + penaltyForce ...
      ...+ 0 * linearDampingForce(z)*(z(1) <= rampUpFraction*physConst.operatingChamberLength) ...
      ...+ 0 * -5e2 * sign(z(2)) * 25 ...
      ...+ constantDampingForce(z) ...
      ...+ 0*coulombDampingForce(z, physConst)...
      );
      -flowL2R;
      -c_p * flowL2R * TMax - p_rear * physConst.shuttle_area_right_rear * z(2); % Energy is approximatied...
      flowL2R;
      c_p * flowL2R * TMax + p_front * A_R * z(2);];
assert(all(isreal(dz)) && all(~isnan(dz)));

% Arbitrary damping model
% Input: state vector z --- [pos; vel]
function dampingForce = linearDampingForce(z)
% dampingForce = -4e4 * z(2);% - 5e4 * z(2);
% dampingForce = -3e4 * z(2);
% dampingForce = -1e4 * z(2);
% dampingForce = -5e3 * z(2);
% dampingForce = -2e3 * z(2);
dampingForce = -5e2 * z(2);

% Constant damping force
% Coefficient of friction (would be empirical). Doesn't really do much!
% Input: state vector z --- [pos; vel]
function dampingForce = coulombDampingForce(z, physConst)
magnitude = 0.3 * 9.8 * physConst.shuttleAssemblyMass;
dampingForce = -sign(z(2)) * magnitude;

function dampingForce = constantDampingForce(z)
dampingForce = -10e3; % 2 kN (overpredicting the INcompressible turbulent 1/7-law)

% Viscoelastic model: is this just linear damping?
% Input: state vector z --- [pos; vel]
function dampingForce = viscoelasticDampingForce(z)
dampingForce = - 1e3 * z(2);
