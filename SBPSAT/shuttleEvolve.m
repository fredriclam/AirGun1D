function dz = shuttleEvolve(z, p_L, physConst)
% F = ma on shuttle, with isentropic compression of closed op. chamber
%   Input: z = [pos; vel]
%   Output: dz = [dpos; dvel] (time derivatives of state)

% Compute pressure based on (effective) length of operating chamber
% hard-coded version:
% p_R = physConst.p_R0 * ...
%     (physConst.operatingChamberLength/...
%       (physConst.operatingChamberLength-z(1)))^physConst.gamma;

% Fraction of op chamber length for pressure ramp up
rampUpFraction = 0.10;
% rampUpFraction = 5/6;

% TODO: Not complete!
p_R_front = physConst.p_R0 * ...
    (physConst.airgunOperatingChamberProfile(z(1)))^physConst.gamma;
p_R_rear = physConst.p_R0 * ...
    min([1,z(1)/(rampUpFraction*physConst.operatingChamberLength)]);
assert(p_R_front >= 0);

penaltyForce = physConst.shuttleBdryPenaltyStrength * ...
    (z(1) < 0) * abs(z(1));
% assert(penaltyForce == 0);
A_L = physConst.shuttle_area_left;
A_R = physConst.shuttle_area_right;

dz = [...(z(1) >= 0) * ... % This line disables motion when to the left
      z(2);
      1/physConst.shuttleAssemblyMass ...
      * (p_L*A_L - p_R_front*A_R ...                                            % Disables right-pressure
      + p_R_rear*physConst.shuttle_area_right_rear ...
      + penaltyForce ...
      ...+ 0 * linearDampingForce(z)*(z(1) <= rampUpFraction*physConst.operatingChamberLength) ...
      ...+ 0 * -5e2 * sign(z(2)) * 25 ...
      + constantDampingForce(z) ...
      ...+ 0*coulombDampingForce(z, physConst)...
      ) ];
assert(all(isreal(dz)));

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
