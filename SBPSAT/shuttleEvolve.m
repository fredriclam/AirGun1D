function dz = shuttleEvolve(z, p_L, physConst)
% F = ma on shuttle, with isentropic compression of closed op. chamber
%   Input: z = [pos; vel]
%   Output: dz = [dpos; dvel] (time derivatives of state)

% Compute pressure based on (effective) length of operating chamber
% hard-coded version:
% p_R = physConst.p_R0 * ...
%     (physConst.operatingChamberLength/...
%       (physConst.operatingChamberLength-z(1)))^physConst.gamma;

p_R = physConst.p_R0 * ...
    (physConst.airgunOperatingChamberProfile(z(1)))^physConst.gamma;

penaltyForce = physConst.shuttleBdryPenaltyStrength * ...
    (z(1) < 0) * abs(z(1));
% assert(penaltyForce == 0);
A_L = physConst.shuttle_area_left;
A_R = physConst.shuttle_area_right;

dz = [...(z(1) >= 0) * ... % This line disables motion when to the left
      z(2);
      1/physConst.shuttleAssemblyMass ...
      * (p_L*A_L - p_R*A_R + penaltyForce ...
      + linearDampingForce(z) )];
assert(all(isreal(dz)));

% Arbitrary damping model
% Input: state vector z --- [pos; vel]
function dampingForce = linearDampingForce(z)
dampingForce = - 1e4 * z(2);

% Viscoelastic model: is this just linear damping?
% Input: state vector z --- [pos; vel]
function dampingForce = viscoelasticDampingForce(z)
dampingForce = - 1e3 * z(2);
