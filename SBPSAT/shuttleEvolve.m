function dz = shuttleEvolve(z, p_L, physConst)
% F = ma on shuttle, with isentropic compression of closed op. chamber
%   Input: z = [pos; vel]
%   Output: dz = [dpos; dvel] (time derivatives of state)

% Compute pressure based on (effective) length of operating chamber
p_R = physConst.p_R0 * ...
    (physConst.operating_chamber_length/...
      (physConst.operating_chamber_length-z(1)))^physConst.gamma;
penaltyForce = physConst.shuttleBdryPenaltyStrength * ...
    (z(1) < 0) * abs(z(1));
% assert(penaltyForce == 0);
A_L = physConst.shuttle_area_left;
A_R = physConst.shuttle_area_right;

dz = [...(z(1) >= 0) * ... % This line disables motion when to the left
      z(2);
      1/physConst.shuttleAssemblyMass ...
      * (p_L*A_L - p_R*A_R + penaltyForce ...
      - 0 * 1e3 * z(2) )];
      ...- 1e-1 * z(2) * abs(z(1) < 0))]; % arbitrary linear damping