%%
% @author xormos00
% @date Feb 2017
%
% power_recei     % Pr % Received power [W]
% power_trans     % Pt % Transmitter peak power [W]
% trans_gain      % Ft % Transmitter gain [Hz]
% recei_gain      % Fr % Receiver gain [Hz]
% RCS             % RCS % Radar Cross Section - object reflection area [m^2]
% loss_factor     % loss % Loss of radar [dB]
% distance_r2p    % d % Range from the antenna to target [m]
%
% @return Power of receiving frequency in Watts
%
%           Pt * Fr * Ft * RCS
%   Pr = ------------------------
%          (4pi)^2 * d^4 * loss
%
	
function power_recei = radar_equation(trans_gain,recei_gain,RCS,loss_factor,distance)
	
	% Static value, anntenna
	power_trans = 5.6485e+03;

	power_recei = ((power_trans * trans_gain * recei_gain * RCS) / ...
				  (power(4*(pi), 2) * power(distance, 4) * loss_factor));
end