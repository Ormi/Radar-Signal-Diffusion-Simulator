%%
% @author xormos00
% @date Feb 2017
%
% power_recei     % Pr % Received power [W]
% power_trans     % Pt % Transmitter peak power [W]
% f_trans      % Ft % Transmitter gain [Hz]
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
	
function power_recei = radar_equation(f_trans,c,RCS,loss_factor,distance)
	
	% Static value, anntenna
	%power_trans = 5.6485e+03;
    
    %18dB
    %power_trans = 6.30957e-2;
    %13dB
    %f_trans = 1.9952e-2;
    %16dB
    %receiv_gain = 3.9810e-2;
	
	lambda = c/f_trans;

	power_recei = ((lambda * RCS * nthroot((10^loss_factor), 10)) / (power(4*(pi), 2) * power(distance, 4)));
				  
	%power_recei = ((power_trans * f_trans * receiv_gain * RCS * 10*log10(loss_factor)) / ...
				  %(power(4*(pi), 2) * power(distance, 4)));				  
end
