%%
% @author xormos00
% @date Feb 2017
% @title Calculatin receiving frequency
% @input speed of object v_object
% @input transmiting frequency F_trans
% Speed of light lightspeed
% @return receiving frequency
%
%               Ft
% Fr = 2 * v * ---- * cos(Alfa)
%               c
%
function res = return_signal_freq(v_object, F_trans, angle)
	c = physconst('Lightspeed');
	res = 2 * v_object * ((F_trans) / c) * cos(angle);
end