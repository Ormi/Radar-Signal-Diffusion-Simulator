%%
% @author xormos00
% @date Feb 2017
% @title Calculatin receiving frequency
% @input csv files
% Interpolation
% @return returning interpolated(180000) radar antenna characteristic
%
function res=parse_csv_file(src_file)
   filename = fullfile(src_file);
   vertical_data = readtable(filename);

   %%
   % Interpolation
   x = vertical_data.Var1;
   v = vertical_data.Var2;
   xq = -90:0.001:90;

   res = spline(x,v,xq);
end