%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

clear;
close all;
%profile on

% create an instance to generate the filter function
target = fbmc_chain();

% calculate the ambiguity inteferente with higher resolution
k_k_dot_diff_max = 3;
m_m_dot_diff_max = 5;
f_res = 50;
t_res = 64/2;                 % must be a devider of target.n_subc/2
[ambiguity_interference, ...
    k_k_dot_diffs, ...
    m_m_dot_diffs] = fbmc_lib.ambiguity_interference(target, k_k_dot_diff_max, m_m_dot_diff_max, f_res, t_res);

% take the logarithm of the absolute
ambiguity_interference_abs = abs(ambiguity_interference);
ambiguity_interference_abs_log = 20*log10(ambiguity_interference_abs);

% limit to a minimum
minimum = -100;
ambiguity_interference_abs_log(ambiguity_interference_abs_log < minimum) = minimum;

% plot 2d picture, we only need one vector
h = pcolor(m_m_dot_diffs(1,:), k_k_dot_diffs(:,1), ambiguity_interference_abs_log);
set(h, 'EdgeColor', 'none');
colorbar()
colormap(jet)

% labels
title('Ambiguity Function')
xlabel('\tau / T/2')
ylabel('\nu / carrier spacing 1/K')

%profile off
%profile viewer
disp('Simulation finished')