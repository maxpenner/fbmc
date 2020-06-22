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

function [weights] = insertion_aux_solver(interference_data, interference_pilot_indiv)

    Q = interference_data;
    p = interference_pilot_indiv;

    % A: Moore-Penrose-Pseudoinverse
    weights = pinv(p*p')*(-Q*p);
    
    % B: Poschadel-Penner-PI
    %weights = -Q/(p'*p)*p;
    
    % SECURITY
    counter_Q = sum(weights.*interference_pilot_indiv);
    if abs(1-abs(Q/counter_Q)) > 1e-6
        error('Q and counter Q are too far apart.');
    end
end