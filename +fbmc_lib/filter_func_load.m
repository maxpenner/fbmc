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

function [filter_func] = filter_func_load(prototype_filter, n_subc, k)

    % For IOTA we need default values.
    % The alpha is coded into the prototype_filter name.
    tau0        = 0.5;
    nu0         = 1.0;
    filterSym   = 'odd';
    delay       = 0;
    alpha       = 2;
    
    % fuhrwerk passes the iota alpha in the string
    if strcmp(prototype_filter,'iota')
        prototype_filter = strcat('iota_', num2str(alpha));
    end

    % create a dummy instance of fuhrwerk's code
    dummy_instance = fbmc_lib.pulseDesign();
    
    % let the dummy instance load the filter function values
    [filter_func, ~, ~, ~, ~] = dummy_instance.generate(prototype_filter, n_subc, k, tau0, nu0, filterSym, delay);
end

