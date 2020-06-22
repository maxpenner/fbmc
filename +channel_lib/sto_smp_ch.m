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

function [samples_ch] = sto_smp_ch(obj, samples)

    if obj.sto >= 0
        samples_ch = [zeros(1,obj.sto) samples(1:end-obj.sto)];
    else
        samples_ch = [samples(1-obj.sto:end) zeros(1,-obj.sto)];
    end
end