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

function filter_func = phydyas(n_subc, k)
    % phydyas prototype
    flt_len = k*n_subc;
    switch k
        case 4
            A       = 0.971960;
            h_coef  = [1 A sqrt(2)/2 sqrt(1-A^2)];
            filter_func_f   = [h_coef zeros(1,flt_len-(2*k-1)) h_coef(k:-1:2)];
            filter_func     = ifftshift(ifft(filter_func_f));
            filter_func     = filter_func/norm(filter_func);
            
        case 3
            A       = 0.971960;
            h_coef  = [1 A sqrt(1-A^2)];
            filter_func_f   = [h_coef zeros(1,flt_len-(2*k-1)) h_coef(k:-1:2)];
            filter_func     = ifftshift(ifft(filter_func_f));
            filter_func     = filter_func/norm(filter_func);       
    end
end % func
