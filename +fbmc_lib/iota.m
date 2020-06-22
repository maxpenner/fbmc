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

function filter_func = iota(n_subc, k)
    % IOTA
    % ref: SIOHAN - 2000 - Cosine-modulated filterbanks based on extended Gaussian function
    flt_len = k*n_subc;
    alpha   = 1;
    nu_0    = 1/sqrt(2);
    tau_0   = 1/sqrt(2);
    big_K   = 14;
    t_sampl = 1/(n_subc*nu_0);
    limit   = 2*k/(4*nu_0);
    n       = -limit : t_sampl : limit; % time index

    sum1 = zeros(size(n));
    sum2 = zeros(size(sum1));
    for c = 0:big_K
        sum1   = sum1 + fbmc_lib.d_c_alpha_nu(c,alpha,nu_0, big_K)*[fbmc_lib.g_alpha(alpha,n+c/nu_0) + fbmc_lib.g_alpha(alpha,n-c/nu_0)];
        sum2   = sum2 + fbmc_lib.d_c_alpha_nu(c,1/alpha,tau_0, big_K)*[cos(2*pi*c*n./tau_0)];
    end
    % time domain prototype
    filter_func = 1/2*(sum1 .* sum2);
    filter_func = filter_func(1:flt_len) / sum(filter_func(1:flt_len));
    filter_func = filter_func/norm(filter_func);
    filter_func_f = fftshift(fft(filter_func));

end % func
