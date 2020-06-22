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

function [constell_pam_norm, constell_pam_norm_fac] = constellation_pam_norm(obj)

    % if M > 2 we are using both the real and imag part, so both parts need half the energy (0.5) of the average QAM symbol power (1)
    if obj.M > 2
        
        % how many real symbols do we have?
        M_real_symbols = sqrt(obj.M);
        
        % we need to split the energy into real and imaginary part
        target_average_power = 0.5;
        
    elseif obj.M == 2
        
        % we have to symbols
        M_real_symbols = 2;
        
        % we are using only the real part, so we can put the entire energy into it
        target_average_power = 1.0;
        
    else
        % Should never happen
        error('Unknown modulation type.');
    end
    
    % make sure it is an integer
    if floor(M_real_symbols) ~= M_real_symbols
        error('Non-integer number of real symbols');
    end

    % generate each symbol index once
    symbol_indices = 0 : 1 : M_real_symbols-1;

    % generate unnormalized symbols
    constellation_PAM_unnormalized = pammod(symbol_indices, M_real_symbols, 0, 'gray');

    % determine normalization factor
    constell_pam_norm_fac = sqrt( target_average_power * M_real_symbols / sum(abs(constellation_PAM_unnormalized).^2) );

    % create normalized constellation diagram
    constell_pam_norm = constellation_PAM_unnormalized*constell_pam_norm_fac;

    % SECURITY
    % check if power is actually 0.5
    difference = abs(sum(abs(constell_pam_norm).^2)/M_real_symbols - target_average_power);
    if difference > 1e-7
        error('PAM constellation does not have the correct average power.');
    end    
end