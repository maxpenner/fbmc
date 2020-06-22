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

function [oqam_frame] = oqam_frame_comp(mod_param)

    % real - imag structure matrix
    re_im_mat_1 = fbmc_lib.re_im_pattern(mod_param.n_subc,mod_param.n_symbols,1);
    re_im_mat_2 = fbmc_lib.re_im_pattern(mod_param.n_subc,mod_param.n_symbols,2);

    % real- imag -interleaved symbol
    branch_1_symb = re_im_mat_1.*real(mod_param.frame);
    branch_2_symb = re_im_mat_2.*imag(mod_param.frame);
    
    % interleaving the two branches
    oqam_frame            = zeros(mod_param.n_subc,2*mod_param.n_symbols);
    oqam_frame(:,1:2:end) = branch_1_symb;
    oqam_frame(:,2:2:end) = branch_2_symb; 
end