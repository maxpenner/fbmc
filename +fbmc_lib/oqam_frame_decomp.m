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

function [frame, oqam_frame_equ_no_pattern_rx] = oqam_frame_decomp(mod_param)
    
    % real - imag structure matrix
    re_im_mat_1 = conj(fbmc_lib.re_im_pattern(mod_param.n_subc,size(mod_param.oqam_frame_equ_rx,2)/2,1));
    re_im_mat_2 = conj(fbmc_lib.re_im_pattern(mod_param.n_subc,size(mod_param.oqam_frame_equ_rx,2)/2,2));

    % inverting re-im matrix 
    branch_1_symb = mod_param.oqam_frame_equ_rx(:,1:2:end) .* re_im_mat_1;
    branch_2_symb = mod_param.oqam_frame_equ_rx(:,2:2:end) .* re_im_mat_2;

    % reconstructing complex QAM symbols by taking real part
    frame = real(branch_1_symb) + 1i*real(branch_2_symb);
    
    % FeelMaTyc
    oqam_frame_equ_no_pattern_rx = zeros(size(mod_param.oqam_frame_equ_rx));
    oqam_frame_equ_no_pattern_rx(:, 1:2:end) = mod_param.oqam_frame_equ_rx(:,1:2:end) .* re_im_mat_1;
    oqam_frame_equ_no_pattern_rx(:, 2:2:end) = mod_param.oqam_frame_equ_rx(:,2:2:end) .* re_im_mat_2;
end
