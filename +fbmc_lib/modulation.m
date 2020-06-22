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

function [samples] = modulation(mod_param)

    % Create a new temporary shifted version so that DC carrier has index 1.
    % This corresponds to matlab's fft or ifft which assumes DC at index 1.
    oqam_frame_shifted = fftshift(mod_param.oqam_frame,1);

    % choose which prototype filter to use
%     switch mod_param.prototype_filter
%         case 'iota'
%             filter_func = fbmc_lib.iota(mod_param.n_subc, mod_param.k);
%         case 'phydyas'
%             filter_func = fbmc_lib.phydyas(mod_param.n_subc, mod_param.k);
%     end
    filter_func = mod_param.filter_func;

    % main

    % ifft
    %branch_1_symb_ifft = ifft(mod_param.oqam_frame(:,1:2:end))*mod_param.n_subc;
    %branch_2_symb_ifft = ifft(mod_param.oqam_frame(:,2:2:end))*mod_param.n_subc;    
    branch_1_symb_ifft = ifft(oqam_frame_shifted(:,1:2:end))*mod_param.n_subc;
    branch_2_symb_ifft = ifft(oqam_frame_shifted(:,2:2:end))*mod_param.n_subc;

    % PPN (poly phase network)
    branch_1_ppn = fbmc_lib.ppn(mod_param.n_subc, filter_func,mod_param.k);
    branch_2_ppn = fbmc_lib.ppn(mod_param.n_subc, filter_func,mod_param.k);

    % expanding the input of the convolution to avoid erros due to cyclic convolution
    branch_1_ppn    = [branch_1_ppn  zeros(mod_param.n_subc,size(branch_1_symb_ifft,2)-1)];
    branch_2_ppn    = [branch_2_ppn  zeros(mod_param.n_subc,size(branch_2_symb_ifft,2)-1)];
    branch_1_symb_ifft   = [branch_1_symb_ifft zeros(mod_param.n_subc,mod_param.k-1)];
    branch_2_symb_ifft   = [branch_2_symb_ifft zeros(mod_param.n_subc,mod_param.k-1)];

    % application of the PPN filter to the IFFT output
    branch_1_sampl = ifft(fft(branch_1_ppn,[],2).*fft(branch_1_symb_ifft,[],2),[],2);
    branch_2_sampl = ifft(fft(branch_2_ppn,[],2).*fft(branch_2_symb_ifft,[],2),[],2);

    % parallel-to-serial conversion
    branch_1_serial = reshape(branch_1_sampl,1,numel(branch_1_sampl));
    branch_2_serial = reshape(branch_2_sampl,1,numel(branch_2_sampl));

    % delaying branch #2 about n_subc/2 samples
    branch_1_serial = [branch_1_serial zeros(1,mod_param.n_subc)];
    branch_2_serial = [zeros(1,mod_param.n_subc/2) branch_2_serial zeros(1,mod_param.n_subc/2)];

    % adding the two branches together
    samples = branch_1_serial + branch_2_serial;

end % func