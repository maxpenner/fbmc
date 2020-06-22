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

function [oqam_frame] = demodulation(mod_param)

    % choose which prototype filter to use
    %switch mod_param.prototype_filter
    %    case 'iota'
    %        filter_func = fbmc_lib.iota(mod_param.n_subc, mod_param.k);
    %    case 'phydyas'
    %        filter_func = fbmc_lib.phydyas(mod_param.n_subc, mod_param.k);
    %end
    filter_func = mod_param.filter_func;

    % splitting input into two branches
    branch_1_sampl = mod_param.samples_rx(1:end-mod_param.n_subc);
    branch_2_sampl = mod_param.samples_rx(mod_param.n_subc/2+1:end-mod_param.n_subc/2);

    % serial-to-parallel conversion
    branch_1_sampl = reshape(branch_1_sampl,mod_param.n_subc,numel(branch_1_sampl)/mod_param.n_subc);
    branch_2_sampl = reshape(branch_2_sampl,mod_param.n_subc,numel(branch_2_sampl)/mod_param.n_subc);

    % filter bank receiver
    branch_1_ppn = fbmc_lib.ppn(mod_param.n_subc, fliplr(filter_func),mod_param.k);
    branch_2_ppn = fbmc_lib.ppn(mod_param.n_subc, fliplr(filter_func),mod_param.k);

    % expanding the input of the convolution to avoid erros due to cyclic convolution
    branch_1_ppn     = flipud([branch_1_ppn zeros(mod_param.n_subc,size(branch_1_sampl,2)-1)]);
    branch_2_ppn     = flipud([branch_2_ppn zeros(mod_param.n_subc,size(branch_2_sampl,2)-1)]);
    branch_1_sampl   = [branch_1_sampl zeros(mod_param.n_subc,mod_param.k-1)];
    branch_2_sampl   = [branch_2_sampl zeros(mod_param.n_subc,mod_param.k-1)];

    % application of the PPN filter to the FFT input
    branch_1_sampl = ifft(fft(branch_1_ppn,[],2).*fft(branch_1_sampl,[],2),[],2);
    branch_2_sampl = ifft(fft(branch_2_ppn,[],2).*fft(branch_2_sampl,[],2),[],2);

    % fft
    branch_1_symb_fft = fft(branch_1_sampl);
    branch_2_symb_fft = fft(branch_2_sampl);

    % removing the transient received symbols at the beginning (k-1) and the end (k-1)
    %oqam_frame            = zeros(mod_param.n_subc,2*mod_param.n_symbols);
    %oqam_frame(:,1:2:end) = branch_1_symb_fft(:,mod_param.k+(0:mod_param.n_symbols-1));
    %oqam_frame(:,2:2:end) = branch_2_symb_fft(:,mod_param.k+(0:mod_param.n_symbols-1));
    oqam_frame_shifted            = zeros(mod_param.n_subc,2*mod_param.n_symbols);
    oqam_frame_shifted(:,1:2:end) = branch_1_symb_fft(:,mod_param.k+(0:mod_param.n_symbols-1));
    oqam_frame_shifted(:,2:2:end) = branch_2_symb_fft(:,mod_param.k+(0:mod_param.n_symbols-1));
    
    % Currently the oqam frame is shifted and DC is at index 1. This is due to matlab's ifft.
    % Remove this shift.
    oqam_frame = fftshift(oqam_frame_shifted, 1);

end
