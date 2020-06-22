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

function afb_pilot = afb(fsync_param, samples_rx)

    % Number of oqam_frame vectors without transient symbols in received samples for real OR imaginary part.
    % The final number of vectors is doubled (real and imag).
    n_subc_multiples = numel(samples_rx)/fsync_param.n_subc - 2*(fsync_param.k-1);

    % splitting input into two branches
    branch_1_sampl = samples_rx(1:end-fsync_param.n_subc);
    branch_2_sampl = samples_rx(fsync_param.n_subc/2+1:end-fsync_param.n_subc/2);

    % serial-to-parallel conversion
    branch_1_sampl = reshape(branch_1_sampl,fsync_param.n_subc,numel(branch_1_sampl)/fsync_param.n_subc);
    branch_2_sampl = reshape(branch_2_sampl,fsync_param.n_subc,numel(branch_2_sampl)/fsync_param.n_subc);

    % filter bank receiver
    %branch_1_ppn = fbmc_lib.ppn(mod_param.n_subc, fliplr(filter_func),mod_param.k);
    %branch_2_ppn = fbmc_lib.ppn(mod_param.n_subc, fliplr(filter_func),mod_param.k);
    branch_1_ppn = fbmc_lib.ppn(fsync_param.n_subc, fliplr(fsync_param.filter_func),fsync_param.k);
    branch_2_ppn = fbmc_lib.ppn(fsync_param.n_subc, fliplr(fsync_param.filter_func),fsync_param.k);    

    % expanding the input of the convolution to avoid erros due to cyclic convolution
    branch_1_ppn     = flipud([branch_1_ppn zeros(fsync_param.n_subc,size(branch_1_sampl,2)-1)]);
    branch_2_ppn     = flipud([branch_2_ppn zeros(fsync_param.n_subc,size(branch_2_sampl,2)-1)]);
    branch_1_sampl   = [branch_1_sampl zeros(fsync_param.n_subc,fsync_param.k-1)];
    branch_2_sampl   = [branch_2_sampl zeros(fsync_param.n_subc,fsync_param.k-1)];

    % application of the PPN filter to the FFT input
    branch_1_sampl = ifft(fft(branch_1_ppn,[],2).*fft(branch_1_sampl,[],2),[],2);
    branch_2_sampl = ifft(fft(branch_2_ppn,[],2).*fft(branch_2_sampl,[],2),[],2);

    % fft
    branch_1_symb_fft = fft(branch_1_sampl);
    branch_2_symb_fft = fft(branch_2_sampl);

    % removing the transient received symbols at the beginning (k-1) and the end
    oqam_frame_metric            = zeros(fsync_param.n_subc,2*n_subc_multiples);
    oqam_frame_metric(:,1:2:end) = branch_1_symb_fft(:,fsync_param.k+(0:n_subc_multiples-1));
    oqam_frame_metric(:,2:2:end) = branch_2_symb_fft(:,fsync_param.k+(0:n_subc_multiples-1));
    
    % oqam frame is shifted, so deshift
    oqam_frame_metric = fftshift(oqam_frame_metric,1);

    % extract subcarriers with pilot tones only
    afb_pilot = oqam_frame_metric(fsync_param.pilot_idx,:);
end