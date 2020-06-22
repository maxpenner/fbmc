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

clear;
close all;
%profile on

% show time
datetime('now')
t_start = tic;

%rng(123);

% define system parameters
samp_rate = 20000000;
fc = 800e6;

% tx
target = fbmc_chain();

% overwrite with profile
samp_rate = profiler_lib.lte(target,'bandwidth_5_0');

% generate frame
samples = target.generate_frame();

% create channel
ch = channel_lib.rf_channel();
ch.index                    = 1;
ch.type                     = 'statistical';
ch.amp                      = 1.0;
ch.noise                    = true;
ch.snr_db                   = 20;
ch.sto                      = 0;
ch.err_phase                = deg2rad(00);
ch.cfo                      = 0.0/100*(1/target.n_subc);
ch.s_samp_rate              = samp_rate;
ch.s_n_subc                 = target.n_subc;
ch.s_max_doppler        	= channel_lib.convert_speed_kmh_2_maxdoppler_hz(500,fc);
ch.s_rayleigh_type          = 'Custom Exponential Decay';
ch.s_rayleigh_tau_mean      = target.n_subc*0.0001*1/samp_rate;
ch.s_rayleigh_tau_rms       = target.n_subc*0.0001*1/samp_rate;
ch.s_rayleigh_gains_yn      = true;
ch.s_random_source          = 'global';
ch.s_rayleigh_random_seed   = 11;
ch.s_noise_random_stream	= RandStream('mt19937ar','Seed', 22);
ch.init_statistical_channel();

% pass samples through channel
samples_ch = ch.pass_samples(samples, 0);

% remove any samples at the end that the channel added
samples_ch = samples_ch(1:numel(samples));

% try decoding the frame
target.decode_frame(samples_ch, ch);

% check for bit errors
[transmitted_bit_errors, BER] = target.check_transmitted_bit_errors();
[user_bit_errors, BER_user] = target.check_user_bit_errors();
[block_errors, BLER_user, throughput, throughput_bpshz] = target.check_block_errors_throughput();

% determine interference for data and pilots
[~, ~, sinr_power_data, ~, ~] = interference_lib.points_real_imag_sinr_hist_db(target,'data',200, -2, 2);
[~, ~, sinr_power_pilot, ~, ~] = interference_lib.points_real_imag_sinr_hist_db(target,'pilot',200, -2, 2);

% show results
fprintf('Signal power before ch:    %f\n', sum(abs(samples).^2)/numel(samples));
fprintf('Signal power after ch:     %f\n', sum(abs(samples_ch).^2)/numel(samples_ch));
fprintf('Points used in d_k_m:      %f\n', target.d_k_m_statistics.usage);
fprintf('Point data SINR:           %f\n', sinr_power_data);
fprintf('Point pilot SINR:          %f\n', sinr_power_pilot);
fprintf('Transmitted bit:           %f\n', numel(target.transmitted_bits));
fprintf('Transmitted bit errors:    %f\n', numel(transmitted_bit_errors));
fprintf('Transmitted BER:           %f %%\n', BER*100.0);
fprintf('User bit:                  %f\n', numel(target.user_bits));
fprintf('User bit errors:           %f\n', numel(user_bit_errors));
fprintf('User BER:                  %f %%\n', BER_user*100.0);
fprintf('User blocks:               %f\n', target.n_blocks);
fprintf('User block errors:         %f\n', numel(block_errors));
fprintf('User BLER:                 %f %%\n', BLER_user*100.0);
fprintf('User throughput:           %f Mbps\n', throughput*samp_rate/1e6);
fprintf('User throughput norm:      %f b/s/hz\n', throughput_bpshz);

% plot time domain signal before and after channel
figure()
    subplot(3,1,1)
    plot(abs(samples));
    xlabel('Time');
    ylabel('Absolute');
    title('Transmitted Samples before Channel');
        subplot(3,1,2)
        plot(abs(samples_ch));
        xlabel('Time');
        ylabel('Absolute');
        title('Transmitted Samples after Channel')
            subplot(3,1,3)
            plot(abs(samples));
            hold on
            plot(abs(samples - samples_ch),'r');
            xlabel('Time');
            ylabel('Absolute');
            title('Transmitted Samples: Difference before minus after Channel')
            legend('Transmitted Samples before Channel', 'Difference due to channel');

% plot scatter
scatterplot(target.frame_rx(target.frame_data_indices.lin));

% restore the frame
if numel(user_bit_errors) == 0
    
    disp('Testing restore function!');
    
    restored_samples = target.generate_frame_restored(target.user_bits, []);

    % plot difference between transmitted samples after channel and restored samples
    figure()
        subplot(3,1,1)
        plot(abs(samples_ch));
        xlabel('Time');
        ylabel('Absolute');
        title('Transmitted Samples after Channel');
            subplot(3,1,2)
            plot(abs(restored_samples));
            xlabel('Time');
            ylabel('Absolute');
            title('Restored samples')
                subplot(3,1,3)
                plot(abs(samples));
                hold on
                plot(abs(samples_ch - restored_samples),'r');
                xlabel('Time');
                ylabel('Absolute');
                title('Difference Sampels after Channel minus restored Samples')
                legend('Transmitted Samples after Channel', 'Difference due to imperfect restore function');
else
    disp('Restore function not tested since there are user bit errors!');
end

% show elapsed time
t_elapsed = toc(t_start);
fprintf('Total elapsed time: %f sec\n', t_elapsed);

%profile off
%profile viewer
disp('Simulation finished')
