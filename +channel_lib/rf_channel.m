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

%% channel parameters from ETSI:
%   page 29:    https://www.etsi.org/deliver/etsi_ts/136100_136199/136116/11.04.00_60/ts_136116v110400p.pdf
%   page 861:   https://www.etsi.org/deliver/etsi_ts/136100_136199/136101/13.03.00_60/ts_136101v130300p.pdf
classdef rf_channel < handle
    
    properties
        
        index;                  % index of the channel
        
        type;                   % deterministic or statistical
        
        % for all channel types
        amp;                    % amplitude as linear factor
        noise;                  % add noise to signal, true or false
        snr_db;                 % signal-to-noise-ratio in dB with assumed signal power of 1 over the full bandwidth of 1 Hz
        
        % for deterministic channel only
        sto;                    % symbol timing offset im samples
        err_phase;              % error phase in radians
        cfo;                    % carrier frequency offset in Hertz
        
        % for statistical channel only
        s_samp_rate;            % system sampling rate in Samples/s
        s_n_subc;               % number of subcarriers
        s_max_doppler;      	% max doppler in Hertz
        s_rayleigh_type;        % several types, e.g. 'Custom Exponential Decay'
        s_rayleigh_tau_mean;    % exponential decay, mean of tau in seconds
        s_rayleigh_tau_rms;     % exponential decay, rms of tau in seconds
        s_rayleigh_gains_yn;    % activate of deactivate output of the path gains
        
        % GLOBAL:
        %   Both the matlab rayleigh object and the awgn function use the global random stream.
        %   A call to reset only resets the filter of the matlab rayleigh object.
        % LOCAL:
        %   Both the matlab rayleigh object and the awgn function use an INDIVIDUAL random stream.
        %
        %       The matlab rayleigh object uses an internal mt19937ar generator with the seed "s_rayleigh_random_seed".
        %
        %       The awgn noise function uses an mt19937ar generator "s_noise_random_stream" with a externally set seed.
        %
        %   A call to reset resets the matlab rayleigh object to the default state as when it was generated (filter and s_rayleigh_random_seed).
        %   A call to reset can also reset the random stream of the noise to the default state as when it was initialized with the externally set seed.
        s_random_source;      	% can be global or local
        s_rayleigh_random_seed;	% if local, this is the initial seed for the matlab rayleigh object
        s_noise_random_stream;	% if local, this is the random stream for the awgn function
        
        % activate of deactivate output of the path gains -> will be set when function pass_samples is called
        s_rayleigh_gains;
        
        % appended samples by channel -> will be set in init function
        s_appendix;
        
        % reference to the matlab object -> will be set in init function
        s_ml_rayleigh;
        
        % the statistics before normalization -> will be set in init function
        s_pathDelays_beforeInterpolation;
        s_avgPathGains_beforeInterpolation;
        
        % statistics
        n_frames_passed;
    end
    
    methods (Static = true, Access = public)
       
       %% initialization is done outside of this function
       function obj = rf_channel()
           
            obj.index               = [];
            
            obj.type                = [];
            
            obj.amp                 = [];
            obj.noise               = [];
            obj.snr_db              = [];
                        
            obj.sto                 = [];
            obj.err_phase           = [];
            obj.cfo                 = [];
            
            obj.s_samp_rate             = [];
            obj.s_n_subc                = [];
            obj.s_max_doppler        	= [];
            obj.s_rayleigh_type         = [];
            obj.s_rayleigh_tau_mean     = [];
            obj.s_rayleigh_tau_rms      = [];
            obj.s_rayleigh_gains_yn     = [];
            
            obj.s_random_source        	= [];
            obj.s_rayleigh_random_seed 	= [];
            obj.s_noise_random_stream	= [];
            
            obj.s_rayleigh_gains = [];
            
            % will be set to different value if statistical channel is initialized
            obj.s_appendix = 0;
            
            obj.s_ml_rayleigh = [];           
            
            obj.s_pathDelays_beforeInterpolation    = [];
            obj.s_avgPathGains_beforeInterpolation  = [];
            
            obj.n_frames_passed = 0;
       end
    end
    
    methods (Static = false, Access = public)
        
        %% initialize the statistical rf channel
        function init_statistical_channel(self)
            
            % calculate sampling time from sampling rate
            Ts = 1/self.s_samp_rate;
            
            % we are using a rayleigh channel, so first load the correct channel type
            if strcmp(self.s_rayleigh_type, 'Custom Exponential Decay')
                
                % SECURITY
                if self.s_rayleigh_tau_mean <= 0 || self.s_rayleigh_tau_rms <= 0
                    error('Custom exponential decay must have a tau mean and rms larger than 0.');
                end
                
                % we do not interpolate
                self.s_pathDelays_beforeInterpolation = [];
                self.s_avgPathGains_beforeInterpolation = [];
                
            elseif strcmp(self.s_rayleigh_type, 'Extended Pedestrian A model')
                
                self.s_pathDelays_beforeInterpolation = [0, 30, 70, 90, 110, 190, 410]*1e-9;
                self.s_avgPathGains_beforeInterpolation = [0, -1, -2, -3, -8, -17.2, -20.8];
                
            elseif strcmp(self.s_rayleigh_type, 'Extended Vehicular A model')
                
                self.s_pathDelays_beforeInterpolation = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510]*1e-9;
                self.s_avgPathGains_beforeInterpolation = [0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];
                
            elseif strcmp(self.s_rayleigh_type, 'Extended Typical Urban model')
                
                self.s_pathDelays_beforeInterpolation = [0, 50, 120, 200, 230, 500, 1600, 2300, 5000]*1e-9;
                self.s_avgPathGains_beforeInterpolation = [-1, -1, -1, 0, 0, 0, -3.0, -5.0, -7.0];
                
            else
                
                % load a standard model from matlab:
                % source: https://de.mathworks.com/help/comm/ref/stdchan.html
                temporary_ml_rayleigh = stdchan(self.s_rayleigh_type, self.s_samp_rate, self.s_max_doppler);
                
                % extract path delays and path gains
                self.s_pathDelays_beforeInterpolation = temporary_ml_rayleigh.PathDelays;
                self.s_avgPathGains_beforeInterpolation = temporary_ml_rayleigh.AvgPathGains;
                
                % release the temporary rayleigh channel
                delete(temporary_ml_rayleigh);
            end
            
            % if we are not using a custom exponential decay, we must extract tau mean and rms from the loaded model
            if strcmp(self.s_rayleigh_type, 'Custom Exponential Decay') == false
                
                % SECURITY
                if numel(self.s_pathDelays_beforeInterpolation) == 0 || numel(self.s_avgPathGains_beforeInterpolation) == 0
                    error('Empty path delays and average path gains. Cannot interpolate.');
                end                
                
                % convert path gains to linear space
                s_avgPathGains_beforeInterpolation_linear = db2pow(self.s_avgPathGains_beforeInterpolation);

                % calculate the weighted mean of the path delays
                self.s_rayleigh_tau_mean = sum(s_avgPathGains_beforeInterpolation_linear.*self.s_pathDelays_beforeInterpolation)/sum(s_avgPathGains_beforeInterpolation_linear);

                % calculate the weighted rms of the path delays
                self.s_rayleigh_tau_rms = sqrt(sum(s_avgPathGains_beforeInterpolation_linear.*(self.s_pathDelays_beforeInterpolation.^2))/sum(s_avgPathGains_beforeInterpolation_linear) - self.s_rayleigh_tau_mean^2);                
            end

            % until which point will we calculate weights?
            weight_limit_dB = -20;
            weight_limit_lin = db2pow(weight_limit_dB);
            
            % how many samples do we need until we reach the limit
            n_points = ceil(log(weight_limit_lin)/(-Ts/self.s_rayleigh_tau_rms));
            
            % SECURITY
            if n_points <= 0
                error('Zero points from interpolation.');
            end

            % create array of relevant points
            points = 0:1:(n_points-1);
            
            % calculate path delays
            pathDelays = points*Ts;

            % calculate path gains
            avgPathGains_lin = exp(-points*Ts/self.s_rayleigh_tau_rms);

            % recalculate path gains in dB
            avgPathGains = pow2db(avgPathGains_lin);
            
            % replace minus infinity with -1000 dB
            avgPathGains(~isfinite(avgPathGains)) = -1000;
            
            % create channel object
            if strcmp(self.s_random_source,'global') == true
                
                self.s_ml_rayleigh = comm.RayleighChannel(	'SampleRate', self.s_samp_rate, ...
                                                            'PathDelays', pathDelays, ...
                                                            'AveragePathGains', avgPathGains, ...
                                                            'MaximumDopplerShift', self.s_max_doppler, ...
                                                            'FadingTechnique', 'Sum of sinusoids', ...
                                                            'InitialTimeSource', 'Input Port');
            
            elseif strcmp(self.s_random_source,'local') == true

                self.s_ml_rayleigh = comm.RayleighChannel(	'SampleRate', self.s_samp_rate, ...
                                                            'PathDelays', pathDelays, ...
                                                            'AveragePathGains', avgPathGains, ...
                                                            'MaximumDopplerShift', self.s_max_doppler, ...
                                                            'FadingTechnique', 'Sum of sinusoids', ...
                                                            'InitialTimeSource', 'Input Port', ...
                                                            'RandomStream','mt19937ar with seed', ...
                                                            'Seed', self.s_rayleigh_random_seed);
            else
                error('Unknown randomness source.');
            end
                                                    
            % we can activate and deactivate the output of path gains
            if self.s_rayleigh_gains_yn == true
                self.s_ml_rayleigh.PathGainsOutputPort = true;
            end
                                                    
            % determine how many samples the channel is adding
            self.s_appendix = ceil(max(pathDelays)/Ts);
        end
        
        %% pass samples to the channel at a certain time
        function [samples_ch] = pass_samples(self, samples, sample_time_u64)
            
            % we always need an input time
            if nargin < 3
                error('Incorrect number of inputs. We need a sample time.');
            end
            
            % first adjust the amplitude of the samples, same for each channel model
            samples_ch = samples*self.amp;            
            
            % pass through channel
            if strcmp(self.type, 'awgn')
                
                % if the channel is awgn, do nothing
                
            elseif strcmp(self.type, 'deterministic')
                
                samples_ch = channel_lib.sto_smp_ch(self, samples);
                samples_ch = channel_lib.channel_coef_ch(self, samples_ch);
                samples_ch = channel_lib.cfo_ch(self, samples_ch);
                
            elseif strcmp(self.type, 'statistical')
                
                % the offsets has to be passed as a multiple of Ts = 1/self.samp_rate
                % source: https://de.mathworks.com/help/matlab/matlab_oop/property-attributes.html
                sample_time = double(sample_time_u64);
                time_multiple_period = sample_time*1/self.s_samp_rate;
                
                % the channel adds samples at the end
                samples_ch = [samples_ch, zeros(1,self.s_appendix)];
                
                % complex transpose
                samples_ch = samples_ch.';
                
                % pass through rayleigh channel
                if self.s_rayleigh_gains_yn == true
                    
                    % we save the gains of the channel -> can be very large array
                    [samples_ch, self.s_rayleigh_gains] = self.s_ml_rayleigh(samples_ch, time_multiple_period);
                    
                    % reshape path gains
                    self.s_rayleigh_gains = self.s_rayleigh_gains.';
                else
                    samples_ch = self.s_ml_rayleigh(samples_ch, time_multiple_period);
                end
                
                % reshape to row vector
                samples_ch = samples_ch.';
                
            else
                error('Channel neither awgn, nor deterministic nor statistical.');
            end            
            
            % we might wanna add noise
            if self.noise == true
                
                if strcmp(self.s_random_source,'global') == true
                    
                    samples_ch = awgn(samples_ch, self.snr_db, pow2db(1));
                    
                elseif strcmp(self.s_random_source,'local') == true
                    
                    samples_ch = awgn(samples_ch, self.snr_db, pow2db(1), self.s_noise_random_stream);
                    
                else
                    error('Unknown random stream type for noise.');
                end
            end                        
            
            % stats
            self.n_frames_passed = self.n_frames_passed + 1;
        end
        
        % Reset random streams for rayleigh and for noise.
        % After reset, passing the same samples at the same time should create equal output.
        function [] = reset_random_streams(self, reset_type)
            
            if isempty(reset_type)
                error('Reset type must be defined.');
            end
                        
            if strcmp(reset_type, 'matlab_channel_object_only')
                
                % matlab rayleigh object
                self.s_ml_rayleigh.reset();
                
            elseif strcmp(reset_type, 'matlab_channel_object_and_noise')
                
                if strcmp(self.s_random_source, 'global')
                    error('Trying to reset random stream for noise, which is possible only if source is local. But random stream type is set to global.');
                end
               
                % matlab rayleigh object
                self.s_ml_rayleigh.reset();
               
                % rest random stream for awgn noise
                self.s_noise_random_stream.reset();
                
            else
                error('Unknown reset option.');
            end
        end
        
        %% debug plot between actual channel and interpolation
        function [] = debug_simulate_plot(self)
            
            %% first plot the path gains
            % calculate sampling and symbol time
            Ts = 1/self.s_samp_rate;
            Tsymbol = Ts*self.s_n_subc;
            
            % extract data from interpolation
            pathDelays_interpolated = self.s_ml_rayleigh.PathDelays;
            avgPathGains_interpolated = self.s_ml_rayleigh.AveragePathGains;
            
            % plot loaded model, interpolation and symbol length in dB
            figure()
            subplot(2,1,1);
            plot(self.s_pathDelays_beforeInterpolation, self.s_avgPathGains_beforeInterpolation, '-*b');
            hold on            
            plot(pathDelays_interpolated, avgPathGains_interpolated, '-*r');
            plot([0, Tsymbol], [avgPathGains_interpolated(1), avgPathGains_interpolated(1)], '-*k', 'LineWidth', 2)
            xlabel('Time in Seconds');
            ylabel('Path gain in dB');
            title('Comparison of loaded model and interpolation in dB');
            if strcmp(self.s_rayleigh_type, 'Custom Exponential Decay') == false
                legend('loaded model', 'interpolation', 'symbol length', 'Location', 'southeast');
            else
                legend('interpolation', 'symbol length', 'Location', 'southeast');
            end
            
            % display values for orientation
            text(0,-12, ['Tau mean:', num2str(self.s_rayleigh_tau_mean)]);
            text(0,-15, ['Tau rms:', num2str(self.s_rayleigh_tau_rms)]);

            % plot loaded model, interpolation and symbol length in linear space
            subplot(2,1,2);
            plot(self.s_pathDelays_beforeInterpolation, db2pow(self.s_avgPathGains_beforeInterpolation), '-*b');
            hold on
            plot(pathDelays_interpolated, db2pow(avgPathGains_interpolated), '-*r');
            plot([0, Tsymbol], [db2pow(avgPathGains_interpolated(1)), db2pow(avgPathGains_interpolated(1))], '-*k', 'LineWidth', 2)
            xlabel('Time in Seconds');
            ylabel('Path gain in linear space');
            title('Comparison model and interpolation in linear space');
            if strcmp(self.s_rayleigh_type, 'Custom Exponential Decay') == false
                legend('loaded model', 'interpolation', 'symbol length', 'Location', 'southeast');
            else
                legend('interpolation', 'symbol length', 'Location', 'southeast');
            end
            
            % display values for orientation
            text(0,0.3, ['Tau mean:', num2str(self.s_rayleigh_tau_mean)]);
            text(0,0.1, ['Tau rms:', num2str(self.s_rayleigh_tau_rms)]);            

            %% next simulate the path gains
            % reset the channel before simulation
            reset(self.s_ml_rayleigh);
            
            % simulate the channel and check parameters
            n_realisations = 200;
            n_samp = ceil(10.0*max(pathDelays_interpolated)/Ts);

            % power average
            power_avg = zeros(n_samp,1);

            % time base
            t = 0:1:(n_samp-1);
            t = t*Ts;

            % iterate over multiple iterations
            for i=1:1:n_realisations

                % create samples
                x = ones(n_samp, 1);
                
                % pass through rayleigh channel
                if self.s_rayleigh_gains_yn == true
                    [y, self.s_rayleigh_gains] = self.s_ml_rayleigh(x, 123456789);
                    
                    self.s_rayleigh_gains = self.s_rayleigh_gains.';
                else
                    y = self.s_ml_rayleigh(x, 123456789);
                end
                
                % tranpose x and y
                x = x.';
                y = y.';
                
                % DEBUGGING: check how matlab calculates the channel output
                y_should_be = 0;
                n_delays = numel(self.s_ml_rayleigh.PathDelays);
                n_random = n_delays + 16;
                for j = n_delays : -1 : 1
                    y_should_be = y_should_be + sum(x(n_random + (j-n_delays))*self.s_rayleigh_gains(j,n_random));
                end
                if abs(y(n_random) - y_should_be) > 10e-6
                    error('Wrong calculation for specific point.');
                end

                % reset the channel
                reset(self.s_ml_rayleigh);

                % save iteration
                power_avg = power_avg + abs(y).^2;
            end

            % not the correct orientation
            power_avg = power_avg.';

            % average the power
            power_avg = power_avg/n_realisations;
            
            % make copy in db
            power_avg_db = pow2db(power_avg);            

            % plot delays and power from input to rayleigh
            figure()
            subplot(2,1,1);
            plot(pathDelays_interpolated, avgPathGains_interpolated, '-*');
            hold on
            plot(t, power_avg_db, '-r');
            plot([0, Tsymbol], [db2pow(avgPathGains_interpolated(1)), db2pow(avgPathGains_interpolated(1))], '-*k', 'LineWidth', 2)
            xlabel('Time in Seconds');
            ylabel('Path gain in dB');
            title('Interpolated Average Path Gains vs. summed up power in dB');
            legend('interpolation', 'simulated summed up power', 'symbol length', 'Location', 'southeast');      

            % plot delays and power from input to rayleigh
            subplot(2,1,2);
            plot(pathDelays_interpolated, db2pow(avgPathGains_interpolated), '-*');
            hold on
            plot(t, power_avg, '-r');
            plot([0, Tsymbol], [db2pow(avgPathGains_interpolated(1)), db2pow(avgPathGains_interpolated(1))], '-*k', 'LineWidth', 2)
            xlabel('Time in Seconds');
            ylabel('Path gain in linear space');
            title('Interpolated Average Path Gains vs. summed up power in linear space');
            legend('interpolation', 'simulated summed up power', 'symbol length', 'Location', 'southeast');
        end
    end
end
