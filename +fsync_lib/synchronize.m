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

% Takes IQ-Samples and returns frames found in the stream.
% Frames is array with samples and offsets of cut outs.

function [frames] = synchronize(fsync_param, samples_rx)

    % SECURITY: length of IQ-stream must be a multiple of n_subc
    if mod(numel(samples_rx), fsync_param.n_subc)
        error('Input IQ stream is not a multiple of n_subc.');
    end
    
    % SECURITY: IQ-stream minimum length
    if numel(samples_rx)/fsync_param.n_subc < 10
        error('Input IQ-Stream too short.');
    end   
    
    % 1st stage afb, only pilot tones are returned
    afb_pilot = fsync_lib.afb(fsync_param, samples_rx);
    
    % coarse detection metric, threshold decision
    detection_metric = fsync_lib.detection(fsync_param, afb_pilot);     
    
    % two concepts from Thein Eurasip
    if fsync_param.fsync_concept == 1
        
        % extract candidates for frames with cfo and sto estimation
        candidate_list = fsync_lib.candidates(fsync_param, afb_pilot, detection_metric, 'cfo_and_sto');
        
        % create container for detected frames
        frames.found        = numel(candidate_list.idx);                                % number of frames found        
        frames.samp_orig    = zeros(numel(candidate_list.idx),fsync_param.frame_len);   % samples originally received
        frames.samp         = zeros(numel(candidate_list.idx),fsync_param.frame_len);   % samples after cfo correction
        frames.offs         = zeros(numel(candidate_list.idx),1);                       % offset in samples_rx

        % perform clean cutouts and correct cfo
        for i = 1:numel(candidate_list.idx)
            tau = (candidate_list.idx(i)-1)*fsync_param.n_subc/2 + 1 + candidate_list.sto(i);
            frames.samp_orig(i,:) = samples_rx(tau:tau+fsync_param.frame_len-1);
            frames.samp(i,:) = fsync_lib.correct_cfo(frames.samp_orig(i,:), candidate_list.cfo(i));
            frames.offs(i) = tau;
        end        
        
    elseif fsync_param.fsync_concept == 2
        
        % extract candidates for frames with cfo estimation only
        candidate_list = fsync_lib.candidates(fsync_param, afb_pilot, detection_metric, 'cfo_only');
        
        % create container for detected frames
        frames.found        = numel(candidate_list.idx);	% number of frames found        
        frames.samp_orig    = [];                           % samples originally received
        frames.samp         = [];                           % samples after cfo correction
        frames.offs         = [];                           % offset in samples_rx        
        
        % for each detected frame
        for i = 1:numel(candidate_list.idx)
            
            % extract relevant samples for this candidate
            tau = (candidate_list.idx(i)-1)*fsync_param.n_subc/2 + 1;
            sec = 4;
            sec_samp = sec*fsync_param.n_subc/2;
            
            % if extract lies outside of reasonable limits skip candidate
            start_idx = tau-sec_samp;
            end_idx = tau+fsync_param.frame_len+sec_samp-1;
            if start_idx < 1 || end_idx > numel(samples_rx)
                continue;
            end          
            
            % cut out samples (eff = effective)
            samples_rx_eff_orig = samples_rx(tau-sec_samp:tau+fsync_param.frame_len+sec_samp-1);
            
            % we create a copy for cfo correction
            samples_rx_eff = samples_rx_eff_orig;
            
            % 1st stage cfo correction
            samples_rx_eff = fsync_lib.correct_cfo(samples_rx_eff, candidate_list.cfo(i));

            % 2nd stage afb
            afb_pilot = fsync_lib.afb(fsync_param, samples_rx_eff);

            % 2nd stage cfo + 1st sto estimation
            candidate_list.cfo(i) = fsync_lib.estim_cfo(fsync_param,afb_pilot(:,[sec+1 sec+3]));
            candidate_list.sto(i) = fsync_lib.estim_sto(fsync_param,afb_pilot(:,[sec+1 sec+3]));
            
            % perform clean cutouts, correct cfo and append to frame list
            tau = sec_samp+1+candidate_list.sto(i);
            frames.samp_orig(end+1,:) = samples_rx_eff_orig(tau:tau+fsync_param.frame_len-1);
            frames.samp(end+1,:) = fsync_lib.correct_cfo(samples_rx_eff(tau:tau+fsync_param.frame_len-1), candidate_list.cfo(i));
            frames.offs(end+1) = (candidate_list.idx(i)-1)*fsync_param.n_subc/2 + 1 + candidate_list.sto(i);
        end
    end
    
    % debugging
    if fsync_param.dbg_fsync > 0
        disp(' ');
        disp('fsync_lib->synchronize: cfo and sto:')        
        for i=1:numel(candidate_list.idx)
            fprintf('Frame Candidate: %d \n', i);
            fprintf('Frame Candidate Index: %d \n', candidate_list.idx(i));
            fprintf('Frame Candidate CFO (perc. subc sp): %d \n', candidate_list.cfo(i)/(1/fsync_param.n_subc)*100);
            fprintf('Frame Candidate STO: %d \n', candidate_list.sto(i));
        end
    end
end