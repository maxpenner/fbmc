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

function candidate_list = candidates(fsync_param, afb_pilot, detection_metric, add_info)

    % create output container
    candidate_list.idx = [];
    candidate_list.cfo = [];
    candidate_list.sto = [];

    % compare with threshold and get indices
    candidate_list.idx = find(detection_metric > fsync_param.coarse_thr);
    
    % for each candidate calculate additional information
    for i = candidate_list.idx
        if strcmp(add_info, 'cfo_only')                        
            candidate_list.cfo(end+1) = fsync_lib.estim_cfo(fsync_param,afb_pilot(:,[i i+2]));
            candidate_list.sto(end+1) = 0;
        elseif strcmp(add_info, 'cfo_and_sto')                        
            candidate_list.cfo(end+1) = fsync_lib.estim_cfo(fsync_param,afb_pilot(:,[i i+2]));
            candidate_list.sto(end+1) = fsync_lib.estim_sto(fsync_param,afb_pilot(:,[i i+2]));
        end       
    end
end