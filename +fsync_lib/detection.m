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

function detection_metric = detection(fsync_param, afb_pilot)
    
    % split up in two delayed branches
    afb_branch_0 = afb_pilot(:,1:end-2);
    afb_branch_1 = afb_pilot(:,3:end);
    
    % symbol wise multiplication
    C_b_A = abs(sum(afb_branch_0.*conj(afb_branch_1),1));
    
    % power
    Q_b = sum(abs(afb_branch_0).^2 + abs(afb_branch_1).^2,1)/2;
    
    % actual detection metric
    detection_metric = C_b_A./Q_b;

    % debugging
    if fsync_param.dbg_fsync > 0
        figure(fsync_param.dbg_fsync)
        clf()
        
        subplot(3,1,1)
        bar(abs(C_b_A));       
        title('C b A');
        
        subplot(3,1,2)
        bar(abs(Q_b));
        title('Q b');
        
        subplot(3,1,3)
        bar(abs(detection_metric));        
        title('detection metric');
        xlabel('symbol index')
    end 
end

