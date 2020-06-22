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

function [edges, bin_centers, bin_distance] = hist_edges(n_bins, edge_l, edge_r)

    edge_lr_distance = edge_r - edge_l;
    
    edges = 0:1:n_bins;
    
    scale_fac = edge_lr_distance/n_bins;
    
    edges = edges*scale_fac;
    
    edges = edges + edge_l;
    
    bin_centers = mean([edges(1:end-1);edges(2:end)]);
    
    bin_distance = bin_centers(2) - bin_centers(1);
end

