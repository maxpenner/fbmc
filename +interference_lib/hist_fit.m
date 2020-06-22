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

function [propabilities_fitted, distribution_specific] = hist_fit(bin_centers, propabilities_unfitted, distribution_type)

    % determine bin distance
    bin_distance = bin_centers(2) - bin_centers(1);      
    
    if strcmp(distribution_type,'gauss')
        
        % we need the mean
        mean_of_dist = sum(propabilities_unfitted.*bin_centers*bin_distance);

        % we need the standart distribution
        std_of_dist = sqrt(sum(bin_distance*propabilities_unfitted.*(bin_centers-mean_of_dist).^2));        
               
        % calculate fitted values at bin centers 
        propabilities_fitted = 1/sqrt(2*pi*std_of_dist^2)*exp(-(bin_centers-mean_of_dist).^2/(2*std_of_dist^2));
           
        % write distribution specific data
        distribution_specific.type = distribution_type;
        distribution_specific.mean_of_dist = mean_of_dist;
        distribution_specific.std_of_dist = std_of_dist;
        distribution_specific.var_of_dist = std_of_dist^2;         
        
    elseif strcmp(distribution_type,'rayleigh')
        
        % we need the mean
        mean_of_dist = sum(propabilities_unfitted.*bin_centers*bin_distance);

        % we need the standart distribution
        std_of_dist = mean_of_dist/sqrt(pi/2);        
               
        % calculate fitted values at bin centers 
        propabilities_fitted = bin_centers/std_of_dist^2.*exp(-(bin_centers).^2/(2*std_of_dist^2));
           
        % write distribution specific data
        distribution_specific.type = distribution_type;
        distribution_specific.mean_of_dist = mean_of_dist;
        distribution_specific.std_of_dist = std_of_dist;
        distribution_specific.var_of_dist = std_of_dist^2;         
        
    else
        error('Unknown distribution to fit.');
    end
end

