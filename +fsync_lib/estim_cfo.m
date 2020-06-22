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

% formula (22) from Thein Eurasip

function cfo_ = estim_cfo(fsync_param, inp_vecs)

    % the cfo correction is performed with 0 -> multiplication with 1
    if fsync_param.cfo_corr_on == false
        cfo_ = 0;
        return;
    end

    % split for readability
    vec_0 = inp_vecs(:,1);
    vec_2 = inp_vecs(:,2);
    
    % cfo unit is Hertz
    cfo_ = sum(conj(vec_0).*vec_2);
    cfo_ = angle(cfo_);
    cfo_ = 1/(2*pi*fsync_param.n_subc)*cfo_;
end