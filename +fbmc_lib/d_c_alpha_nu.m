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

function d = d_c_alpha_nu(c,alpha,nu,big_K)

%% b values
b=[1 3/4 105/64 675/256 76233/16384 457107/65536 12097169/1048576 70545315/4194304;
    -1 -15/8  -219/64 -6055/1024 -161925/16384  -2067909/131072 -26060847/1048576 0;
    3/4 19/16 1545/512 9765/2048 596277/65536 3679941/262144 394159701/16777216  0;
    -5/8 -123/128 -2289/1024  -34871/8192  -969375/131072 -51182445/4194304 0 0 ;
    35/64 213/256 7797/4096 56163/16384 13861065/2097152 87185895/8388608 0 0;
    -63/128 -763/1024  -13875/8192 -790815/262144 -23600537/4194304     0 0 0 ;
    231/512 1395/2048 202281/131072 1434705/524288  85037895/16777216 0 0 0;
    -429/1024 -20691/32768  -374325/262144 -5297445/2097152  0 0 0 0;
    6435/16384 38753/65536 1400487/1048576  9895893/4194304 0 0 0 0;
    -12155/32768 -146289/262144 -2641197/2097152  0 0 0 0 0;
    46189/131072  277797/524288 20050485/16777216 0 0 0 0 0;
    -88179/262144 -2120495/4194304 0 0 0 0 0 0 ;
    676039/2097152 4063017/8388608 0 0 0 0 0 0 ;
    -1300075/4194304 0 0 0 0 0 0 0;
    5014575/16777216 0 0 0 0 0 0 0];

%% init d
d = 0;
for j = 0:floor((big_K-c)/2)
    d = d + b(c+1,j+1)*exp(-[pi*alpha/2/nu.^2]*[2*j+c]);
end

end % func

