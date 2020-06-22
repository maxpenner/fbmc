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

function [convEncoder,vitDecoder] = coder_init(obj)

    % convolutional coder with contraint length 7
    if strcmp(obj.encoding,'convolutional_7') == true
        
        convEncoder = comm.ConvolutionalEncoder(   poly2trellis(7, [171 133]), ...
                                                        'TerminationMethod', 'Truncated');
                                                        %'FinalStateOutputPort', true);
        vitDecoder = comm.ViterbiDecoder(poly2trellis(7, [171 133]), ...
                                                        'InputFormat', 'Hard', ...
                                                        'TerminationMethod', 'Truncated', ...
                                                        'TracebackDepth', max(34,ceil(2.5*obj.denumerator*(7-1))));
        if obj.numerator == 1 && obj.denumerator == 2
            convEncoder.PuncturePatternSource = 'None';
            vitDecoder.PuncturePatternSource = 'None';
        elseif obj.numerator == 2 && obj.denumerator == 3
            convEncoder.PuncturePatternSource = 'Property';
            convEncoder.PuncturePattern = [1;1;1;0];
            vitDecoder.PuncturePatternSource = 'Property';
            vitDecoder.PuncturePattern = [1;1;1;0];
        elseif obj.numerator == 3 && obj.denumerator == 4
            convEncoder.PuncturePatternSource = 'Property';
            convEncoder.PuncturePattern = [1;1;1;0;0;1];
            vitDecoder.PuncturePatternSource = 'Property';
            vitDecoder.PuncturePattern = [1;1;1;0;0;1];
        else
            error('Unknown code rate.');
        end
        
    % convolutional coder with contraint length 9        
    elseif strcmp(obj.encoding,'convolutional_9') == true
        
        % source: https://www.3gpp2.org/Public_html/Specs/C.S0002-D_v2.0_051006.pdf
        convEncoder = comm.ConvolutionalEncoder(   poly2trellis(9, [765 671 513 473]), ...
                                                        'TerminationMethod', 'Truncated');
                                                        %'FinalStateOutputPort', true);
        vitDecoder = comm.ViterbiDecoder(poly2trellis(9, [765 671 513 473]), ...
                                                        'InputFormat', 'Hard', ...
                                                        'TerminationMethod', 'Truncated', ...
                                                        'TracebackDepth', max(34,ceil(2.5*obj.denumerator*(9-1))));
                                                    
        if obj.numerator ~= 1 || obj.denumerator ~= 4
            error('Unknown code rate.');
        end                                                    
    else
        error('Unknown coder type.');
    end
end

