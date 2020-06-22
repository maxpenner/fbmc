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

classdef pulseDesign
    %PULSEDESIGN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access=private, Constant=true)
        egfBKj = [            1,        3/2^2,       105/2^6,       675/2^8,      76233/2^14,    457107/2^16,  12097169/2^20, 70545315/2^22; ...
                             -1,      -15/2^3,      -219/2^6,     -6055/2^10,   -161925/2^14,  -2067909/2^17, -26060847/2^20, 0; ...
                          3/2^2,       19/2^4,      1545/2^9,      9765/2^11,    596277/2^16,   3679941/2^18, 394159701/2^24, 0; ...
                         -5/2^3,     -123/2^7,     -2289/2^10,   -34871/2^13,   -969375/2^17, -51182445/2^22,              0, 0; ...
                         35/2^6,      213/2^8,      7797/2^12,    56163/2^14,  13861065/2^21,  87185895/2^23,              0, 0; ...
                        -63/2^7,     -763/2^10,   -13875/2^13,  -790815/2^18, -23600537/2^22,              0,              0, 0; ...
                        231/2^9,     1395/2^11,   202281/2^17,  1434705/2^19,  85037895/2^24,              0,              0, 0;...
                      -429/2^10,   -20691/2^15,  -374325/2^18, -5297445/2^21,              0,              0,              0, 0;...
                      6435/2^14,    38753/2^16,  1400487/2^20,  9895893/2^22,              0,              0,              0, 0; ...
                    -12155/2^15,  -146289/2^18, -2641197/2^21,             0,              0,              0,              0, 0;...
                     46189/2^17,   277797/2^19, 20050485/2^24,             0,              0,              0,              0, 0; ...
                    -88179/2^18, -2120495/2^22,             0,             0,              0,              0,              0, 0; ...
                    676039/2^21,  4063017/2^23,             0,             0,              0,              0,              0, 0; ...
                  -1300075/2^22,             0,             0,             0,              0,              0,              0, 0; ...
                   5014575/2^24,             0,             0,             0,              0,              0,              0, 0];
        %% PHYDYAS Reference: 
        % [1] Overlapped Complex-Modulated Transmultiplexer Filters With Simplified Design and Superior Stopbands, Shahriar Mirabbasi
        % [2] Small Side-Lobe Filter Design for Multitone Data-Communication Applications, Kenneth W. Martin
        phydyasK = [ 1,          0,          0,          0,          0,           0,          0,           0; ...  [2]
                     1,  sqrt(2)/2,          0,          0,          0,           0,          0,           0; ...  [2]
                     1, 0.91143783, 0.41143783,          0,          0,           0,          0,           0; ...  [2]
                     1, 0.97195983,  sqrt(2)/2, 0.23514695,          0,           0,          0,           0; ...  [2]
                     1, 0.99184131, 0.86541624, 0.50105361, 0.12747868,           0,          0,           0; ...  [1]
                     1, 0.99722723, 0.94136732,  sqrt(2)/2,  0.3373834,  0.07441672,          0,           0; ...  [2]
                     1, 0.99938080, 0.97838560, 0.84390076, 0.53649931,  0.20678881, 0.03518546,           0; ...  [1]
                     1, 0.99988389, 0.99315513, 0.92708081,  sqrt(2)/2, 0.374861854, 0.11680273, 0.01523841]; %  [2]
        % [1] Design of prototype filter for near-perfect-reconstruction overlapped complex-modulated transmultiplexers, S. Mirabbasi and K. Martin
        %% Hermite reference:
        % [1] A Time-Frequency Well-localized Pulse for Multiple Carrier Transmission, Ralf Haas
        % [2] Efficient prototype filter design for Filter Bank Multicarrier (FBMC) System based on Ambiguity function analysis of Hermite
        % [3} Generated coefficients depending on filter length
        % Polynomials, Arun Prakash J 
        hermiteH  = [1, -1.9324881e-3, -7.3110588e-6, -3.1542096e-9, 9.6634138e-13] * 1.1850899; % [1]
        hermiteH2 = [1.1413, -3.1222e-3, -0.4098e-6, 2.8530e-9];  % [2]
        hermiteH3 = {[0.314590334765008 -0.000422016597856764 -4.01193254264472e-06 6.13907744971722e-09], ...  [3]
                     [0.312696647265234 -0.000978735551712819 -9.03966698114626e-07 5.07157212690287e-09], ...  [3]
                     [0.312987153526711 -0.000676326843840724 -1.91521951170748e-06 -3.99863347164311e-10], ...
                     [0.313031823966989 -0.000489223022701810 -2.88258668833100e-06 -1.94727817264764e-09 7.31503407792200e-13],...
                     [0.313031837479309 -0.000489188976996805 -2.88256497353907e-06 -1.94848356690920e-09 7.29710370237617e-13],...
                     [0.313031837479390 -0.000489188976768419 -2.88256497358647e-06 -1.94848357372110e-09 7.29710360640898e-13],...
                     [0.313031837479380 -0.000489188976742743 -2.88256497381146e-06 -1.94848357402318e-09 7.29710360791987e-13],...
                     [0.313031837479380 -0.000489188976742743 -2.88256497381146e-06 -1.94848357402318e-09 7.29710360791987e-13]}; %  [3]
                 
        %% Hanning
        hanningK = [0.5, 0.25];
        %% Hamming
        hammingK = [0.54, 0.23];
        %% Blackman
        blackmanK = [0.42, 0.25, 0.04];
    end
    
    properties
    end
    
    methods(Static=true, Access=public)
        function [impulseResponse, oqam, tauVec, nuVec, pffType] = generate(pulseShape, nCarr, nSymb, tau0, nu0, filterSym, delay)
            if nargin < 4
                tau0 = 0.5;
                nu0 = 1;
            end
            if nargin < 5
                nu0 = 0.5 / tau0;
            end
            if nargin < 6
                filterSym = 'odd';
            end
            if nargin < 7 || isempty(delay)
                delay = 0;
            end
            %% Prepare pulse shape generation
            pulseShape = lower(pulseShape);
            if strcmp(pulseShape(1:3), 'egf')
                alpha = pulseShape(5:end);
                pulseShape = 'egf';
            elseif strcmp(pulseShape(1:4), 'iota')
                alpha = pulseShape(6:end);
                pulseShape = 'iota';
            else
                alpha = 1;
            end
            momentRatios = [];
            if ischar(alpha)
                if strcmp(alpha, 'opt')
                    alpha = -inf;
                    alphaVec = alphaMin:0.05:alphaMax;
                    momentRatios = calcMomentRatio(nCarr, nSymb, alphaVec, tau0, nu0);
                else
                    alpha = str2num(alpha);
                end
            end
            if isempty(alpha)
                alpha = 1;
            end
            
            % setup lattice grid
            if strcmp(pulseShape, 'rect')
                tau_ = 1;
                nu_ = 1;
            else
                tau_ = tau0;
                nu_ = nu0;
            end
            nCarrMax = 32;%nCarr/2-1;
            tauVec = -floor((nSymb+1)/tau_-1)*tau_:tau_:floor((nSymb+1)/tau_-1)*tau_;
            nuVec = -nCarrMax*nu_:nu_:nCarrMax*nu_;
            
            oqam = true;
            if mod(nCarr, 2) == 0 && strcmp(filterSym, 'even')
                t = -nSymb/2+0.5/nCarr:1/nCarr:(nSymb/2-0.5/nCarr);
            else
                t = -nSymb/2:1/nCarr:(nSymb/2-1/nCarr);
            end
            
            t = t - rem(delay, 1/nCarr);
            switch pulseShape
                case 'hamming'
                    impulseResponse = pulseDesign.hamming(nCarr, nSymb, filterSym, delay, t);
                    dispString = 'Generating Hamming pulse shape';
                    pffType = 'timeSpreaded';
                case 'hanning'
                    impulseResponse = pulseDesign.hanning(nCarr, nSymb, filterSym, delay, t);
                    dispString = 'Generating Hanning pulse shape';
                    pffType = 'timeSpreaded';
                case 'blackman'
                    impulseResponse = pulseDesign.blackman(nCarr, nSymb, filterSym, delay, t);
                    dispString = 'Generating Blackman pulse shape';
                    pffType = 'timeSpreaded';
                case 'phydyas'
                    impulseResponse = fbmc_lib.pulseDesign.phydyas(nCarr, nSymb, filterSym, delay, t);
                    dispString = 'Generating Phydyas pulse shape';
                    pffType = 'timeSpreaded';
                case 'hermite'
                    impulseResponse = pulseDesign.hermite(nCarr, nSymb, filterSym, delay, alpha, tau0, t);
                    dispString = ['Generating Hermite pulse shape with alpha = ', num2str(alpha)];
                    pffType = 'isotropic';
                case 'hermite opt'
                    impulseResponse = pulseDesign.hermiteOpt(nCarr, nSymb, filterSym, delay, alpha, tau0, t);
                    dispString = ['Generating optimized Hermite pulse shape with alpha = ', num2str(alpha)];
                    pffType = 'isotropic';
                case 'ofdp'
                    dispString = 'Generating OFDP pulse shape';
                    error('ToDo:Implement', 'Implementation open')
                case 'rect'
                    impulseResponse = pulseDesign.rect(nCarr, nSymb, 0, delay);
                    %             impulseResponse = {pulseDesign.rect(nCarr, nSymb, 0), impulseResponse};
                    dispString = 'Generating Rectangular pulse shape';
                    oqam = false;
                    pffType = 'freqSpreaded';
                case {'egf', 'iota'}
                    dispString = 'Generating EGF/IOTA pulse shape with alpha = ';
                    if alpha == -inf
                        targetRatio = 1/3*delaySpreads(delayIdx) / dopplerShifts(dopplerIdx);
                        %                                 targetRatio = 1/sqrt(2)*delaySpreads(delayIdx) / dopplerSpreads(dopplerIdx);
                        idxLow = find(momentRatios < targetRatio);
                        idxHigh = find(momentRatios >= targetRatio);
                        if isempty(idxHigh)
                            alphaInt = alphaMin;
                        elseif isempty(idxLow)
                            alphaInt = alphaMax;
                        else
                            idxLow = idxLow(1);
                            idxHigh = idxHigh(end);
                            alphaInt = interp1([momentRatios(idxLow), momentRatios(idxHigh)], ...
                                [alphaVec(idxLow), alphaVec(idxHigh)], ...
                                targetRatio);
                        end
                        impulseResponse = pulseDesign.egf(nCarr, nSymb, alphaInt, filterSym, delay, tau0, nu0, t);
                        dispString = [dispString, num2str(alphaInt)];
                    else
                        impulseResponse = fbmc_lib.pulseDesign.egf(nCarr, nSymb, alpha, filterSym, delay, tau0, nu0, t);
                        dispString = [dispString, num2str(alpha)];
                        if alpha < 2
                            pffType = 'timeSpreaded';
                        elseif alpha > 2
                            pffType = 'freqSpreaded';
                        else
                            pffType = 'isotropic';
                        end
                    end
                case 'gaussian'
                    impulseResponse = pulseDesign.gaussian(nCarr, nSymb, alpha, delay, t);
                    dispString = ['Generating Gaussian pulse shape with alpha = ', num2str(alpha)];
                case 'prolate'
                    impulseResponse = pulseDesign.prolate(nCarr, nSymb, delay, t);
                    dispString = 'Generating Prolate pulse shape';
                otherwise
                    error('Input:pulseShape', ['Unknown pulse shape: ', pulseShape])
            end
            %disp([dispString ' with # Carr = ' num2str(nCarr)])
        end
        

        function impulseResponse = egf(nCarr, nSymb, alpha, sym, delay, tau0, nu0, t, K)
            if nargin < 3 || isempty(alpha)
                alpha = 1;
                disp('Applying alpha = 1')
            end
            if nargin < 4 || isempty(sym)
                sym = 'odd';
            end
            if nargin < 5 || isempty(delay)
                delay = 0;
            end
            if nargin < 6 || nargin < 7 || isempty(tau0) || isempty(nu0)
                disp('Applying tau0 = 0.5, nu0 = 1')
                tau0 = 0.5;
                nu0 = 1;
            end
            if nargin < 8 || isempty(t)
                if mod(nCarr, 2) == 0 && strcmp(sym, 'even')
                    t = -nSymb/2+0.5/nCarr:1/nCarr:(nSymb/2-0.5/nCarr);
                else
                    t = -nSymb/2:1/nCarr:(nSymb/2-1/nCarr);
                end
                t = t + rem(delay, 1/nCarr);
            end
            if nargin < 9 || isempty(K)
                K = 14;
            end
            if K > 14
                warning('PulseShape:IOTA:K','Accuracy of calculation of IOTA pulse shape limited to K=14')
                K =14;
            end
            alpha_m = 0.528*nu0^2;
            if alpha_m > alpha || alpha > 1/alpha_m
                %%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%
                %warning('PulseShape:IOTA:alpha',['Calculation of EGF pulse shape is not accurate to match IOTA pulse shape! Alpha out of boundaries (' num2str(alpha) ' not in [' num2str(alpha_m) ', ' num2str(1/alpha_m) '])'])
                %%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%
            end
            dK = fbmc_lib.pulseDesign.calcDk(alpha, nu0, tau0, K);
            dGauss = zeros(size(t));
            dCos = zeros(size(t));
            for k=0:K
                dGauss = dGauss + dK(1, k+1) * fbmc_lib.pulseDesign.sumGAlpha(alpha, t, k, nu0);
                dCos = dCos + dK(2, k+1) * cos(2*pi*k*t/tau0);
            end
            impulseResponse = 1/2 * dGauss .* dCos;
            impulseResponse = impulseResponse / norm(impulseResponse);
            if strcmp(sym, 'odd')
                impulseResponse(1) = 0;
            end
        end
        
        function impulseResponse = hamming(nCarr, nSymbIn, sym, delay, t)
            if nargin < 2 || isempty(nSymbIn)
                nSymbIn = 4;
            end
            nSymb = 1;
            if nargin < 3 || isempty(sym)
                sym = 'odd';
            end
            if nargin < 4 || isempty(delay)
                delay = 0;
            end
            if nargin < 5 || isempty(t)
                if mod(nCarr, 2) == 0  && strcmp(sym, 'even')
                    t = -nSymb/2+0.5/nCarr:1/nCarr:(nSymb/2-0.5/nCarr);
                else
                    t = -nSymb/2:1/nCarr:(nSymb/2-1/nCarr);
                end
                t = t + rem(delay, 1/nCarr);
            end
            impulseResponse = pulseDesign.hammingK(1)*ones(size(t));
            for k=1:1
                impulseResponse = impulseResponse + 2*pulseDesign.hammingK(k+1) * cos(2*pi*k*t/(nSymb));
            end
            if strcmp(sym, 'odd')
%                 impulseResponse(1) = 0;
            end
            figure(999)
            plot(impulseResponse)
            impulseResponse = impulseResponse / norm(impulseResponse);
            if nSymbIn > nSymb
                delta = 0.5 * (nSymbIn - nSymb);
                impulseResponse = [zeros(1, nCarr*delta), impulseResponse, zeros(1, nCarr*delta)];
            end
        end
        
        function impulseResponse = hanning(nCarr, nSymbIn, sym, delay, t)
            if nargin < 2 || isempty(nSymbIn)
                nSymbIn = 4;
            end
            nSymb = 1;
            if nargin < 3 || isempty(sym)
                sym = 'odd';
            end
            if nargin < 4 || isempty(delay)
                delay = 0;
            end
            if nargin < 5 || isempty(t)
                if mod(nCarr, 2) == 0  && strcmp(sym, 'even')
                    t = -nSymb/2+0.5/nCarr:1/nCarr:(nSymb/2-0.5/nCarr);
                else
                    t = -nSymb/2:1/nCarr:(nSymb/2-1/nCarr);
                end
                t = t + rem(delay, 1/nCarr);
            end
            impulseResponse = pulseDesign.hanningK(1)*ones(size(t));
            for k=1:1
                impulseResponse = impulseResponse + 2*pulseDesign.hanningK(k+1) * cos(2*pi*k*t/(nSymb));
            end
            if strcmp(sym, 'odd')
%                 impulseResponse(1) = 0;
            end
            figure(999)
            plot(impulseResponse)
            impulseResponse = impulseResponse / norm(impulseResponse);
            if nSymbIn > nSymb
                delta = 0.5 * (nSymbIn - nSymb);
                impulseResponse = [zeros(1, nCarr*delta), impulseResponse, zeros(1, nCarr*delta)];
            end
        end
        
        function impulseResponse = blackman(nCarr, nSymbIn, sym, delay, t)
            if nargin < 2 || isempty(nSymbIn)
                nSymbIn = 4;
            end
%             if nSymbIn < 3
%                 error('PulseShape:BLACKMAN:length','Filter length out of range (< 2)')
%             end
            nSymb = 1;
            if nargin < 3 || isempty(sym)
                sym = 'odd';
            end
            if nargin < 4 || isempty(delay)
                delay = 0;
            end
            if nargin < 5 || isempty(t)
                if mod(nCarr, 2) == 0  && strcmp(sym, 'even')
                    t = -nSymb/2+0.5/nCarr:1/nCarr:(nSymb/2-0.5/nCarr);
                else
                    t = -nSymb/2:1/nCarr:(nSymb/2-1/nCarr);
                end
                t = t + rem(delay, 1/nCarr);
            end
            impulseResponse = pulseDesign.blackmanK(1)*ones(size(t));
            for k=1:2
                impulseResponse = impulseResponse + 2*pulseDesign.blackmanK(k+1) * cos(2*pi*k*t/(nSymb));
                
            end
            if strcmp(sym, 'odd')
%                 impulseResponse(1) = 0;
            end
            figure(999)
            plot(impulseResponse)
            impulseResponse = impulseResponse / norm(impulseResponse);
            if nSymbIn > nSymb
                delta = 0.5 * (nSymbIn - nSymb);
                impulseResponse = [zeros(1, nCarr*delta), impulseResponse, zeros(1, nCarr*delta)];
            end
        end
        
        function impulseResponse = phydyas(nCarr, nSymbIn, sym, delay, t)
            if nargin < 2 || isempty(nSymbIn)
                nSymbIn = 4;
            end
            if nSymbIn > 8
                warning('PulseShape:PHYDYAS:length','Filter length out of range (>8), using nearest value and appending zeros')
                nSymb = min(nSymbIn, 8);
            else
                nSymb = nSymbIn;
            end
            if nargin < 3 || isempty(sym)
                sym = 'odd';
            end
            if nargin < 4 || isempty(delay)
                delay = 0;
            end
            if nargin < 5 || isempty(t)
                if mod(nCarr, 2) == 0  && strcmp(sym, 'even')
                    t = -nSymb/2+0.5/nCarr:1/nCarr:(nSymb/2-0.5/nCarr);
                else
                    t = -nSymb/2:1/nCarr:(nSymb/2-1/nCarr);
                end
                t = t + rem(delay, 1/nCarr);
            end
            impulseResponse = ones(size(t));
            for k=1:nSymb-1
                impulseResponse = impulseResponse + 2*fbmc_lib.pulseDesign.phydyasK(nSymb, k+1) * cos(2*pi*k*t/(nSymb));
            end
            if strcmp(sym, 'odd')
%                 impulseResponse(1) = 0;
            end
            impulseResponse = impulseResponse / norm(impulseResponse);
            if nSymbIn > nSymb
                delta = 0.5 * (nSymbIn - nSymb);
                impulseResponse = [zeros(1, nCarr*delta), impulseResponse, zeros(1, nCarr*delta)];
            end
        end
        
        function impulseResponse = hermite(nCarr, nSymb, sym, delay, alpha, tau0, t, nH)
            if nargin < 3 || isempty(sym)
                sym = 'odd';
            end
            if nargin < 4 || isempty(alpha)
                alpha = 1;
            end
            if nargin < 5 || isempty(delay)
                delay = 0;
            end
            if nargin < 6 || isempty(tau0)
                tau0 = 0.5;
            end
            if nargin < 7 || isempty(t)
                if mod(nCarr, 2) == 0  && strcmp(sym, 'even')
                    t = -nSymb/2+0.5/nCarr:1/nCarr:(nSymb/2-0.5/nCarr);
                else
                    t = -nSymb/2:1/nCarr:(nSymb/2-1/nCarr);
                end
                t = t + rem(delay, 1/nCarr);
            end
            if nargin < 8 || isempty(nH)
                nH = 4;
            end
            if nH > 4
                warning('PulseShape:Hermite:accuracy','Accuracy of calculation of Hermite pulse shape limited to nH=4')
                nH = 4;
            end
            if alpha ~= 1
                alpha = 1;
                warning('Hermite supports only alpha = 1')
            end
            tau_ = 1/(tau0 * sqrt(2))^2;
            
            hermPoly = pulseDesign.calcHermPoly(t, nH, tau_);
            impulseResponse = zeros(size(t));
            for k=0:nH
                impulseResponse = impulseResponse + pulseDesign.hermiteH(k+1)*hermPoly(4*k+1,:);
            end
            impulseResponse = impulseResponse .* 1/sqrt(2*pi*alpha) .* exp(-pi*t.^2*alpha*tau_);
            impulseResponse = impulseResponse / norm(impulseResponse);
            if strcmp(sym, 'odd')
%                 impulseResponse(1) = 0;
            end
        end
        
        function coeffs = calcHermite2Coeffs(nSymb, tau0, nCarr, nH)
            warning('Debugging needed')
            if nargin < 3 || isempty(nCarr)
                nCarr = 1024;
            end
            if nargin < 4 || isempty(nH)
                nH = 4;
            end
                            
            %% calc Hermite polynomials
            tau_ = 1/(tau0 * sqrt(2))^2;
            alpha = 1;
            t = -nSymb/2+0.5/nCarr:1/nCarr:nSymb/2-0.5/nCarr;
            hermPoly = pulseDesign.calcHermPoly(t, nH, tau_);
            hermPolyIso = zeros(nH + 1, numel(t));
            hermWeight = zeros(nH+1,1);
            for polyIdx = 1:nH+1
%                 hermPolyIso(polyIdx, :) = 1/(2*pi)^(4*(polyIdx-1)/2) * exp(-pi*t.^2*alpha*tau_) .* hermPoly(4*(polyIdx-1)+1,:);
                hermPolyIso(polyIdx, :) = 1/sqrt(2*pi) * exp(-pi*t.^2*alpha*tau_) .* hermPoly(4*(polyIdx-1)+1,:);
                hermWeight(polyIdx) = norm(hermPolyIso(polyIdx, :));
                figure(1); plot(t, hermPolyIso(polyIdx, :)); xlim([-4, 4])
                hermPolyIso(polyIdx, :) = hermPolyIso(polyIdx, :) / hermWeight(polyIdx);
            end
%             hermPolyIso = hermPolyIso / hermWeight(1);
%             figure(1); ambiguity.plotCalc(nCarr, hermPolyIso(1,:))
%             figure(2); ambiguity.plotCalc(nCarr, hermPolyIso(2,:))
%             figure(3); ambiguity.plotCalc(nCarr, hermPolyIso(3,:))
%             figure(4); ambiguity.plotCalc(nCarr, hermPolyIso(4,:))
            %% calc cross ambiguity functions for all hermite polynoms
%             tau = 1 + 0:nSymb-1;
%             nu = 1+ 0:nSymb-1;
            tau = (0:nH-2)*tau0;
            nu = (0:nH-2)*(0.5/tau0);
            % generate RI pattern matrix
            offset = mod((numel(tau) + numel(nu))*0.5, 2);
            riMat = zeros(numel(tau), numel(nu));
            for n = 1:numel(nu)
                for k = 1:numel(tau)
                    riMat(k, n) = 1i^mod(n+k+offset,2);
                end
            end
            crossAmb = zeros(numel(tau), numel(nu), nH+1, nH+1);
            crossAmbClosed = zeros(numel(tau), numel(nu), nH+1, nH+1);
            %% calc cross ambiguity functions in generalized form
            for k = 1:nH+1
                for l = 1:nH+1
                    crossAmb(:,:,k, l) = ambiguity.calc(nCarr, hermPolyIso(k,:), hermPolyIso(l,:), tau, nu, 'unnorm');
                    crossAmb(:,:,k,l) = real(crossAmb(:,:,k,l) ./ riMat);
                end
            end
            %% calc cross ambiguity functions in closed form using generalized Laguerre polynomials
            for k = 4*(0:nH)
                for l = 4*(0:nH)
                    laguerrePoly = LaguerreGen(k, k-l);
                    for tauIdx = 1:numel(tau)
                        for nuIdx = 1:numel(nu)
                            tauVal = tau(tauIdx)*sqrt(2);
                            nuVal = nu(nuIdx)/sqrt(2);
                            pos = tauVal^2 + nuVal^2;
                            L = polyval(laguerrePoly, pi*pos);
                            if abs(tauVal + 1i*nuVal) == 0 && (k-l) < 0
                                res = 0;
                            else
                                res = factorial(l)/2 * exp(-pi/2*pos) * sqrt(2)^(k+l) * sqrt(pi)^(k-l) * ...
                                                                   (tauVal + 1i*nuVal)^(k-l) * exp(1i*pi*tauVal*nuVal) * L;
                            end
                            crossAmbClosed(tauIdx, nuIdx, k/4+1, l/4+1) = res;
                        end
                    end
                    crossAmbClosed(:,:,k/4+1,l/4+1) = real(crossAmbClosed(:,:,k/4+1,l/4+1) ./ riMat);
                end
            end
            %% resort cross ambiguity functions
            crossAmbSort = zeros(nH+1, nH+1, numel(tau), numel(nu));
            for k=1:nH+1
                for l=1:nH+1
                    for tauIdx=1:numel(tau)
                        for nuIdx=1:numel(nu)
                            crossAmbSort(k,l,tauIdx, nuIdx) = crossAmb(tauIdx, nuIdx, k, l);
%                             crossAmbSort(k,l,tauIdx, nuIdx) = crossAmbClosed(tauIdx, nuIdx, k, l);
                        end
                    end
                end
            end
            
            %% initialize iterative process
            d = [1; zeros(numel(tau)*numel(nu)-1, 1)];
            wi = [1, zeros(1, nH)]; % pulseDesign.hermiteH2 .* hermWeight' / hermWeight(1);
            iterationDiff = inf;
            nIter = 1000;
            iteration = 0;
            %% Start iteration
            errorVec = zeros(1, nIter);
            while (iterationDiff > 1e-10 && nIter > iteration)
                iteration = iteration + 1;
                %% calc B = w.'*A
                B = zeros(numel(tau)*numel(nu), nH+1);
%                 B(1, :) = wi * crossAmbSort(:, :, 1, 1);
                rowIdx = 1;
                for m = 1:nH-2
                    for n = m:nH-1
%                 for m = 1:numel(tau)
%                     for n = 1:numel(nu)
%                         [m, n]
                        tmpVec = wi * crossAmbSort(:, :, m, n);
%                         sum(abs(tmpVec))
                        if sum(abs(tmpVec)) <  1e-8
                            continue
                        end
                        B(rowIdx, :) = tmpVec;
                        rowIdx = rowIdx + 1;
                    end
                end
                B(rowIdx:end, :) = [];
                d_ = d(1:size(B,1));
                %% calc w_hat
%                 w_hat = inv(B.' * B) * B.' * d_;
                w_hat = pinv(B) * d_;
                %% update w
%                 wi_new = (wi+3*w_hat.')/4;
                wi_new = (wi+w_hat.')/2;
                iterationDiff = sum(abs(wi_new - wi).^2);
                errorVec(iteration) = norm(B*wi_new.' -d_);
                wi = wi_new;
            end
            figure(98)
            semilogy(errorVec)
            coeffs = wi ./ hermWeight';
            
            %% verification
            pulse = zeros(size(hermPolyIso(1,:)));
            for polyIdx = 1:nH+1
                pulse = pulse + coeffs(polyIdx) * 1/(2*pi)^(4*(polyIdx-1)/2) * exp(-pi*t.^2*alpha*tau_) .* hermPoly(4*(polyIdx-1)+1,:);
%                 pulse = pulse + coeffs(polyIdx) * 1/sqrt(2*pi) * exp(-pi*t.^2*alpha*tau_) .* hermPoly(4*(polyIdx-1)+1,:);
            end
            norm(pulse)
%             pulse = pulse/norm(pulse);
            figure(99)
            ambiguity.plotCalc(nCarr, pulse, [], 'surface', 1)
            
            [ambiguityPart, tau, nu] = ambiguity.calc(nCarr, pulse, [], (-nSymb:nSymb)*tau0, (-8:8)*(0.5/tau0));
            [intPow, sigPow, ~, intMat] = ambiguity.calcInterference(ambiguityPart, nCarr);
            pow2db(sigPow/intPow)
        end
        
        function impulseResponse = hermite2(nCarr, nSymb, alpha, tau0, nH)
            warning('Tidy up this function')
            if nargin < 5
                nH = 3;
            end
            if nH > 3
                warning('PulseShape:Hermite2:accuracy','Accuracy of calculation of Hermite pulse shape limited to nH=4')
                nH = 3;
            end
            tau0 = 1/(tau0 * sqrt(2))^2;
            if mod(nCarr, 2) == 0
                t = -nSymb/2+0.5/nCarr:1/nCarr:(nSymb/2-0.5/nCarr);
            else
                t = -nSymb/2:1/nCarr:(nSymb/2-1/nCarr);
            end
            hermPoly = pulseDesign.calcHermPoly(t, numel(pulseDesign.hermiteH3)-1, tau0);
            impulseResponse = zeros(size(t));
            for k=0:numel(pulseDesign.hermiteH2)-1
                impulseResponse = impulseResponse + pulseDesign.hermiteH2(k+1)*hermPoly(4*k+1,:);
            end
            impulseResponse = impulseResponse .* 1/sqrt(2*pi*alpha) .* exp(-pi*t.^2*alpha*tau0);
%             impulseResponse = impulseResponse .* 1/sqrt(2*pi*alpha) .* exp(-pi*t.^2*alpha*tau0);
            impulseResponse = impulseResponse / norm(impulseResponse);
        end
        
         function impulseResponse = hermiteOpt(nCarr, nSymb, alpha, tau0)
            tau0 = 1/(tau0 * sqrt(2))^2;
            if mod(nCarr, 2) == 0
                t = -nSymb/2+0.5/nCarr:1/nCarr:(nSymb/2-0.5/nCarr);
            else
                t = -nSymb/2:1/nCarr:(nSymb/2-1/nCarr);
            end
            coeffs = pulseDesign.hermiteH3{min(nSymb, 8)};
            hermPoly = pulseDesign.calcHermPoly(t, numel(coeffs)-1, tau0);
            impulseResponse = zeros(size(t));
            for k=0:numel(coeffs)-1
                impulseResponse = impulseResponse + coeffs(k+1)*hermPoly(4*k+1,:);
            end
            impulseResponse = impulseResponse .* 1/sqrt(2*pi*alpha) .* exp(-pi*t.^2*alpha*tau0);
            impulseResponse = impulseResponse / norm(impulseResponse);
        end
        
        function impulseResponse = gaussian(nCarr, nSymb, alpha)
            if nargin < 3 || isempty(alpha)
                alpha = 1;
            end
            
            t = -nSymb/2:1/nCarr:(nSymb/2-1/nCarr);
            impulseResponse = (2*alpha)^(1/4) * exp(-pi*alpha*t.^2);
        end
        
        function impulseResponse = prolate(nCarr, nSymb)
            warning('Tidy up this function')
            
            tau0 = 0.5;
            
            seq_length = ceil(nCarr*nSymb);
            time_halfbandwidth = nSymb;
            num_seq = 2*(time_halfbandwidth)-1;
            %Obtain DPSSs
            [dps_seq, ~] = dpss(seq_length,time_halfbandwidth,num_seq);
            
            dps_seq_freq = util.calcTransferFunc(dps_seq(:, 1).', nCarr);
            
%             t = -nSymb/2:1/nCarr:(nSymb/2-1/nCarr);
            t = (0:numel(dps_seq_freq)-1) / nCarr;
            dps_seq_freq_nyq = sqrt(dps_seq_freq.^2 .* cos(2*pi*t/8));
%             figure(1)
%             plot(t, 10*log10(abs(dps_seq_freq_nyq)))
%             figure(2)
%             plot(abs(ifft(dps_seq_freq_nyq)))
            
            extImpRes = circshift(ifft(dps_seq_freq_nyq), [0, nCarr*nSymb/2]);
%             figure(3)
%             plot(abs(extImpRes))
%             figure(4)
%             plot(abs(extImpRes(1:nCarr*nSymb)))
            impulseResponse = extImpRes(1:nCarr*nSymb);
            impulseResponse = impulseResponse / norm(impulseResponse);
%             k = 1/(8*pi);
%             impulseResponse = dps_seq(:, 1).' .* sqrt(cos(2*pi*k*t/tau0));
        end
        
        function impulseResponse = ofdp(nCarr, nSymb)
            warning('Tidy up this function')
            
            tau0 = 0.5;
            
            seq_length = ceil(nCarr*nSymb*2);
            time_halfbandwidth = 2;
            num_seq = 7;%2*(time_halfbandwidth)-1;
            %Obtain DPSSs
            [dps_seq, lamda] = dpss(seq_length,time_halfbandwidth,num_seq);
            
            coeffs = [0.8589, 0, -0.4142, 0, 0.0355, 0, 0.0482].';
            nyquist_pulse = dps_seq * coeffs;
            extImpRes = circshift(ifft(nyquist_pulse), [0, nCarr*nSymb/2]);
            impulseResponse = ifft(sqrt(fft(nyquist_pulse)));
            impulseResponse = impulseResponse / norm(impulseResponse);            
        end
        
        function impulseResponse = rect(nCarr, nSymb, cp, delay)
            if nargin < 3 || isempty(cp)
                cp = 0;
            end
            cpLeft = 0;%ceil(cp*nCarr/2)/nCarr;
            cpRight = 2*ceil(cp*nCarr/2)/nCarr;
            if mod(nCarr, 2) == 0
                t = -cpLeft-1/2+0.5/nCarr:1/nCarr:(1/2-0.5/nCarr)+cpRight;
            else
                t = -cpLeft-1/2:1/nCarr:(1/2-1/nCarr)+cpRight;
            end
            t = t + rem(delay, 1/nCarr);
            impulseResponse = 1/sqrt(nCarr) * ones(size(t));
            if nSymb > 1
                delta = 0.5 * (nSymb - 1);
                impulseResponse = [zeros(1, nCarr*(delta-cpLeft)), impulseResponse, zeros(1, nCarr*(delta-cpRight))];
            end
        end
        
        function hermPoly = calcHermPoly(t, nH, tau0)
            % Calculate Hermite polynoms up to desired order and with required normalization for time axis (tau0)
            hermPoly = zeros(4*nH+1, numel(t));
            gamma = sqrt(2*pi*tau0);
            tHerm = t*gamma;
            % init recursive loop
            hermPoly(1, :) = 1;
            hermPoly(2, :) = -2*tHerm;
            for loopIdx = 2:4*nH
                hermPoly(loopIdx+1,:) = -2*(tHerm.*hermPoly(loopIdx, :) + (loopIdx-1).*hermPoly(loopIdx-1, :));
            end

        end
    end
    
    methods(Static=true, Access=private)
        function d_k = calcDk(alpha, nu_0, tau_0, K)
            % Calculate the d_k's for EGF generation
            d_k = zeros(2, K+1);
            for k=0:K
                j_k = 0:floor((K-k)/2);
                d_k(1, k+1) = sum(fbmc_lib.pulseDesign.egfBKj(k+1,j_k+1) .* exp(-pi*alpha/(2*nu_0^2)*(2*j_k+k)));
                d_k(2, k+1) = sum(fbmc_lib.pulseDesign.egfBKj(k+1,j_k+1) .* exp(-pi/alpha/(2*tau_0^2)*(2*j_k+k)));
            end
        end
        
        function g = sumGAlpha(alpha, t, k, nu_0)
            g = (2*alpha)^(1/4)*(exp(-pi*alpha*(t+k/nu_0).^2) + exp(-pi*alpha*(t-k/nu_0).^2));
        end
        
    end
    
end

