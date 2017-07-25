% Matlab sub-function for EE233 Lab 4 Filters
%
% This function is used to find
%     center frequency, f0
%     gain at center frequency, G0
%     cutoff frequencies, fc1, fc2
%     plot gain in frequency domain
% Given by 
%     Resistance (R1,R2,R3,R4,R5)
%     Capacitance (C1,C5)
%     Left portion of potentiometer (gamma)
%     frequency range
%     Command: SingleCalculation, VectorCalculation
%
% IF R1 = R2 and  R3 = R4
%   If gamma > 0.5, the gain is always less than 1, so it is a band-stop filter
%   If gamma < 0.5, the gain is always greater than 1, so it is a band-pass filter
%   If gamma = 0.5, the gain is always equal to 1
%
% Notice: Input references should all be considered as numbers


function [fc, Gc, fc1, fc2] = Lab4_filter_Gain(R1, R2, R3, R4, R5, C1, ...
                                      C2, gamma, minfreq, maxfreq, command)
    % angular frequency range
    samples = 1e6;
    freq = linspace(minfreq, maxfreq, samples);
    w = 2 * pi * freq;
    
    % s domain in steady state
    s = 1i * w;
    
    % delta-Y transformation
    Za = gamma .* R5 ./ (1 + s * R5 * C1);
    Zb = (1 - gamma) * R5 ./ (1 + s * R5 * C1);
    Zc = gamma * (1 - gamma) * R5 ^2 * C1 * s ./ (1 + s * R5 * C1);
    
    % transfer function
    Hs = (1 ./ (R3 + Za) + (Zc + 1 ./ (s .* C2)) ./ R1 .* ...
         (1 ./ (R3 + Za) + 1./ (R4 + Zb) + 1./ (Zc + 1 ./ (s .* C2)))) ./ ...
         (1 ./ (R4 + Zb) + (Zc + 1 ./ (s .* C2)) ./ R2 .* ...
         (1 ./ (R3 + Za) + 1./ (R4 + Zb) + 1./ (Zc + 1 ./ (s .* C2))));
    
    % Gain and its plot
    gain = abs(Hs);
    if(strcmp(command, 'SingleCalculation'))
        clf;
        loglog(freq, gain);
        hold on;
        grid on;
        xlabel('frequency (Hz)');
        ylabel('gain');
    end
    
    % find the range of gain
    max_gain = max(gain);
    min_gain = min(gain);
    
    % find the center frequency and its gain
    if (gamma > 0.5)         % band-stop filter
        fc = freq(find(gain == min_gain, 1, 'first'));
        Gc = min_gain;
    elseif (gamma < 0.5)     % band-pass filter
        fc = freq(find(gain == max_gain, 1, 'first'));
        Gc = max_gain;
    else                    % Hs = 1
        err = MException ('Gain:ConstantGain', ...
            'The gain is always 1 and the center frequency cannot be identified. Please set gama not equal to 0.5');
        throw(err);
    end
    
    % find the cutoff frequencies
    if(strcmp(command, 'SingleCalculation'))
        if ((gamma > 0.5 && max_gain >= (1 / sqrt(2))) || (gamma < 0.5 && max_gain <= sqrt(2)))
            errBound = MException('CutoffFrequency:NoCutoffFrequency',...
                'The gain is too small that there is no cutoff frequecies for this filter. Please change to other options.');
            throw(errBound);
        else
            cutoff_gain = max_gain / sqrt(2);
            if (gamma > 0.5)         % band-stop filter
                fc1 = freq(find(gain <= cutoff_gain, 1, 'first'));
                fc2 = freq(find(gain <= cutoff_gain, 1, 'last'));
            else                    % band-pass filter
                fc1 = freq(find(gain >= cutoff_gain, 1, 'first'));
                fc2 = freq(find(gain >= cutoff_gain, 1, 'last'));
            end
            if (isempty(fc1) || isempty(fc2)...       % the frequency range is not enough to cover the cutoff frequency
                    || fc1 == fc2 || fc1 == 0 || fc2 == maxfreq)    
                errBound = MException('CutoffFrequency:OutOfBound',...
                    'The frequency range is not enough to cover two cutoff frequencies. Please increase your frequency range');
                throw(errBound);
            end
        end
    else
        fc1 = 0;
        fc2 = 0;
    end
    
    % title on the plot
    str = sprintf('Center frequency f0 = %2f Hz\n Gain at f0 is G0 = %2f \n Cutoff frequencies fc1 = %2f Hz, fc2 = %2f Hz', ...
        fc, Gc, fc1, fc2);
    title(str);
end








