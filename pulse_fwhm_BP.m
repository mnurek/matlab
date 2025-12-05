function fwhm = pulse_fwhm_BP(t, y, interp_N, reference_levels)
%PULSE_FWHM Calculate FWHM of the pulse
%
%   Calculates FWHM (Full Width Half Maximum) of the provided pulse. The method uses MATLAB's 'midcross' function. 
%   It is possible to use FFT interpolation to increase accuracy of FWHM estimation in case of low sampling densities. 
%   If the method is called with single argument, then the time vector is assumed to be equal to the sample index.
%   The pulse can be of either positive or negative polarity.
%
%   fwhm = pulse_fwhm( y ) 
%   fwhm = pulse_fwhm( t, y ) 
%   fwhm = pulse_fwhm( t, y, interp_N ) 
%   fwhm = pulse_fwhm( t, y, interp_N, reference_levels )
%
%   Arguments:
%       t - time vector
%       y - pulse vector (must be same length as time)
%       interp_N - interpolation factor (i.e. by how much sample density is increased)
%       reference_levels - low- and high-state levels for waveform (if empty, the function estimates them by itself)
%
%   Returns:
%       fwhm - measured FWHM of the pulse.

% Default tolerance for state levels supplied to MATLAB's midcross function (in percents).
DEFAULT_LEVEL_TOLERANCE = 2;

% Test whether time vector is given
if nargin < 2
    % Only a single argument is given - we need to make our 't' vector, but since 't' is first argument (it must be so
    % to maintain compatibility with rest of the code), we first need to move its content into 'y'.
    y = t;
    t = 1:length(y);
end

% Change data type of y to 'double', as it is required by the method calculating state levels
y = double(y);

% See if we are supposed to interpolate
if nargin >= 3 && ~isempty(interp_N) && interp_N > 1
    y = interpft(y, interp_N * length(y));
    t = 1:length(y);
    t = t ./ interp_N;
end

% Check if reference levels were given
if nargin < 4 || isempty(reference_levels)
    % Estimate reference levels ignoring tail of the pulse
    [~, ~, reference_levels] = pulse_find_peak(y);

    if reference_levels(1) > reference_levels(2)
        temp = reference_levels(1);
        reference_levels(1) = reference_levels(2);
        reference_levels(2) = temp;
    end

    if reference_levels(1) == reference_levels(2)
        reference_levels(2) = 4 * reference_levels(2);
    end
end

% Measure pulse FWHM
C = [];
ref_levels_tolerance = DEFAULT_LEVEL_TOLERANCE;
while length(C) < 2 && ref_levels_tolerance < 50   % We run only till 50 % tolerance - further running is pointless
    C = midcross(y, t, 'StateLevels', reference_levels, 'Tolerance', ref_levels_tolerance);
    ref_levels_tolerance = ref_levels_tolerance + 1;
    
    % See if we have to increase the levels tolerance in order to find two mid-level crossing points
    if length(C) < 2
        warning('Insufficient mid-level crossing points, increasing state levels tolerance to %.1f %%', ref_levels_tolerance);
    end
    
    % See if more than two edges were found
    if length(C) > 2
        warning('More than two edges were found, estimate may be inaccurate');
    end
end

fwhm = C(2) - C(1);

end

