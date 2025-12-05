function R = pulse_fall_time_BP( t, y, thr_low, thr_high, interp_N, reference_levels )
%PULSE_FALL_TIME Calculate fall time of the pulse
%
%   Calculates fall time of the provided pulse. The method uses MATLAB function 'falltime'. By default, the method
%   returns 90% to 10% fall time. Alternatively, one can provide low and high thresholds. Finally, it is possible to use
%   FFT interpolation to increase accuracy of fall time estimation in case of low sampling densities. If the method is
%   called with single argument, then the time vector is assumed to be equal to the sample index.
%
%   R = pulse_fall_time( y )
%   R = pulse_fall_time( t, y )
%   R = pulse_fall_time( t, y, thr_low, thr_high )
%   R = pulse_fall_time( t, y, thr_low, thr_high, interp_N )
%   R = pulse_fall_time( t, y, thr_low, thr_high, interp_N, reference_levels )
%
%   Parameters:
%       t - time vector
%       y - pulse vector (must be same length as time)
%       thr_low - low threshold, specified as a fraction of the difference of state levels of the waveform (default: 0.1)
%       thr_high - high thresholds, specified as a fraction of the difference of state levels of the waveform (default: 0.9)
%       interp_N - interpolation factor (i.e. by how much sample density is increased)
%       reference_levels - low- and high-state levels for waveform (if empty, the function estimates them by itself)
%
%   Returns:
%       R - measured fall time of the pulse

% Default tolerance for state levels supplied to MATLAB's risetime function (in percents).
DEFAULT_LEVEL_TOLERANCE = 2;

% Test whether time vector is given
if nargin < 2
    % Only a single argument is given - we need to make our 't' vector, but since 't' is first argument (it must be so
    % to maintain compatibility with rest of the code), we first need to move its content into 'y'.
    y = t;
    t = 1:length(y);
end

% See if thresholds were provided
if nargin >= 4 && ~isempty(thr_low) && ~isempty(thr_high)
    % Use provided thresholds
    percent_of_levels = [thr_low thr_high] * 100;
else
    % Use default thresholds
    percent_of_levels = [10 90];
end

% Change data type of y to 'double', as it is required by the method calculating state levels
y = double(y);

% See if we are supposed to interpolate
if nargin >= 5 && ~isempty(interp_N) && interp_N > 1
    y = interpft(y, interp_N * length(y));
    t_step_interp = (t(2)-t(1))/interp_N;
    t = t(1) : t_step_interp : t(end) + t_step_interp * (interp_N-1);
end

% Check if reference levels were given
if nargin < 6 || isempty(reference_levels)
    % Estimate reference levels
    [~, ~, reference_levels] = pulse_find_peak(y);

    if reference_levels(1) > reference_levels(2)
        temp = reference_levels(1);
        reference_levels(1) = reference_levels(2);
        reference_levels(2) = temp;
    end

    if reference_levels(1) == reference_levels(2)
        reference_levels(2) = 3 * reference_levels(2);
    end

end


% Calculate fall time (ensure that we only measure one edge). If an edge is not found then increase state levels
% tolerance (up to lower crossing level). If we reach lower crossing level, then it means that the pulse was not found.
R = [];
ref_levels_tolerance = DEFAULT_LEVEL_TOLERANCE;

while length(R) < 1 && ref_levels_tolerance < percent_of_levels(1) %&& sum(reference_levels) > 0.01
    R = falltime(y, t, 'PercentReferenceLevels', percent_of_levels, 'StateLevels', reference_levels, 'Tolerance', ref_levels_tolerance);
    ref_levels_tolerance = ref_levels_tolerance + 1;

    % See if we have to increase the levels tolerance in order to find two mid-level crossing points
    if length(R) < 1
        warning('Could not estimate fall time, increasing state levels tolerance to %.1f %%', ref_levels_tolerance);
    end

    % See if multiple pulses were found
    if length(R) > 1
        warning('Multiple edges were found, returning multiple results');
    end
end

end

