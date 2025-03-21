function [cf_time_ns, cr, amplitude_cr, amplitude_max, samples_intergral, cr_idx, cf_time_coarse, cf_time_min] = calculate_time_cf_brb(samples)
% params, samples - data with subtracted pedestal
%          time of first sample in ns
cfd = [];
cfd.delay = 1;
cfd.gain_delayed = 2;
time_base_ns = 8;
cr_lut_size = 33;
cr_lut_data = [0.00,0.07,0.17,0.29,0.38,0.51,0.65,0.79,0.91,1.04,1.14,1.25,1.33,1.38,1.45,1.46,1.50,1.50,1.49,1.44,1.39,1.33,1.27,1.18,1.09,0.97,0.88,0.73,0.62,0.49,0.32,0.16,0.00];
cr_amp_lut_data = [1.03,1.04,1.05,1.06,1.08,1.10,1.13,1.15,1.18,1.22,1.25,1.29,1.33,1.37,1.41,1.45,1.49,1.53,1.57,1.61,1.65,1.68,1.72,1.76,1.79,1.83,1.86,1.89,1.92,1.95,1.98,2.00,2.02];
% Create constant fraction waveforms by inverting, multiplying and shifting source waveforms
delayed = [-cfd.gain_delayed * samples(1:cfd.delay); -cfd.gain_delayed * samples;];
samples = [samples; samples(end-cfd.delay+1:end)];
cf_samples = samples + delayed;

% Find minimum and maximum of constant fraction waveform and deduct whether we have negative or positive pulses
[~, idx_min] = min(cf_samples);
idx = idx_min;

if idx > 1
    % Find first point  that is above zero
    while true
        % We did - decrease its index by one (i.e. move left)
        idx = idx - 1;
        if (idx == 1)
            break;
        end
        if (cf_samples(idx) >= 0) 
            break;
        end
    end

    % Use linear interpolation to get sub-sample index of the zero-crossing point
    x1 = idx;
    x2 = x1+1;
    y1 = cf_samples(sub2ind(size(cf_samples), x1, 1:size(cf_samples,2)));
    y2 = cf_samples(sub2ind(size(cf_samples), x2, 1:size(cf_samples,2)));
    y_ratio = (double(y1) ./ double(y1-y2));
    % linear interpolation
    cf_time_vector = x1 + y_ratio;
    t0_pos = (cf_time_vector-1)*time_base_ns;
    cf_time_coarse = x1;
    cr = y_ratio;
    
    % perform time correction accordingly to cr value
    cr_idx = round(cr*(cr_lut_size-1))+1;
    if cr_idx <= 0
        cr_idx = 1;
    end
    if cr_idx > cr_lut_size
        cr_idx = cr_lut_size;
    end
    cf_time_ns = t0_pos(1) + cr_lut_data(cr_idx);
    cf_time_min = cf_time_ns/8.0 - cf_time_coarse;
    amplitude_cr = samples(x1) * cr_amp_lut_data(cr_idx);
    amplitude_max = max(samples);
   
    % integration start from first sample
    samples_intergral = 0;
    for idx = 1:numel(samples)
        % falling edge - check if below thr
        if idx > x2 && samples(idx) < 0
            break
        end
        samples_intergral = samples_intergral + samples(idx);
    end

else
    cf_time_ns = NaN;
    cr = NaN;
    amplitude_cr = NaN;
    amplitude_max = NaN;
    samples_intergral = NaN;
end

end