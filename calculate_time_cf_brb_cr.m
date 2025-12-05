function [cf_time_ns, cr, amplitude_max, samples_intergral, cf_time_coarse] = calculate_time_cf_brb_cr(samples, correction_flag, cr_lut_data)
% params, samples - data with subtracted pedestal
%          time of first sample in ns

if ~exist('correction_flag','var')
    correction_flag = 1;
    
end

cfd = [];
cfd.delay = 1;
cfd.gain_delayed = 2;
time_base_ns = 8;
cr_lut_size = numel(cr_lut_data);
samples = samples(:);
% Create constant fraction waveforms by inverting, multiplying and shifting source waveforms
delayed = [-cfd.gain_delayed * samples(1:cfd.delay); -cfd.gain_delayed * samples;];
samples = [samples; samples(end-cfd.delay+1:end)];
cf_samples = samples + delayed;

% Find minimum and maximum of constant fraction waveform and deduct whether we have negative or positive pulses
[~, idx_min] = min(cf_samples);

idx = idx_min;

%results
cf_time_ns = NaN;
cr = NaN;
amplitude_max = max(samples);
samples_intergral = NaN;
cf_time_coarse = NaN;
cf_time_min = NaN;

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

    if (y_ratio < 1) && (y_ratio >= 0)
        % linear interpolation
        cf_time_vector = x1 + y_ratio;
        t0_pos = (cf_time_vector)*time_base_ns;
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
        if correction_flag
            cf_time_ns = t0_pos(end) + cr_lut_data(cr_idx);

        else
            cf_time_ns = t0_pos(end);
        end

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
    end
end

end