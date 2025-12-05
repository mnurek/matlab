%%
load('C:\Users\mnurek\OneDrive - Politechnika Warszawska\HK\Testy\brb_spice\pulse_spice.mat')

waveform_length = length(pulse.time);


%%
current_time_s = pulse.time(1);
time_interval_s = 0;
samples_in_one_ns = [];
new_waveform = [];

for n = 1 : waveform_length
   current_time_s = pulse.time(n) - current_time_s;
   time_interval_s = time_interval_s + current_time_s;
   current_time_s = pulse.time(n);
    if(time_interval_s < 1e-9)
        samples_in_one_ns(end+1) = pulse.value(n);
    else
        new_waveform(end+1) = median(samples_in_one_ns);
        time_interval_s = 0;
        samples_in_one_ns = [];
        samples_in_one_ns(end+1) = pulse.value(n);        
    end
    
end

%%
time_ns = 0 : length(new_waveform) - 1;
try
    rise_time = pulse_rise_time_BP(time_ns, new_waveform, 0.1, 0.9, 1, [0 max(new_waveform)]);
    fall_time = pulse_fall_time_BP(time_ns, new_waveform, 0.1, 0.9, 1, [0 max(new_waveform)]);
    fwhm = pulse_fwhm_BP(time_ns, new_waveform, 1, [0 max(new_waveform)]);
catch
    rise_time = 0;
    fall_time = 0;
    fwhm = 0;
end