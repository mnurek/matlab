%
%
%
clearvars  -except HW
close all
data_path = 'C:\Users\Andrzej\Downloads\pulse_injection_sweep_1\02_07_2024_16-00\';
[filepath,name,ext] = fileparts(data_path);

ref_filename = 'pmt_sig_ref.csv';

aver_sig_ref = csvread(fullfile(filepath,ref_filename));


pe_amp_vect = [1];
timing_vect = zeros(1,numel(pe_amp_vect));
amp_vect_lsb = zeros(1,numel(pe_amp_vect));

result = struct();
    result  = struct();
    result.wave = [];
result.qf = [];
result.cr_idx = [];
for amp_idx  = 1 : numel(pe_amp_vect)
    filePattern = fullfile(data_path, strcat('pulse_amplitude_',num2str(pe_amp_vect(amp_idx)),'*.bin'));
    rawFiles = dir(filePattern);

    if numel(rawFiles) == 0
        fprintf('No raw data in directory %s.\n',data_path);
        return;
    end

    path = fullfile(data_path, rawFiles(1).name);
    disp(rawFiles(1).name)
    [ filepath , file_name , ext ] = fileparts( path ) ;
    datafile = wut_lab_daq_read_data_file(path);

    %%
    pre_sample = 4;
    post_sample = 31-pre_sample;
  

    for m = 1:18
        channel_data = datafile.channels{m};
        waveforms = reshape(channel_data, datafile.samplesperframe, []);
        waveforms(:,1:5) = [];
        waveforms = reshape(waveforms, [], 1);

        channel_data_pos = -waveforms;

        % equalize threshold
        channel_data_equ = channel_data_pos - median(channel_data_pos);
        waveforms = reshape(channel_data_equ, 512, []);

        signal_shape_cnt = zeros(1,33);
        for k = 1 : size(waveforms,2)
            samples = waveforms(:,k);

            if max(samples) > 2
                [cf_time, cr, amplitude_cr, amplitude_max, samples_intergral, cr_idx,  cf_time_coarse, cf_time_min] = calculate_time_cf_brb(samples);
                if cf_time_coarse > pre_sample
      
                    ss = samples(cf_time_coarse-pre_sample:cf_time_coarse+post_sample)';
                    reference_signal = aver_sig_ref(cr_idx,:);

                     [qf] = calculate_qf_factor(ss,reference_signal);

                    result.wave(end+1,:) = ss;
                    result.qf(end+1) = qf; 
                    result.cr_idx(end+1) = cr_idx; 
                    %                 % test_plot
                if qf > 120
                    figure(3)
                    clf
                    plot(ss,'k-')
                    hold on
                    plot(aver_sig_ref(cr_idx,:),'r-')
                    plot(sig_resid,'m--')
                    grid on
                    legend('Signal','Reference','Residual')
                    waitforbuttonpress
                end
                end
            end
        end
      
    end
end

ref_histogram_qf = figure(44);
histogram(result.qf)
grid on
xlabel('QF factor')
ylabel('Waveforms')
title('Histogram of quality factor for signals from pulse generator')
fig_name = 'ref_histogram_qf';
savefig(fullfile(filepath, strcat(fig_name, '.fig')))
exportgraphics(ref_histogram_qf, fullfile(filepath, strcat(fig_name, '.png')), 'Resolution',300);
