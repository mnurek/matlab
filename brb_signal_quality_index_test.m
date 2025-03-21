%
%
%
clearvars  -except HW
close all
data_path = 'C:\Users\Andrzej\Downloads\pulse_injection_sweep_1\02_07_2024_16-00\nom_hv__15_04_2024_13-03_1_20240415130524.bin';
[filepath,name,ext] = fileparts(data_path);
grid minor
ref_filename = 'pmt_sig_ref.csv';
aver_sig_ref = csvread(fullfile(filepath,ref_filename));



%% load data
datafile = wut_lab_daq_read_data_file(fullfile(data_path));

%% calculate test data
pre_sample = 4;
post_sample = 31-pre_sample;

result  = struct();
result.wave = [];
result.qf = [];
result.cr_idx = [];
result.raw_samples = [];
result.cf_time_coarse = [];
for m = 1:18
    channel_data = datafile.channels{m};
    waveforms = reshape(channel_data, datafile.samplesperframe, []);
    waveforms(:,1:5) = [];
    waveforms = reshape(waveforms, [], 1);

    channel_data_pos = waveforms;

    % equalize threshold
    channel_data_equ = channel_data_pos - median(channel_data_pos);
    waveforms = reshape(channel_data_equ, 512, []);

    signal_shape_cnt = zeros(1,33);
    for k = 1 : 1e3 %size(waveforms,2)
        samples = waveforms(:,k);
        
        if max(samples) > 2
            [cf_time, cr, amplitude_cr, amplitude_max, samples_intergral, cr_idx,  cf_time_coarse, cf_time_min] = calculate_time_cf_brb(samples);
            if cf_time_coarse > pre_sample && (cf_time_coarse < (512-post_sample))
                
                result.raw_samples(end+1,:) = samples;
                ss = samples(cf_time_coarse-pre_sample:cf_time_coarse+post_sample)';
                reference_signal = aver_sig_ref(cr_idx,:);  % reference signal is now 16 samples long
                [qf] = calculate_qf_factor(ss,reference_signal);

                result.wave(end+1,:) = ss;
                result.qf(end+1) = qf;
                result.cr_idx(end+1) = cr_idx;
                result.cf_time_coarse(end+1) = cf_time_coarse;
   
%                 % test_plot
%                 if qf > 60
%                     figure(3)
%                     clf
%                     plot(ss,'k-')
%                     hold on
%                     plot(aver_sig_ref(cr_idx,:),'r-')
%                     plot(sig_resid,'m--')
%                     grid on
%                     legend('Signal','Reference','Residual')
%                     waitforbuttonpress
%                 end
            end
        end
    end

end

%% plot and save data

disp(strcat('Saving data to: ',filepath))
csvwrite(fullfile(filepath,'qf_vect.csv'),result.qf')
csvwrite(fullfile(filepath,'samples.csv'),result.wave)
csvwrite(fullfile(filepath,'raw_samples.csv'),result.raw_samples)
csvwrite(fullfile(filepath,'cr_idx.csv'),result.cr_idx')
csvwrite(fullfile(filepath,'cf_time_coarse.csv'),result.cf_time_coarse')


histogram_qf = figure(44);
histogram(result.qf)
title('Histogram of quality factor for photon signals')
grid on
xlabel('QF factor')
ylabel('Waveforms')
fig_name = 'photon_sig_histogram_qf';
savefig(fullfile(filepath, strcat(fig_name, '.fig')))
exportgraphics(histogram_qf, fullfile(filepath, strcat(fig_name, '.png')), 'Resolution',300);


qf_vector = 40:20:180;
waves = result.wave(:,:);
qf = result.qf(:);
g_sam_fig = figure(333);
sgtitle('QF samples')
for qf_idx = 1 : numel(qf_vector) 
    disp(qf_idx)
    subplot(3,3,qf_idx)
    plot(waves((qf>qf_vector(qf_idx)) & (qf<(qf_vector(qf_idx)+20)),:)', 'k-')
    title(sprintf('Samples with QF > %d & QF < %d',qf_vector(qf_idx),qf_vector(qf_idx)+20))
    pause(0.1)
    grid on
end
fig_name = 'good_samples';
%savefig(fullfile(filepath, strcat(fig_name, '.fig')))
exportgraphics(g_sam_fig, fullfile(filepath, strcat(fig_name, '.png')), 'Resolution',300);



