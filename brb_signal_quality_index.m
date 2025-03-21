%
%
%
clearvars  -except HW
close all
data_path = 'C:\Users\Andrzej\Downloads\pulse_injection_sweep_1\02_07_2024_16-00\';



pe_amp_vect = [3];
timing_vect = zeros(1,numel(pe_amp_vect));
amp_vect_lsb = zeros(1,numel(pe_amp_vect));

result = struct();

result.amp_vs_cr = [];

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

    %% saved signal reference window length is 32 (4 samples before x1 )
    pre_sample = 4;
    post_sample = 31-pre_sample;
    signal_shape = [];
    result  = struct();
    for m = 1:18
        channel_data = datafile.channels{m};
        waveforms = reshape(channel_data, datafile.samplesperframe, []);
        waveforms(:,1:5) = [];
        waveforms = reshape(waveforms, [], 1);

        channel_data_pos = -waveforms;
        [n, edges] = histcounts(channel_data_pos,'BinWidth',1);

        bin_c = edges(1:end-1)+0.5;
        pedestal_fit = fit(bin_c.',n.','gauss1');

        % equalize threshold
        channel_data_equ = channel_data_pos - median(channel_data_pos);
        waveforms = reshape(channel_data_equ, 512, []);
        cr_vect = [];
        amp_vect = [];
        time_ns = [];
        sampl_vect = [];
        time_vect = 8*(1:512);


        signal_shape_cnt = zeros(1,33);
        for k = 1 : size(waveforms,2)
            samples = waveforms(:,k);

            if max(samples) > 2
                [cf_time, cr, amplitude_cr, amplitude_max, samples_intergral, cr_idx,  cf_time_coarse, cf_time_min] = calculate_time_cf_brb(samples);
                if cf_time_coarse > pre_sample
                    cr_vect(k) = cr;
                    amp_vect(m,k) = amplitude_cr;
                    ss = samples(cf_time_coarse-pre_sample:cf_time_coarse+post_sample);
                    tt = time_vect - cf_time;
                    tt = tt(cf_time_coarse-pre_sample:cf_time_coarse+post_sample);
                    time_ns = [time_ns; tt'];
                    sampl_vect = [sampl_vect; ss];

                    cr_idx_vect(m,k) = cr_idx;
                    signal_shape{m}(k,:) = ss;
                end
            end
        end
        %% plot data
        figure(2323)
        clf
        pmt_sig_ref =[];
        for cr_idx = 1: 33
            sampl = signal_shape{m}(cr_idx_vect(m,:)==cr_idx,:);
            samp_aver = mean(sampl);
            pmt_sig_ref(cr_idx,:) = samp_aver;
            plot((1:33:32*33)- cr_idx,samp_aver,'.')
            hold on
        end
        result.sig_ref{m} = pmt_sig_ref;

        figure(232)
        plot(time_ns,sampl_vect,'.')
    end
end


%% plot references
figure(3434)
clf
 aver_sig_ref= zeros(size(result.sig_ref{1},1),size(result.sig_ref{1},2));
%calc averaged response
for m = 1:18
    aver_sig_ref = aver_sig_ref + result.sig_ref{m};
    plot(result.sig_ref{m}','k.')
    hold on
end
aver_sig_ref = round(aver_sig_ref./18);

for cr_idx = 1: 33
    % subtract pedestal
    aver_sig_ref(cr_idx,:) =   aver_sig_ref(cr_idx,:) -  ones(1,size(aver_sig_ref,2))*aver_sig_ref(cr_idx,1); 
    ss =  aver_sig_ref(cr_idx,:);
    % scale to fixed amplitude 2048
     aver_sig_ref(cr_idx,:)= round(2048*ss/max(ss));
end


plot(aver_sig_ref','mx')
grid minor
csvwrite(fullfile(data_path,'pmt_sig_ref.csv'),aver_sig_ref)
