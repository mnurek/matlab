function [qf, ss, sig_resid] = calculate_qf_factor(signal,reference_signal)
    

    % scale samples to amp 2048
    scale_factor = round(2048/max(signal(1:4)));
    ss = int16(signal * scale_factor);
    % subtract reference from signal
    sig_resid = (ss - int16(reference_signal));
    % sum residuum 
    qf = int16(sum(abs(sig_resid))/numel(sig_resid));

end