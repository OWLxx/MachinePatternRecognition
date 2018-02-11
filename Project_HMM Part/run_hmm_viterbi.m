function states = run_hmm_viterbi(s_seq, s_s, s_e, x)
% check if number of bins is an integer
if abs(x-round(x))~=0
      disp('Argument must be an integer! Stopping.');
      return;
end
nBins = x;
% s_e = ones(length(s_s),x)/x + randn(length(s_s),x)/x*0.1;
% obtain histogram for the sequence array
[freq, bins] = hist(s_seq, nBins);
% obtain the bin distance
bin_dist = (bins(2)-bins(1));
% fprintf('Bin width: %f\n', bin_dist);
quant_error = 0;
% indexing each float value to the nearest bin
for i = 1:length(s_seq)
      for j = 1:nBins
          dist = abs(bins - s_seq(i));
          [minval, indx] = min(dist);
          quant_error = quant_error + (s_seq(i) - bins(indx));
          s_seq_ind(i) = indx;
      end
end
% total quantization error in indexing
% fprintf('Total quantization error in indexing: %f\n', (quant_error));
% obtain HMMESTIMATE
states = hmmviterbi(s_seq_ind, s_s,s_e);
end