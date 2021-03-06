function [transition_mat, emission_mat, bins, em_pdf] = run_hmm(s_seq, s_s, s_e, x)
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
[transition_mat, emission_mat] = hmmtrain(s_seq_ind, s_s,s_e,'maxiterations',1000,'tolerance',1e-15);
% smoothing filter coeffs. to obtain PDF plot
B = (0.3905/4)*[1 2 1];
A = [1 -0.9428 .3333];
% shifting bins to account for filtering delay
bins = bins+(1.1*bin_dist);
% smoothing the histogram
em_pdf = (filter(B, A, emission_mat'))';
% figure;
% hold on;
% legend_str = {};
% for i = 1:size(emission_mat,1)
%       %normalizing to PDF
%       em_pdf(i, :) = em_pdf(i, :)/(bin_dist*sum(em_pdf(i, :)));
%       % to get changing color for each State's PDF
%       color = (i-1)/size(emission_mat,1);
%       % plotting each State's PDF
%       plot(bins, em_pdf(i, :), 'Color', [sqrt(color), color, (color)^2]);
%       % concatenating to the legend string
%       legend_str{i} = strcat('State ', num2str(i));
% end
% grid;
% legend(legend_str);
% title('Estimated Prob. Dist. Functions of State Emission values');
% xlabel('Emission value');
% ylabel('PDF value');
end