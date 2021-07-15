function TLBT_DAE3_Stykel_hankel(problem)
%% Accumulated Hankel singular values plot for Stykel model
%
% Input:
% problem       model: 'stykel_small' or 'stykel_large'
%
% Author: Simon B�se (2021)
%
% See LICENSE.md

%% Initialization
% Set time steps. 999 is a placeholder for infinity.
time_steps = [999,5,2,0.5];
plot_legend = {'t_e = \infty','t_e = 5','t_e = 2','t_e = 0.5'};
len_steps = length(time_steps);
hsv = cell(len_steps,1);

%% Calculate Hankel singular values
for i = 1:len_steps   
    tlbt_opts = [];
    tlbt_opts.tf = time_steps(i);   
    if tlbt_opts.tf == 999; tlbt_opts.tl = 0; else; tlbt_opts.tl = 1; end
    hsv{i} = TLBT_DAE3_Stykel_mor(problem,tlbt_opts);  
end

%% Output
% Streamline data for plot.
lenghts = cellfun(@length, hsv);
max_length = max(lenghts);
hsv_collection = zeros(max_length,len_steps);
for i = 1:len_steps
    hsv_temp = hsv{i};
    temp_length = length(hsv_temp);
    hsv_collection(1:temp_length,i) = hsv_temp(:)'./hsv_temp(1);
end
hsv_collection(hsv_collection < 1e-13) = 0;

% Compile plot.
figure;
semilogy(hsv_collection,'LineWidth',1.5); 
title('Computed Hankel singular values');
xlabel('index');
ylabel('magnitude');
set(gca, 'XLimSpec', 'Tight');
legend(plot_legend);

end