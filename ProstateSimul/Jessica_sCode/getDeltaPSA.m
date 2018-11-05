function [delta_PSA] = getDeltaPSA(x, p)


T_plus = 0.01 * p(2);

TP = 0.006 * p(2) + 0.5;

T_minus = 0.006 * p(2) + 0.5;
delta_PSA = x(1)*T_plus + x(2)*TP + x(3)*T_minus;








%% To plot

% p = 1:100;
% 
% T_plus = 0.01 * p;
% 
% TP = 0.03 * p + 1;
% 
% T_minus = 0.03 * p + 1;
% 
% 
% figure
% hold on
%  
% plot(1:100, T_plus, 'k', 'LineWidth', 4)
% plot(1:100, TP, 'g', 'LineWidth', 4)
% plot(1:100, T_minus, '--r', 'LineWidth', 4)
% 
% legend('T+', 'TP', 'T-')
% 
% xlabel('Frequency of TP', 'FontSize', 16)
% ylabel('PSA produced per cell', 'FontSize', 16)




