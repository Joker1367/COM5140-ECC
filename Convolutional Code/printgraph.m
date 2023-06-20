%% BSC hard decision
SNR = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6];
BER = [2.34E-01 2.07E-01 9.76E-02 6.21E-02 2.97E-02 1.30E-02 4.75E-03 2.21E-03 5.39E-04 1.38E-04 4.04E-05];
figure(1)
hold on 
grid on
plot(SNR, BER, '-O', 'color', "#0072BD", 'linewidth', 2);

% hold off
% set(gca, 'YScale', 'log')
% axis([1 6 1e-05 1])
% ylabel('BER');
% xlabel('Eb / N0(dB)');
% title('BSC under hard decision');
% legend('BER')

SNR = [1 1.5 2 2.5 3 3.5 4];
BER = [3.15E-02 1.70E-02 5.06E-03 1.74E-03 4.28E-04 8.52E-05 1.36E-05];
% figure(1)
% hold on 
% grid on
plot(SNR, BER, '-O', 'color', "#A2142F", 'linewidth', 2);

hold off
set(gca, 'YScale', 'log')
axis([1 6 1e-05 1])
ylabel('BER');
xlabel('Eb / N0(dB)');
title('128 Bit');
legend('hard decision', 'soft decision');

%% BSC hard decision
SNR = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7];
BER = [2.68E-01 2.37E-01 1.66E-01 1.14E-01 5.82E-02 2.89E-02 1.70E-02 6.05E-03 2.48E-03 8.74E-04 2.51E-04 7.71E-05 1.92E-05];
figure(1)
hold on 
grid on
plot(SNR, BER, '-O', 'color', "#0072BD", 'linewidth', 2);

% hold off
% set(gca, 'YScale', 'log')
% axis([1 6 1e-05 1])
% ylabel('BER');
% xlabel('Eb / N0(dB)');
% title('BSC under hard decision');
% legend('BER')

%% BSC soft decision
SNR = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6];
BER = [1.16E-01 6.61E-02 3.47E-02 1.57E-02 6.73E-03 3.13E-03 1.18E-03 5.26E-04 1.99E-04 7.73E-05 2.67E-05];
% figure(1)
% hold on 
% grid on
plot(SNR, BER, '-O', 'color', "#A2142F", 'linewidth', 2);

hold off
set(gca, 'YScale', 'log')
axis([1 7 1e-05 1])
ylabel('BER');
xlabel('Eb / N0(dB)');
title('Fixed Statee');
legend('hard decision', 'soft decision');
%% BSC hard decision
SNR = [1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5];
BER = [3.15E-01 2.27E-01 1.33E-01 7.19E-02 4.45E-02 1.93E-02 7.18E-03 3.15E-03 9.66E-4 3.27E-4 7.69E-05 1.90E-05];
figure(1)
hold on 
grid on
plot(SNR, BER, '-O', 'color', "#0072BD", 'linewidth', 2);

% hold off
% set(gca, 'YScale', 'log')
% axis([1 6 1e-05 1])
% ylabel('BER');
% xlabel('Eb / N0(dB)');
% title('BSC under hard decision');
% legend('BER')

%% BSC soft decision
SNR = [1 1.5 2 2.5 3 3.5 4 4.5];
BER = [7.65E-02 2.91E-02 1.03E-02 4.42E-03 1.35E-03 3.96E-04 8.03E-05 1.64E-05];
% figure(1)
% hold on 
% grid on
plot(SNR, BER, '-O', 'color', "#A2142F", 'linewidth', 2);

hold off
set(gca, 'YScale', 'log')
axis([1 6.5 1e-05 1])
ylabel('BER');
xlabel('Eb / N0(dB)');
title('Majority Vote');
legend('hard decision', 'soft decision');

