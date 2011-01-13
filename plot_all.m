clear;
close all;

mrc=load('mrc.mat');
stc=load('stc.mat');
sm =load('sm.mat');
bf =load('bf.mat');



semilogy(mrc.EsNo_db, mrc.P_err_mrc(2,:), 'r--', ...
         stc.EsNo_db, stc.P_err_stc(1,:), 'b-.');


title('STC 2x1 vs. MRC 1x2 (3 dB advantage)');
ylabel('Pr(Symbol Error)')
xlabel('E_s/N_0 [dB]');

legend('MRC 1x2', 'STC 2x1');
grid on;
axis([0 30 1e-4 1]);

figure;

semilogy(mrc.EsNo_db, mrc.P_err_mrc(4,:), 'r--', ...
         stc.EsNo_db, stc.P_err_stc(2,:), 'b-.');


title('STC 2x2 vs. MRC 1x4 (3 dB advantage)');
ylabel('Pr(Symbol Error)')
xlabel('E_s/N_0 [dB]');

legend('MRC 1x4', 'STC 2x2');
grid on;
axis([0 16 1e-4 1]);

figure;

load bf
semilogy(EsNo_db_bf, P_err_bf(2,:), 'b--',...
    EsNo_db_bf, P_err_bf(4,:), 'k-.');


title('Eigen-beamforming');
ylabel('Pr(Symbol Error)')
xlabel('E_s/N_0 [dB]');
legend('BF 2x2', 'BF 2x4');

grid on;

