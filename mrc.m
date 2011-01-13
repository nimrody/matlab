% Ex. 1
%
% QPSK
% Uncorrelated Rayleigh
% ML detection
%
% SISO and MRC 1xN simulation

clear; 
close all;

Nsym = 2e6;  % Symbols
EsNo_db = 0:40; 


% generate symbols (note: power = 2)
s = signc(randcn(1, Nsym));
Es = 2;

snr = 10.^(EsNo_db/10);

lstr = {}; % store legends

P_err_mrc = zeros(4, length(snr));

% b. MRC 1xN (1x1 [SISO], 1x2, 1x4) %%%%%%%%%%%%%%%%%%%%%%%%%%
for Nrx=[1 2 4]
    
    % Generate channel and unit power noise
    h = randcn(Nrx, Nsym);
    n = randcn(Nrx, Nsym);
    
    ss = repmat(s, Nrx, 1);
    
    for i=1:length(snr)
        r = ss .* h + n/sqrt(snr(i)/Es);
        
        % optimal matched receiver (each antenna x its channel*)
        y = sum(r .* conj(h), 1);
        
        s_hat = signc(y);
        P_err_mrc(Nrx,i) = mean(s_hat ~= s);
    end

end

EsNo_db_mrc = EsNo_db;

save mrc EsNo_db P_err_mrc EsNo_db_mrc

semilogy(EsNo_db_mrc, P_err_mrc(1,:), 'r-', ...
         EsNo_db_mrc, P_err_mrc(2,:), 'b--', ...
         EsNo_db_mrc, P_err_mrc(4,:), 'k-.');

title('MRC');
ylabel('Pr(Symbol Error)')
xlabel('E_s/N_0 [dB]');
legend('SISO', 'MRC 1x2', 'MRC 1x4');
grid on;