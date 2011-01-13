% Ex. 1
%
% QPSK
% Uncorrelated Rayleigh
% ML detection
%
% Beamforming simulation
% Assumes perfect CSI at transmitter and receiver

clear; 
close all;

Nsym = 1e6;  % Symbols
EsNo_db = 0:15; %2:10; 


snr = 10.^(EsNo_db/10);

s = signc(randcn(1, Nsym));
Es = 2;

lstr = {};

Nrx = 2;
n = randcn(Nrx, Nsym);

P_err_bf = zeros(4, length(snr));

for Ntx=[2 4]  
    H = randcn(Nrx, Ntx, Nsym);
    
    
    fprintf('\n %d> ', Ntx);
    
    for i=1:length(snr)
    fprintf('.');
    
        s_hat = zeros(1, Nsym);
        for k=1:Nsym
            % decompose H = u*d*v' or v'*H'*H*v = d^2
            % r = v'*H'*H*v + noise (so we need to transmit H*v and
            % detect by multiplying by v'*H'
            Hk = H(:,:,k);
            [u,d,v] = svd(Hk);
            % columns of v are the eigenvectors of H'H
            % first one correspponds to maximal eigenvalue
            % note that the columns are normalized to unit power
            w = v(:,1);       
            r = Hk*w*s(k) + n(:, k)/sqrt(snr(i)/Es);
            
            % Optimal receiver is matched to H*w
            s_hat(k) = (Hk*w)' * r;
        end
        s_hat = signc(s_hat);
        P_err_bf(Ntx, i) = mean(s_hat ~= s);
    end
end
   
EsNo_db_bf = EsNo_db;
save bf EsNo_db_bf P_err_bf

semilogy(EsNo_db_bf, P_err_bf(2,:), 'b--',...
    EsNo_db_bf, P_err_bf(4,:), 'k-.');

title('Eigen-beamforming');
ylabel('Pr(Symbol Error)')
xlabel('E_s/N_0 [dB]');
grid on;

