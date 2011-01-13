% Ex. 1
%
% QPSK
% Uncorrelated Rayleigh
% ML detection

% Spatial Multiplexing
%
%   Comparing the optimal ML decoder with the simple ZF decoder

clear;
close all;

Nsym = 1e6;  % Symbols
EsNo_db = 0:45;


snr = 10.^(EsNo_db/10);

Nst = 2; % number of streams
s = signc(randcn(Nst, Nsym));
Es = 2;

lstr = {};

Nrx = 2;
Ntx = 2;
n = randcn(Nrx, Nsym);

H = randcn(Nrx, Ntx, Nsym);


% generate all possible symbol combinations (4x4)
qpsk = [ 1+1j  1-1j  -1-1j -1+1j ]; % all 4 QPSK options

s_test = zeros(2, 0);
for i=1:length(qpsk)
    for j=1:length(qpsk)
        s_test(:,end+1) = [qpsk(i) qpsk(j)].';
    end
end


for i=1:length(snr)
    fprintf('.');
    
    s_hat_ml = zeros(Nst, Nsym);
    s_hat_zf = zeros(Nst, Nsym);
    d = zeros(1, 16);
    
    for k=1:Nsym
        Hk = H(:,:,k);
        r = Hk*s(:,k)/sqrt(Ntx) + n(:, k)/sqrt(snr(i)/Es);
        
        % Optimal ML receiver tests all posible s1, s2 
        % (compute the distance d -- then minimize it)
        n_hat = Hk*s_test/sqrt(Ntx) - repmat(r, 1, 16);
        d = sum(abs(n_hat).^2, 1);
        
        [dmin, sym_index] = min(d);
        s_hat_ml(:,k) = s_test(:,sym_index);
        
        % Suboptimal zero-forcing equalizer just inverts the channel
        s_hat_zf(:,k) = signc( pinv(Hk)*r );
        
    end
    P_err_zf(i) = mean(s_hat_zf(:) ~= s(:));
    P_err_ml(i) = mean(s_hat_ml(:) ~= s(:));
end

EsNo_db_sm = EsNo_db;

save sm P_err_ml P_err_zf EsNo_db_sm

semilogy(EsNo_db_sm, P_err_zf, 'r--', ...
         EsNo_db_sm, P_err_ml, 'b-.');


title('Spatial Multiplexing');
ylabel('Pr(Symbol Error)')
xlabel('E_s/N_0 [dB]');

legend('SM 2x2 (ML)', 'SM 2x2 (ZF)');
grid on;

