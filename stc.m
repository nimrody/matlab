% Ex. 1
%
% QPSK
% Uncorrelated Rayleigh
% ML detection

clear; 
close all;

Nsym = 1e6;  % Symbols
EsNo_db = 0:30; 


snr = 10.^(EsNo_db/10);

s = signc( randcn(1, Nsym) );
Es = 2;

% c1. STC 2x1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ntx=2; 
Nrx=1;

h = randcn(Ntx, Nsym);
h(:,2:2:end) = h(:,1:2:end); % STC requires constant channel over two symbols

n = randcn(Nrx, Nsym);

% Alamouti encoding   [s1   -conj(s2)]
%                     [s2    conj(s1)]

% Decoding: h1, h2 - channel from tx1 and tx2
%           y1, y2 - received samples at t1 and t2
%   s1_hat = conj(h1)*y1 + h2*conj(y2)
%   s2_hat = conj(h2)*y1 - h1*conj(y2)
%   

s1(1:2:Nsym) = s(1:2:Nsym);
s1(2:2:Nsym) = -conj(s(2:2:Nsym));

s2(1:2:Nsym) = s(2:2:Nsym);
s2(2:2:Nsym) =  conj(s(1:2:Nsym));

ss = [s1 ; s2] / sqrt(2);  % power divided by 2

s_hat = zeros(1, Nsym);
lstr = {};

for i=1:length(snr)
    r = sum(ss .* h) + n/sqrt(snr(i)/Es);

    % Alamouti decoding
    s_hat(1:2:end) = conj(h(1, 1:2:end)).*r(1:2:end) + (h(2, 1:2:end)).*conj(r(2:2:end));
    s_hat(2:2:end) = conj(h(2, 1:2:end)).*r(1:2:end) - (h(1, 1:2:end)).*conj(r(2:2:end));

    s_hat = signc(s_hat);
    P_err_stc(1,i) = mean(s_hat ~= s);
end
 



% c2. STC 2x2

% Same trasnmit structure, so 'ss' remains unchanged
% Decoding (hij - from antenna j to i)
% [h11' h21'   h12 h22 ][y11]    \ t1
% [h12' h22'  -h11 -h21][y12]    / 
%                       [y21*]   \ t2
%                       [y22*]   /
%                       

h11 = randcn(1, Nsym/2);
h12 = randcn(1, Nsym/2);
h21 = randcn(1, Nsym/2);
h22 = randcn(1, Nsym/2);

n1 =  randcn(1, Nsym);
n2 =  randcn(1, Nsym);


for i=1:length(snr)
    r1 = sum(ss .* kron([h11; h12], [1 1])) + n1/sqrt(snr(i)/Es);
    r2 = sum(ss .* kron([h21; h22], [1 1])) + n2/sqrt(snr(i)/Es);
    
    y11 = r1(1:2:end);
    y12 = r2(1:2:end);
    y21 = r1(2:2:end);
    y22 = r2(2:2:end);
    
    % Alamouti decoding
    s_hat(1:2:end) =  conj(h11).*y11 + conj(h21).*y12 + h12.*conj(y21) + h22.*conj(y22);
    s_hat(2:2:end) =  conj(h12).*y11 + conj(h22).*y12 - h11.*conj(y21) - h21.*conj(y22);

    s_hat = signc(s_hat);
    P_err_stc(2,i) = mean(s_hat ~= s);
end

EsNo_db_stc = EsNo_db;

save stc EsNo_db P_err_stc EsNo_db_stc

semilogy(EsNo_db_stc, P_err_stc(1,:), 'r--', ...
         EsNo_db_stc, P_err_stc(2,:), 'k-.');

legend('STC 2x1', 'STC 2x2');
     
title('STC (Alamouti)')
ylabel('Pr(Symbol Error)')
xlabel('E_s/N_0 [dB]');


grid on;
