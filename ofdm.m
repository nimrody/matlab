clear;

EsNo_db = 0:2:15;

Nfft = 1024;
Nprefix = Nfft/8;   % guard interval, cyclic prefix

Nslots = 200;       % A slot is 1 preamble symbol + 25 data symbols
Nsym = 26;          % ofdm symbols per slot (including preamble)

Nguard = 100;       % empty carriers at each side

Es = 2; % qpsk symbols are +/-1

m=[];

% AWGN channel
h1=1;

% Multipath channel

h2 = [];
h2(2) = 1;
h2(2+9) = 0.9;


h = h1

filtstate = [];

for k=1:length(EsNo_db)
    fprintf('\n%d/%d ', k, length(EsNo_db));
    Nerr = 0; % symbol errors counter
    
    for i=1:Nslots
        if mod(i, 10)==0
            fprintf('.');
        end
        % generate data
        s = signc(randcn(Nfft-Nguard*2, Nsym));
        sg= [zeros(Nguard, Nsym)
            s
            zeros(Nguard, Nsym)];
        
        sg = circshift(sg, Nfft/2);
        
        % generate the transmitted signal - fft + adding cyclic prefix
        y = sqrt(Nfft)*ifft(sg);
        y1 = [y(end-Nprefix+1:end, :) ; y];
        y2 = y1(:).';
        
        % channel
        r = filter(h, 1, y2, filtstate);
        %%% m = [m  r]; % for spectrum estimation
        
        % add noise
        v = r + randcn(1, length(r)) / 10^(EsNo_db(k)/20) * sqrt(Es) * sqrt((Nfft-2*Nguard)/Nfft);
        
        % discard cyclic prefix and go back to freq domain
        vv = reshape(v, Nfft+Nprefix, Nsym);
        sg_hat1 = 1/sqrt(Nfft)*fft(vv(Nprefix+1:end,:));
        
        % discard guard carriers
        sg_hat2 = circshift(sg_hat1, Nfft/2-Nguard);
        s_hat1  = sg_hat2(1:Nfft-2*Nguard,:);
        
        % estimate channel and invert channel effect
        h_est = s_hat1(:, 1)  ./ s(:, 1);
        s_hat =  s_hat1 ./ repmat(h_est, 1, Nsym);
        
        % count errors only in the data (ignore preamble symbol)
        Nerr = Nerr + sum(sum( signc(s_hat(:,2:end)) ~= s(:,2:end) ));
    end
    Perr(k) = Nerr / (Nslots*(Nsym-1)*(Nfft-Nguard*2));
end

%%%estimator = spectrum.welch;
%%%psd(estimator, m, 'CenterDC', true);

figure;
semilogy(EsNo_db, Perr);
