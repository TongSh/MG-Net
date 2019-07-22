close all;
clear;
fs = 1e6;
bw = 125e3;
sf = 7;
Nchirp = 2^sf/bw*fs;
u = LoRaUtils(fs, bw, sf);

% read pkt
% datain = u.rfile('E:\LoRa Data\MG-Net\1tx2rx\pt_1');
% datain = u.ampCut(datain);
% [datain,cfo] = u.sync(datain);
% disp('Packet SYNC');
% sym = datain(14.25*u.chirp_n+1+(1:u.chirp_n));
% % sym(1:200) = 0;
% figure;
%     plot(u.chirp_t,real(sym));

sym = u.genSymbol(100);
sym = [zeros(1,101),sym(1:end-101)];
u.demodePlot(sym);

% FFT
fft_len = length(sym) * 100;
align_len = round(fft_len/(u.fs/u.bw));
freq = 0:u.bw/align_len:u.bw-u.bw/align_len;
z = fft(u.despread(sym),fft_len);
align1 = z(1:align_len);
align2 = z(end-align_len+1:end);

% Peak Phase
for i = 0
out = align1 + align2.*exp(1i*i);
[ma,I] = max(abs(out));
fprintf('peak phase = %.3f\n',angle(out(I)));
figure;hold on;
    plot(abs(align1),'--','Color','Blue');
    plot(abs(align2),'--','Color','m');
    plot(abs(out),'LineWidth',1.5,'Color','RED');
    xlim([I-1000 I+1000]);
    grid on
    grid minor
    box on
    xlabel('Freq');
    ylabel('abs. FFT');
    title(num2str(angle(out(I))));
end    


