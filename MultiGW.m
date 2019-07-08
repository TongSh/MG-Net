close all;
clear;
fs = 1e6;
bw = 125e3;
sf = 12;
Nchirp = 2^sf/bw*fs;
u = LoRaUtils(fs, bw, sf);

tx_to = round(u.chirp_n * 11.45);
path_to = round(fs * 500 / 3e8);

% load two packets at G1
p1 = u.rfile('E:\DataSet\multi_gw\SF12\round10\G1.cfile');
p1 = u.ampCut(p1);
p1 = u.sync(p1);
u.mSpectrogram(p1);
p1 = p1(round(12.25*Nchirp):end);

% load two packets at G2
p2 = u.rfile('E:\DataSet\multi_gw\SF12\round10\G2.cfile');
p2 = u.ampCut(p2);
p2 = u.sync(p2);
u.mSpectrogram(p2);
p2 = p2(round(12.25*Nchirp):end);

pkt_len = floor(size(p1,2) / Nchirp);
fid = fopen('E:\DataSet\multi_gw\SF12\round10\G1.csv','w');
fprintf(fid,'%s\n','win,peak,freq,bin,value');
for base = 0:pkt_len-1
    sig = p1(base*(u.chirp_n)+(1:u.chirp_n));
% 	u.demodePlot(sig);
    [peak_n, heights, freqs, Amps] = peakSearch(u,sig,0);
    bins = freqs/bw * 2^sf;
    fprintf('-----------Window %d--------------------\n',base+1);
    for i = 1:peak_n
        fprintf('height = %g, freq = %g [value = %g]\n',heights(i),freqs(i),bins(i));
        fprintf(fid,'%s\n',[num2str(base+1),',',num2str(heights(i)),',',num2str(freqs(i)),',',num2str(bins(i)),',',num2str(round(bins(i)))]);
    end
    fprintf(fid,'%s\n','');
end

pkt_len = floor(size(p2,2) / Nchirp);
fid = fopen('E:\DataSet\multi_gw\SF12\round10\G2.csv','w');
fprintf(fid,'%s\n','win,peak,freq,bin,value');
for base = 0:pkt_len-1
    sig = p2(base*(u.chirp_n)+(1:u.chirp_n));
% 	u.demodePlot(sig);
    [peak_n, heights, freqs, Amps] = peakSearch(u,sig,0);
    bins = freqs/bw * 2^sf;
    fprintf('-----------Window %d--------------------\n',base+1);
    for i = 1:peak_n
        fprintf('height = %g, freq = %g [value = %g]\n',heights(i),freqs(i),bins(i));
        fprintf(fid,'%s\n',[num2str(base+1),',',num2str(heights(i)),',',num2str(freqs(i)),',',num2str(bins(i)),',',num2str(round(bins(i)))]);
    end
    fprintf(fid,'%s\n','');
end


   
    
