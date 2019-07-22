close all;
clear;
fs = 1e6;
bw = 125e3;
sf = 7;
Nchirp = 2^sf/bw*fs;
u = LoRaUtils(fs, bw, sf);

%% Packet Generation
% base_up = u.genSymbol(0);
data = u.genPacket(round(rand(1,40)*2^sf));
% data = u.genPacket(zeros(1,40));
disp('Packet Generated!');
dataraw = data;

%% Write Files
dataout = [zeros(1,u.chirp_n*20),data];
u.mSpectrogram(dataout);
u.wfile('E:\LoRa Data\MG-Net\genPkt\upchirp2',dataout);
fclose all;

%% Read Files
datain = u.rfile('E:\LoRa Data\MG-Net\1tx2rx\pt_1');
datain = u.ampCut(datain);
[datain,cfo] = u.sync(datain);
disp('Packet SYNC');    

%% Delta Phase
close all;
p = zeros(2,1);
figure;
for file = 1:2
    datain = u.rfile(['E:\LoRa Data\MG-Net\1tx2rx\pt_' num2str(file)]);
    % datain = datain(5e5:end);
    datain = u.ampCut(datain);
    [datain,cfo] = u.sync(datain);
    disp('Packet SYNC');

    % Despread and Decode
    a = angle(u.genSymbol(0));

    for i = 12.25:50.25
        sig = datain(round(u.chirp_n*i+0)+(1:u.chirp_n));
        % sig = sig.*exp(1i*2*pi*(-cfo)*u.chirp_t);
        fft_len = length(sig) * 100;
        align_len = round(fft_len/(u.fs/u.bw));
        freq = 0:u.bw/align_len:u.bw-u.bw/align_len;
        z = fft(u.despread(sig),fft_len);
        % 分别提取两段频率处的输出结果
        align1 = z(1:align_len);
        align2 = z(end-align_len+1:end);
        sum_align = align1 + align2;
        [ma,I] = max(abs(sum_align));
        fprintf('【peak%d】phase = %.3f, height = %.3f, I = %d\n',i-11.25,angle(sum_align(I)),abs(sum_align(I)),I);
        idx = round(i - 11.25);
        p(file,idx) = angle(sum_align(I));
%         while idx > 1 && p(file,idx) > p(file,idx-1)
%             p(file,idx) = p(file,idx) - 2*pi;
%         end
    end
    hold on;plot(p(file,:),'-o');
end
grid on
grid minor
box on
xlabel('Symbols');
ylabel('Phase (radian)');
title('两个RX收同一数据包（波峰相位）');

d = mod(p(1,:) - p(2,:),2*pi);
figure;
hold on
plot(d);
grid on
grid minor
box on
xlabel('Symbols');
ylabel('Delta Phase (radian)');
title('波峰相位差');

    
%% 两段信号的相位差，与编码的内容是否有关
close all;
sym = u.genSymbol(50);
% for i = 1:50
% sym = [zeros(1,i),sym(1:end-i)];
sym(1:80) = 0;
u.demodePlot(sym);
% FFT
fft_len = length(sym) * 100;
align_len = round(fft_len/(u.fs/u.bw));
freq = 0:u.bw/align_len:u.bw-u.bw/align_len;
z = fft(u.despread(sym),fft_len);
align1 = z(1:align_len);
align2 = z(end-align_len+1:end);
% Peak Phase
[~,I] = max(abs(align1));
p1 = angle(align1(I));
[~,I] = max(abs(align2));
p2 = angle(align2(I));
fprintf('p1 = %.3f, p2 = %.3f, delta = %.3f\n',p1,p2,p1-p2);

out = align1 + align2;
figure;hold on
plot(abs(align1),'--');
plot(abs(align2),'--');
plot(abs(out),'LineWidth',1.5);
xlim([I-1000 I+1000]);

grid on
grid minor
box on
xlabel('Frequency');
ylabel('abs. FFT');
% end

%% 噪声通过什么方式影响波峰测量
u = LoRaUtils(1e6, 125e3, 12);
sym = u.genSymbol(30);
sym = awgn(sym,-30);
% FFT
fft_len = length(sym) * 100;
align_len = round(fft_len/(u.fs/u.bw));
freq = 0:u.bw/align_len:u.bw-u.bw/align_len;
z = fft(u.despread(sym),fft_len);
align1 = z(1:align_len);
align2 = z(end-align_len+1:end);
% Peak Phase
out = align1 + align2;
[~,I] = max(abs(out));
figure;hold on
plot(abs(align1),'--');
plot(abs(align2),'--');
plot(abs(out),'LineWidth',1.5);
% xlim([I-1000 I+1000]);
grid on
grid minor
box on
xlabel('Frequency');
ylabel('abs. FFT');






   
    
