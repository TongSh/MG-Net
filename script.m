din = [10 22 24 42 37 77 89 22 63 9 0 0 0 0];
len_win = 5;
dout = movmean(din,len_win);
% the kth output uses the kth input as the center of the window
dout = dout(ceil(len_win/2):end);
disp(dout);
% Display the Result of Slide Average
figure; hold on
plot(len_win + (0:length(dout)-1),dout,'-o');
plot(1:length(din),din,'-o');
xlim([1 len_win+length(dout)-1]);

%%
% 单位脉冲信号
din = zeros(1,1e4);
din(5e3) = 1;
figure;plot(din);
% 脉冲信号包含所有频率的谐波分量
z = fft(din,100*length(din));
figure;
    plot(abs(z));
% 使用滑动平均滤波器对脉冲信号滤波    
len_win = 5;
dout = movmean(din,len_win);
% 对滤波后结果做傅里叶变换，在频域上观察滤波效果
z = fft(dout,100*length(dout));
hold on
    plot(abs(z));
% 调整滑动窗口的大小
len_win = 15;
dout = movmean(din,len_win);
z = fft(dout,100*length(dout));
    plot(abs(z));
xlabel('Frequency');
ylabel('abs. FFT');
grid on;
grid minor
box on
  