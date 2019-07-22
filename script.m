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
% ��λ�����ź�
din = zeros(1,1e4);
din(5e3) = 1;
figure;plot(din);
% �����źŰ�������Ƶ�ʵ�г������
z = fft(din,100*length(din));
figure;
    plot(abs(z));
% ʹ�û���ƽ���˲����������ź��˲�    
len_win = 5;
dout = movmean(din,len_win);
% ���˲�����������Ҷ�任����Ƶ���Ϲ۲��˲�Ч��
z = fft(dout,100*length(dout));
hold on
    plot(abs(z));
% �����������ڵĴ�С
len_win = 15;
dout = movmean(din,len_win);
z = fft(dout,100*length(dout));
    plot(abs(z));
xlabel('Frequency');
ylabel('abs. FFT');
grid on;
grid minor
box on
  