%% USAGE
% sf = 12;
% data = [2540 1152 672 2396 1188 3508 40 3088 3236 3916 2728 2764 1416 2832 1388 800 3196 344 960 744 3100 296 1892 468 2699 2263 2061 2277 2164 2157 1578 2453];
% d = LoRaDecoder(fs, bw, sf);
% disp(d.demodulate(data)); % CODE: 09 90 40 01 02 03 04 05 06 07 08 09 ad 01

%%
% Matlab defines two types of class: Value Class and Handle Class.
% Value Class does not support modifing class variables in class methods.
classdef LoRaUtils < handle
    %M_ZIG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fs                      % Sampling frequency
        bw                      % Bandwidth (125kHz 250kHz 500kHz) 
        sf                      % Spreading factor (6-12)
        chirp_n                 % Number of sampling points 
        chirp_len               % Duration length
        chirp_t                 % Time index
    end
    
    methods
        function obj = LoRaUtils(fs, bw, sf)
            %M_ZIG Construct an instance of this class
            %   Detailed explanation goes here
            obj.fs = fs;
            obj.bw = bw;
            obj.sf = sf;
            obj.chirp_len = 2^sf/bw;
            obj.chirp_n = fs*obj.chirp_len;
            obj.chirp_t = 0:1/fs:obj.chirp_len-1/fs;
        end
        
        function [freq,align_sig] = loraFFT(obj, sig, len_fft)
            %LORAFFT FFT on the despreaded signal
            %   step1:FFT with Nchirp points
            %   step2:alias to 0~BW
            if nargin < 3 || isempty(len_fft)
                len_fft = size(sig,2)*10;
            end
            z = abs(fft(sig,len_fft));
            align_n = obj.fs/obj.bw;
            align_len = round(len_fft/align_n);
            freq = 0:obj.bw/align_len:obj.bw-obj.bw/align_len;
            % frequency aligning
            align_sig = z(1:align_len);
            for j = 1:align_n-1
                align_sig(1,:) = align_sig(1,:) + z(j*align_len + (1:align_len));
            end
            align_sig = abs(align_sig);
        end
              
        function sig = genSymbol(obj, k, is_down)
            %LORASYMBOL generate lora symbol
            %   Detailed explanation goes here
            if nargin < 3 || isempty(is_down) || is_down == 0
                f0 = -obj.bw/2; % start freq
                f1 = obj.bw/2;  % end freq
            else
                f0 = obj.bw/2;  % start freq
                f1 = -obj.bw/2; % end freq
            end    
            chirpI = chirp(obj.chirp_t, f0, obj.chirp_len, f1, 'linear', 90);
            chirpQ = chirp(obj.chirp_t, f0, obj.chirp_len, f1, 'linear', 0);
            baseline = complex(chirpI, chirpQ);
            baseline = repmat(baseline,1,2);
            clear chirpI chirpQ
            
            offset = round((2^obj.sf - k) / 2^obj.sf * obj.chirp_n);
            sig = baseline(offset+(1:obj.chirp_n));
        end
        
        function spectrum(obj,sig)
            %SPECTRUM plot the spectrum of a chirp signal
            %   Detailed explanation goes here
            x = 0:1/obj.sf:2^obj.sf/obj.bw;
            y = -obj.bw/2:obj.bw/2;
            win_length = 2^(obj.sf-2);
            s = spectrogram(sig,win_length,win_length-1,1024);
            data_len = 64;
            b = abs(s(data_len:-1:1,:));
            c = abs(s(end:-1:end-data_len+1,:));
            d = [b;c];
            %     figure;
                imagesc(x,y,d);
                title('Spectrogram');
                xlabel('Time');
                ylabel('Frequency');
        end
        
        function d = mSpectrogram(obj,sig)
            %SPECTRUM plot the spectrum of a chirp signal
            %   Detailed explanation goes here
            x = 0:1/obj.sf:2^obj.sf/obj.bw;
            y = -obj.bw/2:obj.bw/2;
            win_length = 2^(obj.sf-2);
            s = spectrogram(sig,win_length,round(win_length*0.7),1024);
            valid_data_len = 64;
            b = abs(s(valid_data_len:-1:1,:));
            c = abs(s(end:-1:end-valid_data_len+1,:));
            d = [b;c];
            figure;
                imagesc(x,y,d);
                title('Spectrogram');
                xlabel('Time');
                ylabel('Frequency');
                set(gcf,'unit','normalized','position',[0.2,0.2,0.95,0.1]);
        end
        
        function demodePlot(obj, sig, win)
            %DEMODEPLOT plot the demodulation result
            %   Detailed explanation goes here
            figure;
                subplot(3,1,1);
                    obj.spectrum(sig);
                    if nargin < 3 || isempty(win)
                        title('Demodulation');
                    else
                        title(['win',num2str(win)]);
                    end
                subplot(3,1,2);
                    x = obj.despread(sig);
                    [fplot,mgplot] = obj.loraFFT(x);
                    plot(fplot,mgplot);
                    xlim([0,125e3]);
                    ylim([0,max(mgplot)*1.1]);
                subplot(3,1,3);
                    x = obj.despread(sig, 1);
                    [fplot,mgplot] = obj.loraFFT(x);
                    plot(fplot,mgplot);
                    xlim([0,125e3]);
        end
        
        function x = despread(obj, sig, is_down)
            %DESPREAD despread signal of a reception window
            %   Detailed explanation goes here
           if nargin < 3 || isempty(is_down) || is_down == 0
                mChirp = obj.genSymbol(0,1);
           else
                mChirp = obj.genSymbol(0);
           end
            x = sig.*mChirp;
%             disp(['Anlge of First Point:',num2str(angle(sig(1)))]);
        end
        
        function [Data,N] = rfile(obj, filename, samples)
            %RFILE read row signal from files
            %   [Data,N] = rfile(filename) returns the signal array and the number of samples
            %   [Data,N] = rfile(filename,samples) returns the first samples signal points
            fileID = fopen(filename, 'r');
            if fileID == -1, error('Cannot open file: %s', filename); end
            % gr_complex<float> is composed of two float32
            format = 'float';
            row = fread(fileID, Inf, format);
            fclose(fileID);

            if nargin < 3 || isempty(samples)
                N = floor(size(row,1)/2);
            else
                N = samples;
            end
            Data = zeros(1,N);
            Data(1:N) = row(1:2:2*N) + row(2:2:2*N)*1i;
        end
        
        function wfile(obj, filename, data)
            % write binary file
            fileID = fopen(filename, 'wb');
            if fileID == -1, error('Cannot open file: %s', filename); end
            % gr_complex<float> is composed of two float32
            fid=fopen(filename,'wb');
            for idx = 1:size(data,2)
                real_part = real(data(1,idx));
                im_part = imag(data(1,idx));
                fwrite(fid,single(real_part),'single');
                fwrite(fid,single(im_part),'single');
            end
            disp('File Write Finish!');
            fclose(fid);
        end
        
        function B = ampCut(obj, data)
            %AMPCUT extract the useful signal based on the amplitude
            %   Detailed explanation goes here 
            mwin = obj.chirp_n/2;
            A = movmean(abs(data),mwin);
            B = data(A >= max(A)/4);
%             figure;
%                 plot(1:size(A,2),A);
        end
        
        function [snr_db,nPow,sPow] = calSNR(obj, data)
            %CALSNR estimate the SNR of a received LoRa packet
            %   Detailed explanation goes here 
            sig_n = data(obj.chirp_n+(1:obj.chirp_n));
            d = obj.despread(sig_n);
%             snr(real(d),obj.fs);
            [snr_db,nPow] = snr(real(d),obj.fs);
            
            totalNoise = 10^(nPow/10);
            totalSig = 10^(snr_db/10) * totalNoise;
            sPow = 10*log10(totalSig);
            
            fprintf("SNR = %g,nPow = %gdB,sPow = %gdB(%g)\n",snr_db,nPow,sPow,totalSig);
        end
        
        function [cfo,to] = calCFO(obj, upsig, downsig)
            %CALCFO estimate CFO using upsig and downsig
            %   Detailed explanation goes here 
            x = obj.despread(upsig);
            [f,upz] = obj.loraFFT(x);
            [~,I] = max(abs(upz));
            fup = f(I);
            
            y = obj.despread(downsig, 1);
            [~,downz] = obj.loraFFT(y);
            [~,I] = max(abs(downz));
            fdown = f(I);
            
            cfo = (fup + fdown) / 2;
            if abs(cfo) > 50e3
                if cfo < 0
                    cfo = cfo + obj.bw/2;
                else
                    cfo = cfo - obj.bw/2;
                end
            end
            f_to = fdown - cfo;
            to = f_to / obj.bw * obj.chirp_len;
        end
        
        function [sdata,cfo,to] = sync(obj, data)
            %SYNC synchronize the packet
            %   Detailed explanation goes here
            [cfo,to] = obj.calCFO(data(2*obj.chirp_n+(1:obj.chirp_n)), data(11*obj.chirp_n+(1:obj.chirp_n)));
            fprintf('CFO = %g, Time Offset = %g\n',cfo,to);
            offset_n = ceil(to * obj.fs); 
            sdata = data(offset_n:end);
%             figure;
%                 plot(obj.chirp_t,real(data(1:obj.chirp_n)),obj.chirp_t,[zeros(1,offset_n-1),real(data(offset_n:obj.chirp_n))])
        end
        
        function [real_sig,len] = genPacket(obj, codeArray, CFO)
            %GENPAKCKET generate raw signal data
            %   Detailed explanation goes here
            if nargin < 3 || isempty(CFO)
                CFO = 0;
            end
            off_sig = exp(1i*2*pi*CFO*obj.chirp_t);
            upChirp = obj.genSymbol(0).*off_sig;
            downChirp = obj.genSymbol(0,1).*off_sig;
            real_sig = repmat(upChirp,1,8);
            real_sig = [real_sig,obj.genSymbol(2^obj.sf-0).*off_sig,obj.genSymbol(2^obj.sf-0).*off_sig];
            real_sig = [real_sig,downChirp,downChirp,downChirp(1:end/4)];
            for i = codeArray(1:end)
                real_sig = [real_sig,obj.genSymbol(2^obj.sf-i).*off_sig];
            end
            len = size(real_sig,2);
        end
        
        function code_array = loraDecoder(obj, data, outfile)
            %LORADECODER a standard decoder for LoRa packet
            %   Detailed explanation goes here
            if nargin > 2 && ~isempty(outfile)
                fid = fopen(outfile,'w');
                fprintf(fid,'%s\n','win,peak,freq,bin,value');
            end
%             cfo = obj.calCFO(data(2*obj.chirp_n+(1:obj.chirp_n)), data(11*obj.chirp_n+(1:obj.chirp_n)));
            data = [data(1:10*obj.chirp_n),data(round(12.25*obj.chirp_n)+1:end)];
            pkt_len = floor(size(data,2) / obj.chirp_n);
            code_array = zeros(1,pkt_len);
            cfo = 0;
            for base = 0:pkt_len-1
                sig = data(base*(obj.chirp_n)+(1:obj.chirp_n)).*exp(1i*2*pi*(-cfo).*obj.chirp_t);
%                 obj.demodePlot(sig);
                d = obj.despread(sig);
                [freq,z] = obj.loraFFT(d);
                [ma,I] = max(z);
                f = freq(I);
                bin = f/obj.bw*2^obj.sf;
                code_array(base+1) = round(bin);
                fprintf('\n window %d\n',base+1);
                fprintf("peak=%d,  freq=%d,  bin=%.2f[%d]\n",ma,f,bin,mod(round(bin),2^obj.sf));
                
                if nargin > 2 && ~isempty(outfile)
                    fprintf(fid,'%s\n',[num2str(base+1),',',num2str(ma),',',num2str(f),',',num2str(bin),',',num2str(mod(round(bin),2^obj.sf))]);
                end
            end
            if nargin > 2 && ~isempty(outfile)
                fclose(fid);
            end
        end
        
        function dout = decode(obj, data, cfo)
            %LORADECODER a standard decoder for LoRa packet
            %   Detailed explanation goes here
            if nargin < 3 || isempty(cfo)
                cfo = 0;
            end
            symbols_needed = floor(size(data,2) / obj.chirp_n);
            dout = zeros(1,symbols_needed);
            for i = 0:symbols_needed-1
                sig = data(i*(obj.chirp_n)+(1:obj.chirp_n)).*exp(1i*2*pi*(-cfo).*obj.chirp_t);
                [freq,z] = obj.loraFFT(obj.despread(sig));
                [~,I] = max(z);
                dout(i+1) = freq(I)/obj.bw*2^obj.sf;
            end
        end
        
        function coef = corrCoef(obj, sig1, sig2)
            %CORRCOEF calculate the correlation coefficient
            %   (normalization)
            %   Detailed explanation goes here
            dot_product = dot(sig1,sig2);
            magsq1 = abs(sig1).^2;
            magsq2 = abs(sig2).^2;
            energy1 = sum(magsq1);
            energy2 = sum(magsq2);
            coef = abs(dot_product / sqrt(energy1 * energy2));
        end
        
        function [peak_n, heights, freqs, Amps] = peakSearch(obj, data, is_down)
            %PEAKSEARCH search the frequency peaks
            %   Detailed explanation goes here
            heights = [];
            freqs = [];
            Amps = [];
            d = obj.despread(data,is_down);
            [freq,mg] = obj.loraFFT(d);
            threshold = max(max(mg) / 20,1);
            cnt = 1;
            while 1
                % Find a peak.
                [h,I] = max(mg);
                cnt = cnt + 1;
                if h < threshold || cnt > 3
                    break;
                end
                
                % Location of the peak.
                Quarter = round(obj.chirp_n/16);
                %---Left Quarter--------------------
                sigL = [data(1:Quarter),zeros(1,obj.chirp_n - Quarter)];
                d = obj.despread(sigL, is_down);
                [~,mgL] = obj.loraFFT(d);
                %---Right Quarter-------------------
                sigR = [zeros(1,obj.chirp_n - Quarter),data(obj.chirp_n - Quarter+1:end)];
                d = obj.despread(sigR, is_down);
                [~,mgR] = obj.loraFFT(d);                

                % Reconstruct the peak.
                fft_bin = (obj.bw - freq(I))/obj.bw*2^obj.sf;
                if is_down
                    fft_bin = 2^obj.sf - fft_bin;
                end
                e_sig = obj.genSymbol(fft_bin, is_down);
                rc_sig = zeros(1,obj.chirp_n);
                if abs(mgL(I)) > abs(mgR(I))
                    Amp = abs(mgL(I)) / Quarter;               
                    K = min(round(h/Amp),obj.chirp_n);
                    rc_sig(1:K) = Amp * e_sig(1:K);
                else
                    Amp = abs(mgR(I)) / Quarter;
                    K = min(round(h/Amp),obj.chirp_n);
                    rc_sig(obj.chirp_n-K+1:end) = Amp * e_sig(obj.chirp_n-K+1:end);
                end
                d = obj.despread(rc_sig, is_down);
                [~,z] = obj.loraFFT(d);
                mg = mg - z;
                mg(mg<0) = 0;
                
                 % Filter the detected peaks.
                if size(freqs,2) == 0 || min(abs(freqs-freq(I))) > (obj.bw/2^obj.sf)
                    heights = [heights,h];
                    freqs = [freqs,freq(I)];
                    Amps = [Amps,Amp];
                end
            end
            peak_n = size(heights,2);
        end
        
        function [data,len] = mixPkt(obj, pkt1, pkt2)
            len = max(size(pkt1,2),size(pkt2,2));
            if size(pkt1,2) < len
                pkt1 = [pkt1,zeros(1,len-size(pkt1,2))];
            else
                pkt2 = [pkt2,zeros(1,len-size(pkt2,2))];
            end
            data = pkt1 + pkt2;
        end
        
        function dataout = anoise(obj,datain,snr)
        %ANOISE Summary of this function goes here
        %   Detailed explanation goes here
            amp_sig = mean(abs(datain));
            org_snr = obj.calSNR(datain);
            amp_noise = amp_sig / 10^(org_snr/20);

            amp_awgn = amp_sig/10^(snr/20) - amp_noise;
            if amp_awgn < 0
                amp_awgn = 0;
            end
            len_data = length(datain);
            dataout  = datain + (amp_awgn/2 * randn([1 len_data]) + 1i*amp_awgn/2 * randn([1 len_data]));
        end

    end
end