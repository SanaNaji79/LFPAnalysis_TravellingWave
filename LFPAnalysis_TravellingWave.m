%% README
%Investigation of Cortical Traveling Waves in Array dataset
%array includes 49 single electrodes distributed in an area of 12𝑚𝑚2. Distance between neighboring electrodes are 400 𝜇𝑚.
%the first partis the LFP analysis(filtering the pink noise finding the dominant frequency of each electrode and draw the travelling wave)
%then having a butterworth filter and investigate the frequency and phase
%content of the data
%%
clear ;
clc ;
%% loading data
load 'C:\Users\Sana\OneDrive\Desktop\semester8\AdvanceNeurosience_ghazizadeh\homework\problems\ArrayData\ArrayData.mat' ;
load 'C:\Users\Sana\OneDrive\Desktop\semester8\AdvanceNeurosience_ghazizadeh\homework\problems\ArrayData\CleanTrials.mat' ;
%% LFP analysis part a
%% initialization
t = Time ;
freq = 200 ;
len = 641 ;
f = freq*(0:(len/2))/len ;
%% plot the fourier transform for a random trial of a random neuron
neuron = [3 , 9 , 15 , 20 , 24 , 31 , 42 , 48] ;
trial = 819 ; 
for i = neuron
    input = chan(i).lfp ;
    subplot(4 , 2 , find(neuron == i)) ;
    plot(Time , input(: , trial)) ;
    xlim([-1.2 , 2]) ;
    xlabel('time (s)') ;
    ylabel('voltage') ; 
    title(['Raw data unit' num2str( i )]) ;
end
input = chan(9).lfp ;
p = 1 ;
b = psd(input(: , trial) , f , len , p) ;
output = reject_pink_noise(b , f , p) ;
output = reject_EC_frequency(output , f , p) ;
[output , f1] = reject_dc(output , f , p) ;
%% plotting the spectrum for a random unit
p = 0 ; 
clean = Intersect_Clean_Trials ;
filter = 1 ;
f = freq*(0:(len/2))/len ;
input = chan(9).lfp ;
[output1 , output2 , f] = dominant(input ,clean , filter , f , len , p) ;
%%
figure ;
hold on ;
title('Final Signal Single-Sided Amplitude Spectrum') ;
xlabel('f (Hz)') ;
ylabel('power') ;
[x , y] = size(output2) ;
for i = 1 : y
    plot(f , output2(: , i)) ;
end
plot(f , mean(output2 , 2) , 'LineWidth' , 2 , 'Color'  ,[0 , 0 , 0]) ;
xlim([f(1) , f(end)]) ;
%% LFP analysis
p = 0 ;
[output] = reject_noise(chan(1).lfp , f , len , p) ;
%%
FrequencyLimit = [1 80] ;
ylimPlots = [2 50];
LFP_s = 0 ;
for i = 1:48
    [LFP_s1 , LFP_f] = LFP_analysis(chan(i).lfp(: , Intersect_Clean_Trials) , freq , FrequencyLimit) ;
    LFP_s = LFP_s + LFP_s1 ;
end
LFP_s = LFP_s/48 ;
%%
surface(t , LFP_f(1:55) , abs(LFP_s(1:55 , :)))
axis tight
shading interp
xlabel("Time [ms]")
ylabel("Frequency [Hz]")
ylim(ylimPlots)
%clim([0 1])
colorbar
colormap turbo
%%
output_psd = zeros(48 , 321) ;
for i = 1:48
    output_psd(i , :) = abs(psd_average(chan(i).lfp(: , Intersect_Clean_Trials)  , f , len , p)) ;
end
%%
figure
surface(1:48 , f(17:146) , output_psd(: , 17:146)') ;
axis tight
shading flat
colormap turbo
xlabel('Channel') ;
ylabel('frequency (Hz)') ;
%% finding the dominant frequency
freq_dominant = zeros(3 , 48) ;
f = freq*(0:(len/2))/len ;
p = 0 ; 
clean = Intersect_Clean_Trials ;
filter = 1 ;
y = [] ;
n= 5 ;
for i = 1:48
    input = chan(i).lfp ;
    [freq_dominant(2 , i) , out2 , f1] = dominant(input ,clean , filter , f , len , p) ;
    out2 = mean(out2 , 2) ;
    [a , b] = max(out2) ;
    freq_dominant(3 , i) = f1(b) ;
    subplot(6 , 8 , i) ;
    plot(f1 , out2) ;
    title(['Average PSD for unit' num2str( i )]) ;
end
%%
for i = 1:48
    input = chan(i).lfp ;
    [x , freq_dominant(1 , i) , f1] = dominant_freq(input ,clean , filter , f , len , p) ;
    subplot(6 , 8 , i) ;
    plot(f1 , x) ;
    title(['Average PSD for unit' num2str( i )]) ;
end
figure ;
plot(freq_dominant(1 , :)) ;
hold on ;
plot(freq_dominant(2 , :)) ;
plot(freq_dominant(3 , :)) ;
legend('average of time series' , 'average of dominant frequencies' , 'average of frequency domians') ;
%% LFP analysis part b
p = 1 ;
[output1] = num2matrix(freq_dominant(1 , :) , ChannelPosition , p , 'topography of dominant frequency') ;
[output2] = num2matrix(freq_dominant(2 , :) , ChannelPosition , p , 'topography of dominant frequency') ;
[output3] = num2matrix(freq_dominant(3 , :) , ChannelPosition , p , 'topography of dominant frequency') ;
%% LFP analysis part c
%% the contour of frequency channels
f = freq*(0:(len/2))/len ;
p = 0 ; 
clean = Intersect_Clean_Trials ;
filter = 1 ;
output_new = [] ;
bandpass = 1 ;
for i = 1:48
    input = chan(i).lfp ;
    [out2 , f1] = average_frequency(input ,clean , filter , f , len , p , bandpass) ;
    output_new = [output_new ; mean(out2)] ;
end
%% plotting
colormap jet ;
imagesc((1:48) , f1+5 , output_new') ;
xlabel('channels') ;
ylabel('frequency') ;
%ylim([]) 
colorbar ;
%%
colormap jet ;
% input_new = ifft(output_new(1 , :)) ;
% plot(input_new)
% figure
stft(output_new(1 , :) , 1/200) ;
%%
x = chan(2).lfp(: , 1) ;
fs = 200 ;
y = bandpass(x , [10 45] , fs) ;
subplot(3 , 1 , 1) ;
plot(t , x) ;
xlabel('time (s)') ;
ylabel('voltage (mV)') ;
title('Raw Data') ;
xlim([-1.2 , 2]) ;
subplot(3 , 1 , 2) ;
plot(t , y) ;
xlabel('time (s)') ;
ylabel('voltage (mV)') ;
title('Bandpass Filtered Data') ;
xlim([-1.2 , 2]) ;
z = lowpass(x , 10 , fs) ;
subplot(3 , 1 , 3) ;
plot(t , z) ;
xlabel('time (s)') ;
ylabel('voltage (mV)') ;
title('Lowpass Filtered Data') ;
xlim([-1.2 , 2]) ;
%% find the bandpass filtered signal
newchan = zeros(48 , 641 , 823) ;
for i =1:48
    for j = 1:823
        newchan(i , : , j) = bandpass(chan(i).lfp(: , j) , [10 45] , fs) ;
    end
end
%%
colormap jet ;
fs = 200 ;
a = 0 ;
for j = 1:48
for i = 1:Intersect_Clean_Trials
    %stft(y , 200 , 'Window',kaiser(50,2)) ;
    %x = bandpass(chan(i).lfp(: , j) , [10 45] , fs) ;
    x1 = chan(j).lfp(: , i) ;
    x = bandpass(x1 , [10 45] , fs) ;
    win = hamming(18,'periodic');
    s = stft(x1,fs,'Window',win,'OverlapLength',16,'FFTLength',350);
    a = s + a ;
end
end
%%
colormap jet ;
a = a./(48*490) ;
y_l = f(1:312) + 5;
x_l = (1:350/2) ;
imagesc(abs(a(1:350/2 , :))) ;
ylim([1 350/2]) ;
xlim([1 312]) ;
xlabel('time') ;
ylabel('frequency') ;
% figure ;
% colormap jet ;
% imagesc(abs(a)) ;
% ylim([0 100]) ;
% xlim([0 62]) ;
%% LFP analysis part d
%% Phase propagation(Traveling waves) part a
[b,a] = butter(1,[10 45]/100);
freqz(b , a) ;
x = chan(1).lfp(: , 1) ;
y = filter(b,a,x);
figure ;
plot(t,x)
hold on
plot(t,y) ; 
xlim([-1.2 , 2]) ;
xlabel('time(s)') ;
ylabel('amplitude') ;
legend('Input Data','Filtered Data') ; 
p = 1 ;
output = psd(x , f , len , p) ;

output = psd(y , f , len , p) ;
%%
p = 0 ;
[b,a] = butter(1,[10 45]/100);
filtered = zeros(48 , 490 , 641) ;
filtered_psd = zeros(48 , 490 , 321) ;
for i = 1:48
    for r = 1:490
        j = Intersect_Clean_Trials(r) ;
        g = filter(b,a,chan(i).lfp(: , j)) ;
        filtered(i , r , :) = reshape(g , [1 , 1 , 641]) ;
        filtered_psd(i , r , :) = psd(filtered(i , r , :) , f , len , p) ;
    end
end
%%
dom_av = [] ;
dom_ide = [] ;
dom_time = [] ;
for i = 1:48 
    av = mean(filtered(i , : , :)) ;
    av_psd = psd(av , f , len , p) ;
    [m1 , m2] = max(av_psd) ;
    dom_time = [dom_time , f(m2)] ;
end
for i = 1:48 
    av = mean(filtered_psd(i , : , :)) ;
    [m1 , m2] = max(av) ;
    dom_av = [dom_av , f(m2)] ;
end
for i = 1:48
    s = [] ;
    for j = 1:490 
        [m1 , m2] = max(filtered_psd(i , j , :)) ; 
        s = [s , f(m2)] ;
    end
    dom_ide = [dom_ide , mean(s)] ;
end
plot(dom_av) ;
hold on 
plot(dom_ide) ;
plot(dom_time) ;
legend('average over frequency domain' , 'average over time domain' , 'average over dominant frequency') ;
xlim([1 , 48]) ;
xlabel('unit number') ;
ylabel('frequency (Hz)') ;
title('dominant frequency') ;
figure ;
colormap jet ;
pspectrum(y,fs,'spectrogram')
%%
p = 1 ;
[output1] = num2matrix(dom_av, ChannelPosition , p , 'topography of dominant frequency') ;
title('topography of dominant frequency (average over frequency domain)') ; 
[output2] = num2matrix(dom_time , ChannelPosition , p , 'topography of dominant frequency') ;
title('topography of dominant frequency (average over time domain)') ; 
[output3] = num2matrix(dom_ide , ChannelPosition , p , 'topography of dominant frequency') ;
title('topography of dominant frequency (average over dominant frequency)') ;
%% Phase propagation(Traveling waves) part b

Hilbert_filtered1 = zeros(48 , 490 , 641) ;
for i = 1:48
    for j = 1:490
        y = filtered(i , j , :) ;
        z = hilbert(y);
        int = y + 1i*z ;
        w = unwrap(angle(z)) ;
        Hilbert_filtered1(i , j , :) = w ;
    end
end
%% Phase propagation(Traveling waves) part c
Hilbert_filtered2 = mean(Hilbert_filtered1) ;
p = 1 ;
Hilbert_filtered = num2matrix2(Hilbert_filtered1 , ChannelPosition) ;
Hilbert_filtered_new = num2matrix3(Hilbert_filtered1 , ChannelPosition) ;
Hilbert_filtered_new_mean = mean(Hilbert_filtered_new , 3) ;
%%
Hilbert_filtered3 = mean(Hilbert_filtered , 3) ;
h = figure ;
folder = 'C:\Users\Sana\OneDrive\Desktop\semester8\AdvanceNeurosience_ghazizadeh\homework\assignments\HW4_SanaAminnaji_98104722' ;
for i = 1:641
    ttt = reshape(Hilbert_filtered3(: , : , 1 , i) , [5 , 10]) ;
    colormap gray ;
    imagesc(ttt) ;
    shading interp ;
    pause(0.01) ;
    make_animation( h,i,'ggg.gif',folder ) ;
end
%% Phase propagation(Traveling waves) part d
grad = zeros(5 , 10 , 641) ;
grad_mean = zeros(5 , 10 , 641) ;
for i = 1:641
    avv1 = 0 ;
    rt = reshape(Hilbert_filtered_new_mean(: , : , 1 , i) , [5 , 10]);
    [r1 , r2] = gradient(rt) ;
    grad(: , : , i) = sqrt(r1.*r1 + r2.*r2) ;
    for j = 1:480
        r = reshape(Hilbert_filtered_new(: , : , j , i) , [5 , 10]);
        [g11 , g22] = gradient(r) ;
        avv = sqrt(g11.*g11 + g22.*g22) ;
        avv1 = avv1 + avv ;
    end
    grad_mean(: , : , i) = avv1/480 ;
end
pgd = grad./grad_mean ;
%%
y = 1:5 ;
x = 1:10 ;
h = figure ;
for i = 1:641
    quiver(x,y,reshape(grad(1 , : , : , i), [5,10]),reshape(grad(2 , : , : , i),[5, 10]))
    pause(0.01) ;
    make_animation( h,i,'gg.gif',folder ) ;
    
end
%% Phase propagation(Traveling waves) part e

%% Phase propagation(Traveling waves) part f

%% Phase propagation(Traveling waves) part g

%% function declaration
function [output] = psd(input , f , len , p) 
y = fft(input) ;
P2 = abs(y/len) ;
output = P2(1:len/2+1) ;
output(2:end-1) = 2*output(2:end-1) ;
if p
    figure ;
    plot(f,output) ; 
    title('Single-Sided Amplitude Spectrum of raw data') ;
    xlabel('f (Hz)') ;
    ylabel('power') ;
end
end
function [output] = reject_pink_noise(input , f , p)
m = polyfit(f' , log(input) , 1) ;
out1 = log(input)' - (m(1).*f + m(2)) ;
output = exp(log(input)' - (m(1).*f + m(2))) ;
if p
    figure ;
    plot(f , log(input)) ;
    hold on ;
    plot(f , m(1).*f + m(2)) ;
    plot(f , out1) ;
    legend('log(Power)' , 'Pink noise' , 'log(Normalized power)') ;
    xlabel('f (Hz)') ;
    ylabel('log(power)') ;
    title('Rejecting pink noise') ;
    figure ;
    plot(f , input) ;
    hold on ;
    plot(f , output) ;
    legend('Power' , 'Normalized power') ;
    xlabel('f (Hz)') ;
    ylabel('Power') ;
end
end
function [output] = reject_EC_frequency(input , f , p)
loc = find(f>57 & f<63) ;
input(loc) = 0 ;
output = input ;
if p
    figure ;
    plot(f,output) ; 
    title('Single-Sided Amplitude Spectrum of proceesed data') ;
    xlabel('f (Hz)') ;
    ylabel('power') ;
end
end
function [output , f] = reject_dc(input , f , p) 
loc = find(f>2);
f = f(loc(1):end) ;
output = input(loc(1):end) ;
if p
    figure ;
    plot(f,output) ; 
    title('Final Signal Single-Sided Amplitude Spectrum') ;
    xlabel('f (Hz)') ;
    ylabel('power') ;
    xlim([f(1) , f(end)]) ;
end
end
function [output] = average(input , clean , filter , p)
if filter
    input = input(: , clean) ;
end
input = normalize(input) ; %%%%
%input = input - mean(input) ;
output = mean(input , 2) ;
output = normalize(output) ; %%%%
if p 
    figure ;
    plot(output) ;
end
end
function [output , f] = reject_dc1(input , f , dc) 
loc = find(f>dc);
f = f(loc(1):end) ;
output = input(loc(1):end) ;
end
function [output1 , output3 , f1] = dominant_freq(input ,clean , filter , f , len , p)
input1 = average(input , clean , filter , p) ;
input2 = psd(input1 , f , len , p) ;
output1 = reject_pink_noise(input2 , f , p) ;
output1 = reject_EC_frequency(output1 , f , p) ;
[output1 , f1] = reject_dc1(output1 , f , 5.5) ;
%output2 = output1(34:end) ;
[a , output3] = max(output1) ;
output3 = f1(output3) ;
end
function [output , out , f] = dominant(input ,clean , filter , f , len , p)
if filter
    input = input(: , clean) ;
end
[a , b] = size(input) ;
out = zeros(314 , b) ;
for i = 1:b
    x = psd(input(: , i) , f , len , p) ;
    x = reject_pink_noise(x , f , p) ;
    x = reject_EC_frequency(x , f , p) ;
    [out(: , i) , f1] = reject_dc(x , f , p) ;
end
f = f1 ;
[a , out1] = max(out) ;
out2 = out1 ;
for i = 1:length(out1)
    out2(i) = f(out1(i)) ;
end
output = mean(out2) ;
end
function [output] = num2matrix(input , matrix , p , ti)
output = matrix ;
[a , b] = size(matrix) ;
for i = 1:a
    for j = 1:b
        if isnan(matrix(i , j)) 
            output(i , j) = 6 ; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        else
            output(i , j) = input(matrix(i , j)) ;
        end
    end
end
if p
    figure ;
    title('ti') ;
    imagesc(output) ;
    shading interp
    colorbar ;
end
end
function [out , f1] = average_frequency(input ,clean , filter , f , len , p , bandpass)
if filter
    input = input(: , clean) ;
end
[a , b] = size(input) ;
%out = zeros(314 , b) ;
out = [] ;
for i = 1:b
    x = psd(input(: , i) , f , len , p) ;
    x = reject_pink_noise(x , f , p) ;
    x = reject_EC_frequency(x , f , p) ;
    [x , f2] = reject_dc(x , f , p) ;
    if bandpass
        loc = find(f2>45) ;
        x = x(1:loc(1)) ;
        f1 = f2(1:loc(1)) ;
    end
    out = [out ; x] ;
end
end

function [output] = num2matrix2(input , matrix)
output = matrix ;
[a , b] = size(matrix) ;
output = zeros(a , b , 480 , 641) ;
for t = 1:641
    for tr = 1:420
        for i = 1:a
            for j = 1:b
                if isnan(matrix(i , j))
                    output(i , j , tr , t) = -2 ;
                else
                    output(i , j , tr , t) = cos(input(matrix(i , j) , tr , t)) ;
                end
            end
        end
    end
end
end
function make_animation( h,index,filename,folder )
drawnow
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
if index == 1
    imwrite(imind,cm,fullfile(folder,filename),'gif', 'Loopcount',inf);
else
    imwrite(imind,cm,fullfile(folder,filename),'gif','WriteMode','append');
end
end
function [output] = num2matrix3(input , matrix)
output = matrix ;
[a , b] = size(matrix) ;
output = zeros(a , b , 480 , 641) ;
for t = 1:641
    for tr = 1:420
        for i = 1:a
            for j = 1:b
                if isnan(matrix(i , j))
                    output(i , j , tr , t) = -2 ;
                else
                    output(i , j , tr , t) = input(matrix(i , j) , tr , t) ;
                end
            end
        end
    end
end
end
function [output] = reject_noise(input , f , len , p)
[a b] = size(input) ;
output = zeros(a , b) ;
for i = 1:b
    x = psd(input(: , i) , f , len , p) ;
    x = reject_pink_noise(x , f , p) ;
    x = reject_EC_frequency(x , f , p) ;
    output(: , i) = ifft(x , a) ; %%%%%%%%%%%%%%%%%%%%%%
end
end
function [LFP_s , LFP_f] = LFP_analysis(input , freq , FrequencyLimit)
[a b] = size(input) ;
[LFP_s , LFP_f] = cwt(input(: , 1), freq, FrequencyLimits=FrequencyLimit);
for i = 2:b
    %[b1, a1] = iirnotch(w0, bw);
    %LFPSig = filtfilt(b1, a1, input(: , 2));
    [LFP_s1, LFP_f] = cwt(input(: , i), freq, FrequencyLimits = FrequencyLimit);
    LFP_s = LFP_s + LFP_s1 ;
end
LFP_s = LFP_s/b ;
end
function [output] = psd_average(input , f , len , p)
[a b] = size(input) ;
output = 0 ;
for i = 1:b
    x = psd(input(: , i) , f , len , p) ;
    x = reject_pink_noise(x , f , p) ;
    %x = reject_EC_frequency(x , f , p) ;
    output = output + x ;
end
output = output/b ;
end