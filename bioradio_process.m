%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Akshay Paul
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for processing BioRadio .csv output
clear all; close all; clc;
%% load in data from BR output file
filename1 = ['BioRadio Day6'...
    '.csv'];
%
channels = [1:1];               % **set which BR channels being used**
numchans = length(channels);

d_in = readmatrix(filename1);   % BR .csv format is 1 col of time + 8 cols of data + 1 col of events 
data = d_in(:,channels+1);      % isolate the EEG channels only
events = d_in(:,3);           % isolate event marker column    
events_idx = find(events==1);   % just the indices for positive events
%clear d_in;
% 
t_in = readmatrix(filename1, 'Range', 'A:A','OutputType', 'duration');
t_in = t_in(2:end,1);
Fs_calc = 1/seconds(t_in(2)-t_in(1));
fprintf('Is the expected Fs = %d-Hz?\n', Fs_calc)
disp('The first 4 time points look like...')
disp(t_in(1:4))

%% Find NaN in EEG data
[row, col] = find(isnan(data(:,1)));  % **find the drops in the first channel only; should match
% row, col side by sidefind
%drops = [row col];  % the drops seem to occur at the same time across channels, need to check
drops = row;
dlen = length(drops);
fprintf('There were %d data drops detected.\n',dlen);

figure(1);
plot(t_in,'.'); 
hold on
for i = 1:dlen
    xline(drops(i,1));
end
yyaxis right; 
plot(data(:,1));
hold off
title(['timestamps, datadrops (',num2str(dlen),') and a CH of EEG']);

%% datadrop smoothing (only if there were drops)

%set missing value to last known good value
trig = 0;
lstgood_cur =  drops(1)-1;
lstgood_idx = zeros(dlen,1);
lstgood_idx(1) = lstgood_cur;

for i = 2:dlen
    if (drops(i) - drops((i-1))) == 1
        lstgood_idx(i) = lstgood_cur; 
    else
        lstgood_cur = drops(i)-1;
        lstgood_idx(i) = lstgood_cur;
    end
end

%% create a fixed dataset (optional, only if there were drops)
data_fix = data; %** adjust **
for i = 1:dlen
    data_fix(drops(i),:) = data(lstgood_idx(i),:);
end

figure(1);
hold on;
plot(data_fix(:,3))
title('datadrops fixed with last good value')
hold off;

%% select which data to proceed with
eegdata = data';    %could be data_fixed'; ** select **
%eegdata = data_fix';  
%% eeglab section, open eeglab

Fs = 1000;
EEG = pop_importdata('data', eegdata, 'nbchan', numchans, 'srate', Fs);

%     Bandpass filter the data from 2 to 100 Hz
%EEG = pop_eegfiltnew(EEG, 2, 58, 1650, 0, [], 0);

%EEG = pop_eegfiltnew(EEG, 58, 62, 1650, 1, [], 0);
%EEG.setname= strcat(EEG.setname,'_2to100BP_removed2s');
%EEG.filename = strcat(EEG.setname,'.set');
%EEG.datfile = strcat(EEG.setname,'.fdt');    
%pop_saveset(EEG,'filename',EEG.filename,'filepath',pwd);

%% PSD Pop_Spectopo EEGLAB - skip if using preferred
win_size = 10;  % in seconds  ** adjust **
freqs_plt = [1 70];
time_plt(1) = events_idx(2) * 1000/EEG.srate;  % convert event idx to ms (start)   
time_plt(2) = events_idx(3) * 1000/EEG.srate;  % convert event idx to ms (end)

figure(201);
pop_spectopo(EEG, 1, ...   % 1 is dataflag to process EEG in (0: components)   
    time_plt, ...
    'EEG' , ...
    'freqrange',freqs_plt,...
    'electrodes','off',...
    'winsize',1000*win_size);
title('test title')

%% PSD Spectopo EEGLAB - Preferred 
str_mark = 2;
end_mark = 3;
toplot = 'ON';       %** Turn Plot ON or OFF

win_size = 10;  % in seconds  ** adjust **
freqs_plt = [1 80];
time_plt(1) = events_idx(str_mark)+5000 * 1000/EEG.srate;  % convert event idx to ms (start)
%time_plt(1) = 5000 * 1000/EEG.srate;  % convert event idx to ms (start)   
time_plt(2) = events_idx(end_mark)-5000 * 1000/EEG.srate;  % convert event idx to ms (end)

if time_plt(1)/1000~=EEG.xmin | timerange(2)/1000~=EEG.xmax
   %posi = 100000;
   posi = round( (time_plt(1)/1000-EEG.xmin)*EEG.srate )+1;
   posf = round( (time_plt(2)/1000-EEG.xmin)*EEG.srate )+1;
   pointrange = posi:posf;
end
if exist('pointrange') == 1, SIGTMP = EEG.data(:,pointrange,:); totsiz = length(pointrange);
else
    SIGTMP = EEG.data; totsiz = EEG.pnts;
end

SIGTMP = reshape(SIGTMP, size(SIGTMP,1), size(SIGTMP,2)*size(SIGTMP,3));

f = figure(30);
[spectra, freqs] = spectopo(SIGTMP, totsiz, EEG.srate, ...
    'freqrange', freqs_plt, ...
    'winsize', 1000*win_size, ...
    'limits', [0 60 NaN NaN], ... %xmin xmax ymin ymax
    'plot', toplot);                                    
%% calculate noise Vrms
noise_br_only = rms(spectra(1:701))

noise_vrms = sqrt(sum((10.^(spectra/10))*freqs)/length(spectra))

%% plot additional spectra
f = figure(30); 
hold on
plot(freqs, spectra, 'LineWidth', 2);
axis([0 70 -190 -130])
hold off
box on
set(gca,'FontSize',14);
set(gca,'linewidth',2);
xlabel('Frequency (Hz)','FontWeight', 'bold','FontSize',19);
ylabel('PSD (10log10(V^{2}/Hz))','FontWeight', 'bold','FontSize',19);
title("Day6 EyesOpened");

%% add plot information 
figure(30)
title('PodStamp EEG channels ASSR 40 Hz')
set(gca,'XGrid', 'on')
dim = [0.55 0.8 0.1 0.1];
%str = {'3Ch 3M Ext. RMst Ref'};
%annotation('textbox',dim,'String',str,'FitBoxToText','on');
legend('EEG #1', ...
        'EEG #2',...
        'EEG #3'),...
set(gcf,'position',[0,0,1440,720]);
%% FFT 

%Fs = 83333; 
gain = 1;
%data = data(26500:271800);
N = length(SIGTMP);

dFT = fft(SIGTMP(1,:));    % fft
dFT = dFT(1:N/2+1); % DC to Nyquist 
dFT = dFT/gain;
dPD = (1/(Fs*N)) * abs(dFT).^2;  % periodogram
dPD(2:end-1) = 2*dPD(2:end-1);   % factor of 2 to conserve power
f = 0:Fs/N:Fs/2;
%%
figure(1);
hold on;
plot(f, 10*log10(dPD));
grid on
title('Periodogram  of Bioradio Recording'); 
xlabel('Frequency (Hz)');
ylabel('PSD (10log10(V^{2}/Hz))');
set(gca, 'XScale', 'log')

figure(2);
hold on;
plot(f(119:end), abs(dFT(119:end)));
grid on
title('FFT of Bioradio Recording'); 
xlabel('Frequency (Hz)');
%ylabel('V');
xlim([0 70])
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')

figure(3);
hold on;
plot(f(:), angle(dFT(:))*180/pi);
grid on
title('Phase of Bioradio Recording'); 
xlabel('Frequency (Hz)');
%ylabel('V');
xlim([0 70])
%set(gca, 'YScale', 'log')
%set(gca, 'XScale', 'log')
%% Rectified signal FFT

%% load the .csv sound file and run FFT rectified
y = data(:,numchans);

%%%%%%%%%%%%%% Min: 60 seconds mark %%%%%%%%%%%%%%
str_mark = 1;
str_idx = events_idx(str_mark);
end_idx = str_idx + Fs*60;  % Fs*60: number of samples for 60seconds

y = y(str_idx:end_idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


offset = 1;                                         % Offset setting
T = 1/Fs;                                           % Sampling period
l = length(y);                                      % Length of recorded signal
L = length(y(offset:l));                            % Length of signal with offset
t = (0:L-1)*T;                                      % Time vector
X = y(offset:L);

f = Fs*(0:(L/2))/L;
%Y = fft(hann(L).*X');
Y = fft(abs(X));
myphase = angle(Y/L);
myphase = myphase(1:L/2+1);
myphase(2:end-1) = 2*myphase(2:end-1);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

figure(8);
subplot(2,1,1); 
plot(f,P1,'Color',[0.75 0.75 0.75]); 
title('.WAV load FFT Mag and Phase BM 40Hz-AM WN Rectified')
ylabel('|P1(f)|')
axis([0 500 -2e-6 5e-5])
subplot(2,1,2);
plot(f, myphase);
ylabel('Phase')
xlabel('f (Hz)');


%% save figures
f = figure(22);
save_name = 'YX_R_Airpod_BR_Alpha_floating';
saveas(gca, save_name,'epsc')
saveas(gca, save_name,'svg')
saveas(gca, save_name,'jpg')
savefig([save_name,'.fig'])

%% Spectrogram Test
spectrogram(eegdata(1,:),'yaxis ');
