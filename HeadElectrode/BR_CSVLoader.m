%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: Min Suk Lee
% Date: 8/16/2021
% Description: Load and format BioRadio csv file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
%% load in data from BR output file
% filenames %
filename1 = "Bleach_Electrode_ASSR(5).csv";
filename2 = "Electro_Electrode_ASSR(5).csv";
filename3 = "Ink_Electrode_ASSR(4).csv";
filename4 = "Wet_Electrode_ASSR.csv";

%FILES = [filename1, filename2, filename3, filename4];
FILES = [filename1,filename2]; % Compile filenames %

LEGEND = FILES; % Create Filenames as Legends for graphing

% Extract data from file
for f = 1:length(FILES)
    file = readtable(FILES(f)); % get file as table
    if sum(sum(ismissing(file))) > 1 % check for any package drops
        disp("Package Drop: " + sum(sum(ismissing(file))));
        break;
    end
    NUMCHAN(f) = {width(file)-2}; % get number of channels
    EVENT_IDX(f) = {find(file.(width(file))==1)}; % get event index
    DATA(f) = {file.(2)}; % get data
    TIMESTAMP(f) = {file.(1)}; % get time stamp
end
%% Initialize Variables %%
Fs = 1000; % Sampling Rate %

%% EEGlab Section
for f = 1:length(FILES)
    EEG(f) = pop_importdata('data', DATA{f}', 'nbchan', NUMCHAN{f}, 'srate', Fs);
end
%% Trim Data and PSD
% Marker Index for each file %
str_mark = [1,1];
%end_mark = [2,2];
str_offset = 5*Fs; % offset from the start index
duration = 60*Fs-1; % duration in seconds
win_size = 2;
freqs_plt = [1 70];
toplot = 'OFF';


for f = 1:length(FILES)
    if length(str_mark)~=length(FILES)
        disp("str_mark length does not match the number of files")
        break;
    end
    totsiz = EEG(f).pnts;
    str_idx = EVENT_IDX{f}(str_mark(f));
    sigtmp = DATA{f}(str_idx:(str_idx+duration))';
    [power, freqs] = spectopo(sigtmp, 0, EEG(f).srate, ...
        'freqrange', freqs_plt, ...
        'winsize', Fs*win_size, ...
        'limits', [0 60 NaN NaN], ... %xmin xmax ymin ymax
        'plot', toplot);
    POWER(f) = {power};
    FREQ(f) = {freqs};
end

%% Plot Figure
close all;
figure(20);
hold on;
for f = 1:length(FILES)
    plot(FREQ{f},POWER{f},'LineWidth',2)
end
axis([0 70 -190 -130])
box on
set(gca,'FontSize',14);
set(gca,'LineWidth',2);
xlabel('Frequency [Hz]','FontWeight', 'bold','FontSize',15);
ylabel('Log Power Spectral Density [10*log_{10}(\muV^{2}/Hz)]','FontWeight', 'bold','FontSize',12);
lgd = legend(LEGEND,'Interpreter','none');
lgd.FontSize = 15;
lgd.FontWeight = 'bold';
lgd.Location = 'northwest';
ylim([-140 -60]);
xlim([1 70]);