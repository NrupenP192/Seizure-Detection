clc
clear
close all

%Loading Day Data
data = readmatrix('Day1.txt');
fs = 200; %sampling frequency in hz after downsampling
t = (0:1/fs:(length(data)-1)/fs)'; %time vector

CortexIpsilateral = data(:,2);
CortexContralateral = data(:,1);
Hippocampus = data(:,3);

%%Loading for test data, only uses the first 1500 seconds, leave suppressed
%%for actual analysis
CortexIpsilateral = CortexIpsilateral(1:600000);
t = t(1:600000);
CortexIpsi = CortexIpsilateral; %making a copy of the raw data


CortexIpsilateral = wdenoise(CortexIpsilateral,14, ...
    'Wavelet', 'db10', ...
    'DenoisingMethod', 'BlockJS', ...
    'ThresholdRule', 'James-Stein', ...
    'NoiseEstimate', 'LevelDependent'); %wavelet denoising


%% Calculating Teager
% Define the window size and overlap for moving windows
window_size = 300; % 200hz = 1s window
window_overlap = round(window_size/1.1111); % 90% overlap

% Compute the windows
windows = buffer(CortexIpsilateral, window_size, window_overlap, 'nodelay');
windowsT = buffer(t, window_size, window_overlap, 'nodelay');

% Compute the Teager non-linear energy values
%Teager Modifiers
Tn = 3;
%Teager discrete formula where n can be changed, and is vectorized: y = x(n+1:end-n).^2 - x(1:end-2n) . x(2*n+1:end); %general formula
teager_values = (sum(windows(Tn+1:end-Tn, :).^2 - windows(1:end-2*Tn,:).*windows(2*Tn+1:end,:),1));
teager_avg = teager_values./window_size;


%teager_values = (sum(windows(2:end-1, :).^2 - windows(1:end-2, :).*windows(3:end, :), 1));%x^2 - x(n-1)*x(n+1), %time step = 0.1s, %f = 10hz
time_axis = windowsT(100,:); %taking the midpoint of each window as the index for the teager value
teager_avg = smoothdata(teager_avg,'movmean','SmoothingFactor',0.1,'SamplePoints',time_axis);
%% Calculating the Step Function Based on the Teager values
minteager = mink(abs(teager_avg), 50);
avg_min = mean(minteager);
SLEThreshold = avg_min*12;
NoiseThreshold = avg_min*100;
%Threshold Plot
SLEThresholdPlot = zeros(1,length(teager_avg));
SLEThresholdPlot = SLEThresholdPlot + SLEThreshold;

indexSLE = find(teager_avg >= SLEThreshold & teager_avg <= NoiseThreshold);
StepEnvelope = zeros(1, length(teager_avg));

StepFunctionAmplitude = 5;
StepEnvelope(indexSLE) = StepFunctionAmplitude;
SLEStartCode = [0, StepFunctionAmplitude];
SLEEndCode = [StepFunctionAmplitude, 0];

%Removing Low duration events (<0.2s)
%This step should eliminate momentary increases to teager because of motion
%artifacts or events in-between SLEs that arent SLEs but have high enough
%frequencies at a low duration to be picked up
[SLEStartTimeNew, SLEEndTimeNew, SLEStartIndexNew, SLEEndIndexNew] = findStartEndPoints(StepEnvelope, time_axis, SLEStartCode, SLEEndCode);

%Removing events that are under 0.5 seconds long
StepEnvelopeReduced = StepEnvelope;
for a = 1:length(SLEEndIndexNew)
    if SLEEndTimeNew(a)-SLEStartTimeNew(a) < 0.5
        StepEnvelopeReduced( SLEStartIndexNew(a): SLEEndIndexNew(a)) = 0;
    end
end

[SLEStartTime, SLEEndTime, SLEStartIndex, SLEEndIndex] = findStartEndPoints(StepEnvelopeReduced, time_axis, SLEStartCode, SLEEndCode);

%Combining events that are 2.5 seconds or less apart from each other
StepEnvelopeCombined = StepEnvelopeReduced;
for b = 1:length(SLEEndIndex)-1
    if SLEStartTime(b+1) - SLEEndTime(b) <= 1
        StepEnvelopeCombined( SLEEndIndex(b): SLEStartIndex(b+1) ) = StepFunctionAmplitude;
    end
end

%Check for overall event duration after combining events, to identify
%events that are over 5s in duration
[SLEStartTimeNew, SLEEndTimeNew, SLEStartIndexNew, SLEEndIndexNew] = findStartEndPoints(StepEnvelopeCombined, time_axis, SLEStartCode, SLEEndCode);

%Removing events that are under 5 seconds long
StepEnvelopeFinal = StepEnvelopeCombined;
for a = 1:length(SLEEndIndexNew)
    if SLEEndTimeNew(a)-SLEStartTimeNew(a) < 5
        StepEnvelopeFinal( SLEStartIndexNew(a): SLEEndIndexNew(a)) = 0;
    end
end

%Creating Labels for Events Picked Up
[SLEStartTimeLabel, SLEEndTimeLabel, SLEStartIndexLabel, SLEEndIndexLabel] = findStartEndPoints(StepEnvelopeFinal, time_axis, SLEStartCode, SLEEndCode);
%% Plot the results against time

figure(1)
ax11 = subplot(3,1,1);
plot(t, CortexIpsilateral);
hold on
plot(time_axis, StepEnvelopeFinal*100, 'linewidth', 2);
hold on
for i = 1:length(SLEStartIndexLabel)
    text(time_axis(SLEStartIndexLabel(i)), 600, num2str(i))
    hold on
end
title('Step Function after combining events within 5s of eachother, and removing events under 5s');
title('Pre-Processed Data with Step Function');
xlabel('Time (s)');
legend('Neural Recordings', 'Step Function')
ylim([-1500, 1500])

ax12 = subplot(3,1,2);
plot(time_axis, teager_avg, 'linewidth', 2);
hold on
plot(time_axis, SLEThresholdPlot)
title('Average Teager Value For a Window Size of 200 samples and an Overlap of 90%');
xlabel('Time (s)');
ylim([0, 50000])

ax13 = subplot(3,1,3);
plot(t, CortexIpsi);
title('Step Function Threshold = 3.5 x Avg(minimum 50 teager values');
xlabel('Time (s)');
ylim([-1500, 1500])

% ax14 = subplot(5,1,4);
% plot(time_axis, StepEnvelopeReduced, 'linewidth', 2);
% title('Removing Low Duration Events (<0.1s)');
% xlabel('Time (s)');
% ylim([-0.5, 5.5])
% 
% ax15 = subplot(5,1,5);
% plot(time_axis, StepEnvelopeFinal, 'linewidth', 2);
% title('Step Function after combining events within 5s of eachother, and removing events under 5s');
% xlabel('Time (s)');
% ylim([-0.5, 5.5])

linkaxes([ax11, ax12, ax13], 'x')
xlim([1, 80])

%Writing to text file to export to labchart

%% Frequency Analysis of Events
% Creating the windows to perform Power Spectrum Analysis
%Calculating Start and End Points on finalized Envelope
[SLEWindowStart, SLEWindowEnd, SLEIndexStart, SLEIndexEnd] = findStartEndPoints(StepEnvelopeFinal, time_axis, SLEStartCode, SLEEndCode);
for i = 1:length(SLEWindowStart)
    Start(i) = find(t == SLEWindowStart(i));
    timeStart(i) = t(find(t == SLEWindowStart(i)));
    
    End(i) = find(t == SLEWindowEnd(i));
    timeEnd(i) = t(find(t == SLEWindowEnd(i)));
end
NFFT = 1024;
Frequency = (0:NFFT/2) *(fs/NFFT); %Just picking 1-40Hz in this array later on
MagnitudeMatrix = zeros(length(Frequency), length(SLEIndexStart)); %Pre-allocating the array


% matrixTeager = StepEnvelopeFinal';
% matrixRaw = CortexIpsilateral;
% filename = 'Day1Teager.txt';  % Name of the output file
% % Write the matrix to a text file
% dlmwrite(filename, matrixTeager, 'delimiter', ' ');
% filename = 'Day1RawIp.txt';  % Name of the output file
% % Write the matrix to a text file
% dlmwrite(filename, matrixRaw, 'delimiter', ' ');


for c = 1:length(Start)
    time = t((Start(c) : End(c)));
    CortexIpSLE = CortexIpsilateral(Start(c) : End(c));
    CortexIpRawSLE = CortexIpsi(Start(c) : End(c));
    
    CortexIpFFT = fft(CortexIpSLE, NFFT);
    CortexIpFFT = abs(CortexIpFFT(1:NFFT/2 + 1));
    CortexIpFFT(2:end-1) = 2*CortexIpFFT(2:end-1);
    
        CortexIpRawFFT = fft(CortexIpRawSLE, NFFT);
    CortexIpRawFFT = abs(CortexIpRawFFT(1:NFFT/2 + 1));
    CortexIpRawFFT(2:end-1) = 2*CortexIpRawFFT(2:end-1);
    
    MagnitudeMatrix(:,c) = CortexIpFFT;
    
    SignalLength = length(CortexIpSLE);
    
   %FFT Plots of SLEs
%         figure
%     subplot(2,1,1)
%     plot(time, CortexIpSLE)
%     subplot(2,1,2)
%     plot(Frequency,CortexIpFFT)
%     hold on
%     plot(Frequency, CortexIpRawFFT)
%     title('Frequency Spectra of SLE')
%     ylabel('|P1(f)|')
%     xlabel('Frequency(hz)')
    
end

%% Functions
function [StartTime, EndTime, StartIndex, EndIndex] = findStartEndPoints(Envelope, Time, StartCode, EndCode)
%Given an envelope, a time array, the associated neural data, a start code
%and an end code, the function outputs regions of the step function where
%values change from 0 and to 0 (inflection points) as an index.


%Find Start Index and Time
StartIndex = strfind(Envelope,StartCode);
%Start Conditions, if the first value in the dataset is a 5, add a start
%point to the data. This is in case the file starts on a seizure, so a
%value of 5 instead of starting on baseline (0) and then picking up into a
%seizure.
if Envelope(1) == 5
    StartIndex = [1, StartIndex];
end
StartTime = Time(StartIndex);


%Find End Index
EndIndex = strfind(Envelope,EndCode);
%End Conditions, if the last value in the dataset is a 5, add an end point
%to the data, this is in case the file ends on a seizure and does not drop
%down to a zero before the file ends.
if Envelope(end) == 5
    EndIndex(end+1) = length(Envelope);
end
EndTime = Time(EndIndex);
end