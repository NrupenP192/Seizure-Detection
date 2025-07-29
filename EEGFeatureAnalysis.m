clc
clear
close all

%% User Inputs
%File Loading and Preprocessing
myFolder = 'C:\Users\Nrupen Pakalapati\OneDrive\Documents\MATLAB\PreSeizure Analysis\10m Before and After Inj\After';
FileToBeAnalyzed = 1;
Fs_original = 20000;        %(Hz)                                          % Original sampling frequency (40,000 Hz)
Fs_new = 2000;              %(Hz)                                           % New sampling frequency (1,000 Hz) for downsampling, make sure this is at least 2x LPF cutoff
Fc_low = 50;               %(Hz)                                          % Lowpass filter cutoff frequency
Fc_High = 5;                %(Hz)                                          % Highpass filter cutoff frequency
GainMultiplier = 15;        %(Unit Less)                                   % Gain = 100; want things in mV, so 1000/100 = 10

%Teager Inputs
TeagerChannel = 1;                                                         % What channel is being used for teager event detection
window_sizeTeager = 500;      %(Samples)                                   % Fs_new = 1s of window, this changes resolution of teager

NoiseMultiplier = 200000000000;  %(Unit Less)                                 % Threshold for noise, needs graphs to visualize
SLEMultiplier = 5;            %(Unit Less)                                 % Threshold for events, needs graphs to visualize; ex: 2 = 2x baseline

Teager_Smoothing = 1;         %(Binary)                                    % Binary Value, 1 = perform smoothing, 0 = ignore smoothing
Teager_SmoothingFactor = 0.4; %(Unit Less)

CombinedEventsTime = 1;       %(s)                                         % How far apart can two consequtive events be, before they are considered the same event
MinimumEventDuration = 0.5;     %(s)                                         % Removes momentary spikes in teager, which usually are artifacts


%Time Series
WindowSeconds = 100;                                  %(s)                % Moving Windows across the whole file, window width is set to window seconds, with 90% overlap
Frequency_Bands = [1 4; 5 12; 8 13; 13 30; 30 100];   %(Hz)               % Power of frequency bands; example: [band_1_low band_1_high; band_2_low band_2_high...]

%Detected Event Analysis
%Cross Correlation of Teager Detected Events
CrossCorrelationEnvelopeWindow = 80;                   %(Samples)         % Uses envelope before computing cross correlation of detected events to improve accuracy

%Information Loss: Cross Correlation of Teager Detected Events, specifically high duration events
ComputeInformationLoss = 1;                             % (Binary)
SeizureMinimumEventDuration = 20;               % (s)                      % Modifies teager function to enable a second calculation for large duration events, which are likely to be seizures

%Change in Information: Measure change in information during high energy events, assuming propagation direction changes mid seizure
ComputeChangeInInformation = 1;                 % (Binary)                 
SeizureMovingWindowDuration = 1;                % (s)                      % Moving window to analyze changes in delays within a seizure with 0% overlap
SeizureCrossCorrelationEnvelope = 20;          % (Samples)                 % Uses envelope before computing cross correlation of seizures to improve accuracy

%Filters cross correlation values based on R values
CrossCorrelation_RFilter = 0.7;   %Integer (0-1)                           % Minimum value of R before it is considered for delays
FocusChannel = 1;                                                          % Channel Number for reference

% Baseline Analysis
PerformBaselineAnalysis = 0;
Pre_Post_SeizureWindowLength = 5;
% Plot Controls
Plot_PreprocessedData_Teager = 1;   %(Binary)
Plot_TimeSeriesAnalysis = 0;        %(Binary)
Plot_EventDetectionAnalysis = 1;    %(Binary)
Plot_XCorrAnalysis = 1;             %(Binary)
%% Working Code
%Saving User Inputs
UserInputs.PreProcessing = [Fs_original, Fs_new, Fc_low, Fc_High, GainMultiplier]';
UserInputs.Teager = [TeagerChannel,window_sizeTeager,NoiseMultiplier,SLEMultiplier,Teager_Smoothing,Teager_SmoothingFactor,CombinedEventsTime,MinimumEventDuration]';
UserInputs.TimeSeriesAnalysis = {WindowSeconds, Frequency_Bands, CrossCorrelationEnvelopeWindow}';
UserInputs.DetectedEventAnalysis = [ComputeInformationLoss, SeizureMinimumEventDuration, FocusChannel, SeizureMovingWindowDuration, SeizureCrossCorrelationEnvelope, ComputeChangeInInformation, CrossCorrelation_RFilter]';
UserInputs.BaselineAnalysis = [Pre_Post_SeizureWindowLength, PerformBaselineAnalysis];
UserInputs.Plots = [Plot_PreprocessedData_Teager, Plot_TimeSeriesAnalysis, Plot_EventDetectionAnalysis]';


[FilePath, BaseNames] = FileLoader(myFolder);
for jj = FileToBeAnalyzed:FileToBeAnalyzed   %Change This for the file                 

    %File Loading and Preprocessing
    filename = fullfile(FilePath(jj), BaseNames(jj));
    fprintf('Currently Reading: %s\n', cell2mat(BaseNames(jj)));

    tic;
    Data = readmatrix(filename{1});                                        %reads all channels on file, you need to know what each channel is
    [PreProcessedData] = PreProcess(Data, UserInputs);
    elapsedTime = toc;
    disp(['Pre Processing: ', num2str(elapsedTime), ' seconds']);

    %Time Series (TS) Windows without Event Detection
    tic;
    [TimeSeriesFeatures] = TimeSeriesAnalysis(PreProcessedData, UserInputs);
    elapsedTime = toc;
    disp(['Time Series Analysis ', num2str(elapsedTime), ' seconds']);

    %Teager Event Identification
    tic;
    [TeagerFeatures, DetectedEventTimes] = Teager(PreProcessedData, UserInputs);
    elapsedTime = toc;
    disp(['Teager Detection ', num2str(elapsedTime), ' seconds']);

    %Analysis using Teager Event Detection (ED)
    tic;
    [DetectedEventFeatures] = DetectedEventAnalysis( ...
        DetectedEventTimes, PreProcessedData, UserInputs);
    elapsedTime = toc;
    disp(['Detected Event Analysis ', num2str(elapsedTime), ' seconds']);
    
    %Loss of information during seizures (Hypothesis: CC does not have high
    %frequency components of seizures)
    tic;
    [LossOfInformation] = LossCalculation(UserInputs, DetectedEventFeatures);
    elapsedTime = toc;
    disp(['Loss of Information ', num2str(elapsedTime), ' seconds']);

    %Change in information during seizures (Hypothesis: Different delays
    %occur during the whole seizure event)(Hypothesis: Different delays are
    %related to RMS)
    tic;
    [ChangeInInformationDuringSeizures] = ChangeInInformation_SLES(UserInputs, DetectedEventTimes, PreProcessedData);
    elapsedTime = toc;
    disp(['Change In Information ', num2str(elapsedTime), ' seconds']);

    %Filters CrossCorrelation Values from Event Detection and Change in
    %information based on user defined R value.
    tic;
    [FilteredXCorr] = FilterXcorr(DetectedEventFeatures, ChangeInInformationDuringSeizures, UserInputs);
    elapsedTime = toc;
    disp(['Filtered XCorr ', num2str(elapsedTime), ' seconds']);

    [BaselineEventFeatures] = BaselineAnalysis(UserInputs, DetectedEventTimes, PreProcessedData);
    
    % %Band Power for whole file for normalization purposes
    % for k = 1:4
    %  [BandPowerWholeFile] = CalculateFrequencyBandPower(PreProcessedData.AllChannels(:, k), UserInputs.PreProcessing(2), UserInputs.TimeSeriesAnalysis{2});
    %  FrequencyBandPower.(sprintf('Ch%d', k)) = BandPowerWholeFile;
    % end
end
clearvars -except PreProcessedData TeagerFeatures TimeSeriesFeatures DetectedEventFeatures UserInputs DetectedEventTimes LossOfInformation ChangeInInformationDuringSeizures FilteredXCorr BaselineEventFeatures FrequencyBandPower
PlotEvents(PreProcessedData, TeagerFeatures, TimeSeriesFeatures, DetectedEventFeatures, UserInputs);
PreProcessedData.Time(end)

%% Function 10: Baseline Analysis. Definition: Everything in between seizures
% seizures is considered a baseline. Hypothesis, RMS and other features increases between
% seizures
function [BaselineEventFeatures] = BaselineAnalysis(UserInputs, DetectedEventTimes, PreProcessedData)
if UserInputs.BaselineAnalysis(2) == 0
    BaselineEventFeatures = [];
    return
else
    PreProcessedData_AllChannels = PreProcessedData.AllChannels;
    timeVector_original = PreProcessedData.Time;
    Fs = UserInputs.PreProcessing(2);
    WindowLength = UserInputs.BaselineAnalysis(1); %looking at 10s of baseline right after a seizure ends, and 10s before the start of a seizure
    WindowLength = WindowLength*Fs;

    BaselineEvent_SeizureEnd_Start = DetectedEventTimes{3};
    BaselineEvent_SeizureEnd_End = BaselineEvent_SeizureEnd_Start + WindowLength; %seizure cant be at the end of the file

    BaselineEvent_SeizureStart_End = DetectedEventTimes{4};
    BaselineEvent_SeizureStart_Start = BaselineEvent_SeizureStart_End - WindowLength; %seizure cant begin at time 0
    %SeizureEnd_Start:SeizureEnd_End should theoritically be the first 10s
    %right after a seizure, and vice versa for the other one

    NumEvents = length(BaselineEvent_SeizureEnd_Start);
    NumChannels = size(PreProcessedData_AllChannels, 2);

    % Preallocate Matrices
    LineLength_SeizureEnd = zeros(NumEvents, NumChannels);
    RMS_SeizureEnd = zeros(NumEvents, NumChannels);
    LineLength_SeizureStart = zeros(NumEvents, NumChannels);
    RMS_SeizureStart = zeros(NumEvents, NumChannels);
    DetectedEventDuration_SeizureEnd = BaselineEvent_SeizureEnd_End' - BaselineEvent_SeizureEnd_Start'; %dummy check, should be equal to window length


    %Analyzing the first 10s after a seizure
    for i = 1:NumEvents
        Channel_Events_SeizureEnd = PreProcessedData_AllChannels(BaselineEvent_SeizureEnd_Start(i):BaselineEvent_SeizureEnd_End(i), :);%all channels
        Channel_Time_SeizureEnd = timeVector_original(BaselineEvent_SeizureEnd_Start(i):BaselineEvent_SeizureEnd_End(i));

        Channel_Events_SeizureStart = PreProcessedData_AllChannels(BaselineEvent_SeizureStart_Start(i):BaselineEvent_SeizureStart_End(i), :);%all channels

        % Line Length Calculation
        dV_End = diff(Channel_Events_SeizureEnd);  % Differences between consecutive points in Channel_Events
        dV_Start = diff(Channel_Events_SeizureStart);

        dt = diff(Channel_Time_SeizureEnd);    % Differences between consecutive points in Channel_Time
        dt = dt(1);  % Take the first value of dt (since it's constant)
        LineLength_SeizureEnd(i,:) = sum(sqrt(dV_End.^2 + dt^2), 1);  % Compute distances between consecutive points for each channel
        LineLength_SeizureStart(i,:) = sum(sqrt(dV_Start.^2 + dt^2), 1);  % Compute distances between consecutive points for each channel

        % RMS Calculation
        RMS_SeizureEnd(i,:) = rms(Channel_Events_SeizureEnd);
        RMS_SeizureStart(i,:) = rms(Channel_Events_SeizureStart);
    end
    % Store all calculated features in the output structure
    BaselineEventFeatures.LineLength_SLE_End = LineLength_SeizureEnd;
    BaselineEventFeatures.RMS_SLE_End = RMS_SeizureEnd;
    BaselineEventFeatures.LineLength_SLE_Start = LineLength_SeizureStart;
    BaselineEventFeatures.RMS_SLE_Start = RMS_SeizureStart;
    BaselineEventFeatures.SLE_End = DetectedEventDuration_SeizureEnd/Fs; %in s
end
end
%% Function 9: Filter data based on XCorr R values
function [FilteredXcorr] = FilterXcorr(DetectedEventFeatures, ChangeInInformationDuringSeizures, UserInputs)
DetectedXCorr_R = DetectedEventFeatures.CrossCorrelation_RMax;
DetectedXCorr_Delays = DetectedEventFeatures.CrossCorrelation_Delays;
ChangeInformationXCorr_R = ChangeInInformationDuringSeizures.CrossCorrelation_RMax;
ChangeInformationXCorr_Delay = ChangeInInformationDuringSeizures.CrossCorrelation_Delays;
CrossCorrelation_RFilter = UserInputs.DetectedEventAnalysis(7);
[NumEvents, NumCombinations] = size(DetectedXCorr_Delays);

% Process Detected Events Delays
DetectedXCorr_Delays_Filtered_Cell = cell(1, NumCombinations); % Initialize as cell array
maxLengthDetected = 0;
for i = 1:NumCombinations
    Indices = find(DetectedXCorr_R(:,i)>CrossCorrelation_RFilter);
    DetectedXCorr_Delays_Filtered_Cell{1,i} = DetectedXCorr_Delays(Indices,i);
    maxLengthDetected = max(maxLengthDetected, length(DetectedXCorr_Delays_Filtered_Cell{1,i}));
end

% Pad Detected Events Delays to create a matrix
PaddedDetectedDelays = NaN(NumCombinations, maxLengthDetected);
for i = 1:NumCombinations
    currentLength = length(DetectedXCorr_Delays_Filtered_Cell{1,i});
    PaddedDetectedDelays(i, 1:currentLength) = DetectedXCorr_Delays_Filtered_Cell{1,i};
end

FilteredXcorr.DetectedEvents = PaddedDetectedDelays';

% Process Change In Information Delays
fieldList = fieldnames(ChangeInformationXCorr_R);
numChannelCombinations = numel(fieldList);
Delay_ChangeInInformation = struct(); % Initialize as an empty struct

for i = 1:numChannelCombinations
    CurrentCombination = fieldList{i};
    R_TimeWindows = ChangeInformationXCorr_R.(CurrentCombination);
    Delay_TimeWindows = ChangeInformationXCorr_Delay.(CurrentCombination);
    [NumEventsTW, EventTimeWindow] = size(R_TimeWindows);
    TimeWindow_Delays_Filtered_Cell = cell(1, EventTimeWindow); % Initialize as cell array
    maxLengthTW = 0;
    for j = 1:EventTimeWindow
        Indices_TimeWindow = find(R_TimeWindows(:,j)>CrossCorrelation_RFilter);
        TimeWindow_Delays_Filtered_Cell{1,j} = Delay_TimeWindows(Indices_TimeWindow, j);
        maxLengthTW = max(maxLengthTW, length(TimeWindow_Delays_Filtered_Cell{1,j}));
    end

    % Pad Time Window Delays to create a matrix
    PaddedTimeWindowDelays = NaN(EventTimeWindow, maxLengthTW);
    for j = 1:EventTimeWindow
        currentLengthTW = length(TimeWindow_Delays_Filtered_Cell{1,j});
        PaddedTimeWindowDelays(j, 1:currentLengthTW) = TimeWindow_Delays_Filtered_Cell{1,j};
    end
    Delay_ChangeInInformation.(CurrentCombination) = PaddedTimeWindowDelays';
end
FilteredXcorr.ChangeInformation = Delay_ChangeInInformation;

end
%% Function 8: Change in Information During Seizures
function [ChangeInInformationDuringSeizures] = ChangeInInformation_SLES(UserInputs, DetectedEventTimes, PreProcessedData)
if UserInputs.DetectedEventAnalysis(5) == 0
    return
else
    % User Inputs
    MinimumEventDuration = UserInputs.DetectedEventAnalysis(2); % Seconds
    Fs = UserInputs.PreProcessing(2);
    CrossCorrelationMovingWindow = UserInputs.DetectedEventAnalysis(4); % Seconds
    CrossCorrelationEnvelopeWindow = UserInputs.DetectedEventAnalysis(5); % For crosscorrelation calculation

    % Pre-Processed Data
    PreProcessedData_AllChannels = PreProcessedData.AllChannels;
    timeVector_original = PreProcessedData.Time;

    % Teager Detected Event Start and End Times
    TeagerEvent_Start = DetectedEventTimes{1};
    TeagerEvent_End = DetectedEventTimes{2};
    EventDuration = TeagerEvent_End - TeagerEvent_Start;
    MinimumEventDuration_Samples = MinimumEventDuration * Fs; % User Input converted to minimum sample difference

    % Find SLEs based on Minimum Event Duration
    SLE_Indices = EventDuration > MinimumEventDuration_Samples;
    SLE_Start = TeagerEvent_Start(SLE_Indices);
    SLE_End = TeagerEvent_End(SLE_Indices);

    CrossCorrelationWindow_SLE = CrossCorrelationMovingWindow * Fs; % Samples

    % Finding combinations, events and channels
    NumEvents = length(SLE_Start);
    NumChannels = size(PreProcessedData_AllChannels, 2);
    NumCombinations = NumChannels * (NumChannels - 1) / 2;

    % Generate channel combinations (all pairs)
    channelCombinations = nchoosek(1:NumChannels, 2); % Get all combinations of channel pairs
    ChangeInInformationDuringSeizures.ChannelCombinations = channelCombinations; % Store combinations in the output structure

    % Preallocate arrays based on Minimum EventDuration and Chosen CrossCorrelation Window
    NumWindows = floor(MinimumEventDuration_Samples / CrossCorrelationWindow_SLE);

    % Preallocate structures for Rmax, Delays, and RMS
    for combIdx = 1:NumCombinations
        ch1 = channelCombinations(combIdx, 1);
        ch2 = channelCombinations(combIdx, 2);
        fieldName = sprintf('Ch%d_Ch%d', ch1, ch2);
        ChangeInInformationDuringSeizures.CrossCorrelation_RMax.(fieldName) = zeros(NumEvents, NumWindows);
        ChangeInInformationDuringSeizures.CrossCorrelation_Delays.(fieldName) = zeros(NumEvents, NumWindows);
    end

    for ch = 1:NumChannels
        channelName = sprintf('Ch%d', ch);
        ChangeInInformationDuringSeizures.RMS.(channelName) = zeros(NumEvents, NumWindows);
    end

    for EventNum = 1:NumEvents % for each event
        Channel_Events = PreProcessedData_AllChannels(SLE_Start(EventNum):min(SLE_Start(EventNum) + MinimumEventDuration_Samples - 1, SLE_End(EventNum)), :); % first MinimumEventDuration or whole event if shorter
        Channel_Time = timeVector_original(SLE_Start(EventNum):min(SLE_Start(EventNum) + MinimumEventDuration_Samples - 1, SLE_End(EventNum)));
        Channel_Envelopes = envelope(Channel_Events, CrossCorrelationEnvelopeWindow, 'peak');

        for windowIdx = 1:NumWindows % Loop through each window in each combination for each event for cross correlation
            startSample = (windowIdx - 1) * CrossCorrelationWindow_SLE + 1;
            endSample = min(windowIdx * CrossCorrelationWindow_SLE, size(Channel_Envelopes, 1)); % prevent going out of bounds
            
            for combIdx = 1:NumCombinations
                ch1 = channelCombinations(combIdx, 1);
                ch2 = channelCombinations(combIdx, 2);
                fieldName = sprintf('Ch%d_Ch%d', ch1, ch2);
                windowDataCh1 = Channel_Envelopes(startSample:endSample, ch1);
                windowDataCh2 = Channel_Envelopes(startSample:endSample, ch2);
                [c, lags] = xcorr(windowDataCh1, windowDataCh2, 'coeff');
                [Rmax, maxIdx] = max(abs(c));
                dt = diff(Channel_Time(1:2));
                bestLag = 1000 * dt * lags(maxIdx);
                ChangeInInformationDuringSeizures.CrossCorrelation_RMax.(fieldName)(EventNum, windowIdx) = Rmax;
                ChangeInInformationDuringSeizures.CrossCorrelation_Delays.(fieldName)(EventNum, windowIdx) = bestLag;
            end

            for ch = 1:NumChannels
                channelName = sprintf('Ch%d', ch);
                windowDataCh = Channel_Events(startSample:endSample, ch);
                RMS_window = rms(windowDataCh);
                ChangeInInformationDuringSeizures.RMS.(channelName)(EventNum, windowIdx) = RMS_window;
            end
        end
    end
end
end
%% Function 7: Loss of information in detected events
    function [LossOfInformation] = LossCalculation(UserInputs, DetectedEventFeatures)
        if UserInputs.DetectedEventAnalysis(1) == 0
            LossOfInformation = 'NA';
        else
            %User Inputs
            MinimumEventDuration = UserInputs.DetectedEventAnalysis(2); %Seconds

            %Data Extraction
            Time = DetectedEventFeatures.Time;
            LineLength = DetectedEventFeatures.LineLength;
            Entropy = DetectedEventFeatures.Entropy;
            RMS = DetectedEventFeatures.RMS;
            BandPowers = DetectedEventFeatures.BandPowers;
            EventDuration = DetectedEventFeatures.EventDuration;
            Fs = UserInputs.PreProcessing(2);
            

            % Identify SLEs based on Minimum Event Duration
            SLE_Indices = EventDuration > MinimumEventDuration;

            % Extract features based on identified SLEs
            SLE_RMS = RMS(SLE_Indices, :);
            SLE_Entropy = Entropy(SLE_Indices, :);
            SLE_LineLength = LineLength(SLE_Indices, :);

            % Extract BandPowers for each channel based on identified SLEs
            NumChannels = length(fieldnames(BandPowers));
            for ch = 1:NumChannels
                channelField = sprintf('Ch%d', ch);  % Channel field name (e.g., 'Ch1', 'Ch2', ...)
                % Extract the band powers for the current channel for all events
                SLE_BandPowers.(channelField) = BandPowers.(channelField)(SLE_Indices, :);
                % Store the extracted band powers in an array with frequency bands as rows and events as columns
                % The matrix is automatically constructed when you index the relevant events using SLE_Indices
            end
            LossOfInformation.RMS = SLE_RMS;
            LossOfInformation.Entropy = SLE_Entropy;
            LossOfInformation.LineLength = SLE_LineLength;
            LossOfInformation.BandPowers = SLE_BandPowers;
        end
    end
%% Function 6: Frequency Spectra Calculation using FFT
function [ChannelBandPower] = CalculateFrequencyBandPower(data_matrix, fs, Frequency_Bands)
[N, M] = size(data_matrix); % Get matrix size
[NumFrequencyBands, ~] = size(Frequency_Bands);

% Compute one-sided FFT for all windows
fft_result = fft(data_matrix);
f = (0:(N/2)) * (fs / N);

 if M == 1
        fft_magnitude = abs(fft_result(1:floor(N/2+1))) / N; % Normalize (no colon needed for single column)
    else
        fft_magnitude = abs(fft_result(1:N/2+1, :)) / N; % For Multiple columns in Time Series Analysis
 end

fft_magnitude(2:end-1, :) = 2 * fft_magnitude(2:end-1, :);
idx_limit = f <= 100;  % Find indices where frequency is <= 100 Hz
frequencies = f(idx_limit);  % Limit frequency axis
fft_magnitude = fft_magnitude(idx_limit, :);  % Limit magnitude spectrum
ChannelBandPower = zeros(M, NumFrequencyBands);

for k = 1:NumFrequencyBands
    Band = (frequencies >= Frequency_Bands(k,1) & frequencies < Frequency_Bands(k,2));
    ChannelBandPower(:, k) = sum(fft_magnitude(Band, :), 1);
end
end
%% Function 5: Detected Event Analysis
function[DetectedEventFeatures] = DetectedEventAnalysis(DetectedEventTimes, PreProcessedData, UserInputs)
UserInputs_Analysis = UserInputs.TimeSeriesAnalysis;
PreProcessedData_AllChannels = PreProcessedData.AllChannels;
timeVector_original = PreProcessedData.Time;

TeagerEvent_Start = DetectedEventTimes{1};
TeagerEvent_End = DetectedEventTimes{2};
CrossCorrelationEnvelopeWindow = UserInputs_Analysis{3};
FrequencyBands = UserInputs_Analysis{2};
Fs = UserInputs.PreProcessing(2);

NumEvents = length(TeagerEvent_Start);
NumChannels = size(PreProcessedData_AllChannels, 2);
NumCombinations = NumChannels*(NumChannels-1)/2;

% Generate channel combinations (all pairs)
channelCombinations = nchoosek(1:NumChannels, 2);  % Get all combinations of channel pairs
DetectedEventFeatures.ChannelCombinations = channelCombinations;  % Store combinations in the output structure

% Preallocate Matrices
RmaxMatrix = zeros(NumEvents, NumCombinations);
DelayMatrix = zeros(NumEvents, NumCombinations);
LineLength = zeros(NumEvents, NumChannels);
Entropy = zeros(NumEvents, NumChannels);
RMS_DetectedEvents = zeros(NumEvents, NumChannels);
Time_TeagerMidPoint = timeVector_original((TeagerEvent_Start + (TeagerEvent_End - TeagerEvent_Start))');
DetectedEventDuration = TeagerEvent_End' - TeagerEvent_Start';
BandPowers = cell2struct(cell(NumChannels, 1), ...
arrayfun(@(ch) sprintf('Ch%d', ch), 1:NumChannels, 'UniformOutput', false));  % Create field names dynamically

for i = 1:NumEvents
    Channel_Events = PreProcessedData_AllChannels(TeagerEvent_Start(i):TeagerEvent_End(i), :);%all channels
    Channel_Time = timeVector_original(TeagerEvent_Start(i):TeagerEvent_End(i));
    Channel_Envelopes = envelope(Channel_Events, CrossCorrelationEnvelopeWindow, 'peak');

    % Line Length Calculation
    dV = diff(Channel_Events);  % Differences between consecutive points in Channel_Events
    dt = diff(Channel_Time);    % Differences between consecutive points in Channel_Time
    dt = dt(1);  % Take the first value of dt (since it's constant)
    LineLength(i,:) = sum(sqrt(dV.^2 + dt^2), 1);  % Compute distances between consecutive points for each channel

    % RMS Calculation
    RMS_DetectedEvents(i,:) = rms(Channel_Events);

    % Cross Correlation and Entropy
    combinationIndex = 1;
    for ch = 1:NumChannels
        % Entropy Calculation
        ChannelData = Channel_Events(:,ch); %takes a single channel
        NormalizedChannelData = (ChannelData-min(ChannelData)) / (max(ChannelData) - min(ChannelData));
        [counts, ~] = histcounts(NormalizedChannelData, 'Normalization', 'probability');
        counts = counts(counts > 0); % Removing zero probabilities to avoid log(0)
        Entropy(i, ch) = -sum(counts .* log2(counts));

        % Spectra Calculation
        %for every event, look at every channel and do the spectra and
        %store the spectra in an array for all events, where the columns are events, rows
        %are frequency bands. Then store these arrays in a struct for each
        %channel. Then store all the channel structs in another struct
        %called BandPowers

        % Spectra Calculation using CalculateFrequencyBandPower
        bandPowers = CalculateFrequencyBandPower(ChannelData, Fs, FrequencyBands);

        % Store the spectra in a struct for the current channel
        BandPowers.(sprintf('Ch%d', ch))(i,:) = bandPowers;  % Store spectra (frequency band powers) for the current channel and event

        for jj = ch+1:NumChannels  % Nested loop for all channel combinations for cross-correlation
            [c, lags] = xcorr(Channel_Envelopes(:, ch), Channel_Envelopes(:, jj), 'coeff');
            [Rmax, maxIdx] = max(abs(c));  % Get the max correlation (absolute value)
            bestLag = 1000 * dt * lags(maxIdx); % Get the corresponding lag for Rmax, converted to ms using dt
            RmaxMatrix(i, combinationIndex) = Rmax;
            DelayMatrix(i, combinationIndex) = bestLag;
            combinationIndex = combinationIndex + 1;
        end
    end
end

% Store all calculated features in the output structure
DetectedEventFeatures.Time = Time_TeagerMidPoint;
DetectedEventFeatures.CrossCorrelation_RMax = RmaxMatrix;
DetectedEventFeatures.CrossCorrelation_Delays = DelayMatrix;
DetectedEventFeatures.ChannelCombinations = channelCombinations;
DetectedEventFeatures.LineLength = LineLength;
DetectedEventFeatures.Entropy = Entropy;
DetectedEventFeatures.RMS = RMS_DetectedEvents;
DetectedEventFeatures.EventDuration = DetectedEventDuration/Fs; %in s
DetectedEventFeatures.BandPowers = BandPowers;  % Store BandPowers struct in the output structure
end
%% Function 4: Time Series Analysis
function[TimeSeriesFeatures] = TimeSeriesAnalysis(PreProcessedData, UserInputs)
UserInputs_Analysis = UserInputs.TimeSeriesAnalysis;
UserInputs_PreProcessing = UserInputs.PreProcessing;

PreProcessedData_AllChannels = PreProcessedData.AllChannels;
timeVector_original = PreProcessedData.Time;
WindowSeconds = UserInputs_Analysis{1};
Frequency_Bands = UserInputs_Analysis{2};
Fs_new = UserInputs_PreProcessing(2);

window_size = WindowSeconds*Fs_new;
window_overlap = round(window_size/1.1111); % 90% overlap
dt = diff(timeVector_original);
dt = dt(1);  % Take the first value of dt (since it's constant)
NumChannels = width(PreProcessedData_AllChannels);
ChannelWindows_TimeSeries = struct();
ChannelTime_TimeSeries = buffer((timeVector_original), window_size, window_overlap, 'nodelay');
ChannelTime_TimeSeries = ChannelTime_TimeSeries((window_size/2),:); %taking the midpoint of each

for i = 1:NumChannels
    windowsPerChannel = buffer((PreProcessedData_AllChannels(:,i)), window_size, window_overlap, 'nodelay');
    ChannelWindows_TimeSeries.(sprintf('Ch%d', i)) = windowsPerChannel;

    %Time Series RMS
    RMS_TimeSeries.(sprintf('Ch%d' ,i)) = (rms(windowsPerChannel, 1))';

    %Time Series Line Length
    dV = diff(windowsPerChannel);  % Differences between consecutive points in Channel_Events (4 channels)
    % Assuming dt is constant, calculate the line length for each channel
    LineLength_TimeSeries.(sprintf('Ch%d' ,i)) = (sum(sqrt(dV.^2 + dt^2), 1))';  % Compute distances between consecutive points for each channel

    %Time Series Teager
    Tn = 2; %probably dont have a reason to change this value
    %Teager discrete formula where n can be changed, and is vectorized: y = x(n+1:end-n).^2 - x(1:end-2n) . x(2*n+1:end); %general formula
    teager_values = (sum(windowsPerChannel(Tn+1:end-Tn, :).^2 - windowsPerChannel(1:end-2*Tn,:).*windowsPerChannel(2*Tn+1:end,:),1));
    Teager_TimeSeries.(sprintf('Ch%d', i)) = normalizeData((teager_values./window_size)');

    % Time Series Entropy
    EntropyValues = zeros(1, size(windowsPerChannel, 2));  % Initialize the entropy array for each window
    for j = 1:size(windowsPerChannel, 2)
        data = windowsPerChannel(:, j);
        data = (data - min(data)) / (max(data) - min(data));  % Normalize to [0, 1]
        [counts, ~] = histcounts(data, 'Normalization', 'probability', 'BinMethod', 'auto');
        counts = counts(counts > 0); %Removing zeros, to perform log
        EntropyValues(j) = -sum(counts .* log2(counts));
    end
    Entropy_TimeSeries.(sprintf('Ch%d', i)) = EntropyValues';

    %Time Series Spectra
    [ChannelBandPower] = CalculateFrequencyBandPower(windowsPerChannel, Fs_new, Frequency_Bands);
    FrequencyBandPower_TimeSeries.(sprintf('Ch%d', i)) = ChannelBandPower;
end

TimeSeriesFeatures.FeatureTime =  (ChannelTime_TimeSeries/60)';
TimeSeriesFeatures.RMS = RMS_TimeSeries;
TimeSeriesFeatures.Teager = Teager_TimeSeries;
TimeSeriesFeatures.LineLength = LineLength_TimeSeries;
TimeSeriesFeatures.Entropy = Entropy_TimeSeries;
TimeSeriesFeatures.FrequencyBandPower = FrequencyBandPower_TimeSeries;
end
function normData = normalizeData(data)
% Normalize data between 0 and 1
minVal = min(data(:));
maxVal = max(data(:));
normData = (data - minVal) / (maxVal - minVal);
end
%% Function 3: Teager Windowing and Detection
function[TeagerFeatures, DetectedEventTimes] = Teager(PreProcessedData, UserInputs)
UserInputs_Teager = UserInputs.Teager;
Data = PreProcessedData.AllChannels;
timeVector_original = PreProcessedData.Time;

TeagerChannel = UserInputs_Teager(1);
windowsize_Teager = UserInputs_Teager(2);
NoiseMultiplier = UserInputs_Teager(3);
SLEMultiplier = UserInputs_Teager(4);
Teager_Smoothing = UserInputs_Teager(5);
Teager_SmoothingFactor = UserInputs_Teager(6);
CombinedEventsTime = UserInputs_Teager(7);
MinimumEventDuration = UserInputs_Teager(8);

window_overlap = round(windowsize_Teager/1.1111); % 90% overlap; example:
% if window size is 100 samples, this moves 10 samples to the right
% before computing teager and so on till the end of the dataset.
% But to do the calculation for teager, it uses the 100 samples
% per window. So 0:100, 10:110, 20:120 and so on. To change overlap,
% change the 1.1111 repeating and recalculate.

% Compute the windows for teager analysis, each window is in a column
windows = buffer((Data(:, TeagerChannel)), windowsize_Teager, window_overlap, 'nodelay'); %creates a matrix where each column is the example above
windowsT = buffer(timeVector_original, windowsize_Teager, window_overlap, 'nodelay');

% Compute the Teager non-linear energy values
%Teager Modifiers
Tn = 2; %probably dont have a reason to change this value
%Teager discrete formula where n can be changed, and is vectorized: y = x(n+1:end-n).^2 - x(1:end-2n) . x(2*n+1:end); %general formula
teager_values = (sum(windows(Tn+1:end-Tn, :).^2 - windows(1:end-2*Tn,:).*windows(2*Tn+1:end,:),1));
teager_avg = teager_values./windowsize_Teager;
time_axis = windowsT((windowsize_Teager/2),:); %taking the midpoint of each
% window as the index for the teager value. Using the same example
% above, 0:100 would identify 50 as the midpoint

%Teager smoothing
%For calculations like crosscorrelation, you need to have sufficient
%padding on either side of the chosen event for accurate measurements.
%To do this, an envelope is created and then the envelope is smoothed
%to remove more abrupt changes. This widens the area of interest. For
%example, if a seizure is observed from 10s:30s, this will pick out
%around 8s:32s. This is unnecessary and probably gets in the way for
%seizure duration measurements. For those, just remove the envelope.
%The smoothing still helps remove abrupt changes. Smoothing factors
%might need to be tweaked.
if Teager_Smoothing == 1
    teager_avg = smoothdata(teager_avg,'movmean','SmoothingFactor',Teager_SmoothingFactor,'SamplePoints',time_axis);
end

% Calculating the Step Function for seizure detection Based on the Teager values
minteager = mink(abs(teager_avg), round(length(teager_avg)*0.1)); % Takes the minimum 10% of teager values in the dataset and uses them to threshold
avg_min = mean(minteager);
SLEThreshold = avg_min*SLEMultiplier;
NoiseThreshold = avg_min*NoiseMultiplier;

%Creating SLE and Noise Threshold Line Graphs
SLEThresholdPlot = SLEThreshold + zeros(1,length(time_axis));
NoiseThresholdPlot = NoiseThreshold + zeros(1,length(time_axis));

%Thresholding Teager
indexSLE = teager_avg >= SLEThreshold & teager_avg <= NoiseThreshold;
StepEnvelope = zeros(1, length(teager_avg));


StepFunctionAmplitude = 5;
StepEnvelope(indexSLE) = StepFunctionAmplitude;
SLEStartCode = [0, StepFunctionAmplitude];
SLEEndCode = [StepFunctionAmplitude, 0];

% (1)
% IF TWO EVENTS OCCUR WITHIN 0.5 SECONDS OF EACHOTHER, THEY ARE CONSIDERED TO
% BE THE SAME EVENT
[SLEStartTime, SLEEndTime, SLEStartIndex, SLEEndIndex] = findStartEndPoints(StepEnvelope, time_axis, SLEStartCode, SLEEndCode);
StepEnvelopeCombined = StepEnvelope;
for b = 1:length(SLEEndIndex)-1
    if SLEStartTime(b+1) - SLEEndTime(b) <= CombinedEventsTime
        StepEnvelopeCombined( SLEEndIndex(b): SLEStartIndex(b+1) ) = StepFunctionAmplitude;
    end
end

%(2)
% IF AN EVENT'S DURATION AFTER BEING COMBINED IS UNDER 0.01 SECONDS, IT IS
% CONSIDERED TO NOT BE OF INTEREST. THIS ALSO GETS RID OF SOME
% ARTIFACTS
[SLEStartTimeNew, SLEEndTimeNew, SLEStartIndexNew, SLEEndIndexNew] = findStartEndPoints(StepEnvelopeCombined, time_axis, SLEStartCode, SLEEndCode);
StepEnvelopeReduced = StepEnvelopeCombined;
for a = 1:length(SLEEndIndexNew)
    if SLEEndTimeNew(a)-SLEStartTimeNew(a) < MinimumEventDuration
        StepEnvelopeReduced( SLEStartIndexNew(a): SLEEndIndexNew(a)) = 0;
    end
end

% Teager Finalized Events
SLE_DetectedEvents = StepEnvelopeReduced;
[SLEWindowStart, SLEWindowEnd, ~, ~] = findStartEndPoints(StepEnvelopeReduced, time_axis, SLEStartCode, SLEEndCode);
TeagerEvent_Start = arrayfun(@(i) find(timeVector_original == SLEWindowStart(i)), 1:length(SLEWindowStart));
TeagerEvent_End = arrayfun(@(i) find(timeVector_original == SLEWindowEnd(i)), 1:length(SLEWindowEnd));

% All events in between teager marked events (assumed to be baseline)
Baseline_DetectedEvents = abs(StepEnvelopeReduced-5);
[SLEWindowStartRMS, SLEWindowEndRMS, ~, ~] = findStartEndPoints(Baseline_DetectedEvents, time_axis, SLEStartCode, SLEEndCode);
Baseline_Start = arrayfun(@(i) find(timeVector_original == SLEWindowStartRMS(i)), 1:length(SLEWindowStartRMS));
Baseline_End = arrayfun(@(i) find(timeVector_original == SLEWindowEndRMS(i)), 1:length(SLEWindowEndRMS));

TeagerFeatures.TeagerAvg = teager_avg';
TeagerFeatures.TeagerTime = time_axis';
TeagerFeatures.DetectedEvents = SLE_DetectedEvents';
TeagerFeatures.BaselineEvents = Baseline_DetectedEvents';
TeagerFeatures.SLEThresholdPlot = SLEThresholdPlot';
TeagerFeatures.NoiseThresholdPlot = NoiseThresholdPlot';

DetectedEventTimes = {TeagerEvent_Start, TeagerEvent_End, Baseline_Start, Baseline_End}; %Input to Detected Event Analysis
end
function [StartTime, EndTime, StartIndex, EndIndex] = findStartEndPoints...% Helper Function for Teager Event Detection
    (Envelope, Time, StartCode, EndCode)
%Given an envelope, a time array, the associated neural data, a start code
%and an end code, the function outputs regions of the step function where
%values change from 0 (inflection points) as an index.


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
%% Function 2: Preprocessing Data, LPF-> Downsample -> HPF
function[PreProcessedData] = PreProcess(Channel, UserInputs)
Fs_original = UserInputs.PreProcessing(1);
Fs_new = UserInputs.PreProcessing(2);
Fc_low = UserInputs.PreProcessing(3);
Fc_High = UserInputs.PreProcessing(4);
GainMultiplier = UserInputs.PreProcessing(5);

%LPF; Second Order
Wn_low = Fc_low / (Fs_original / 2);  % Normalized cutoff frequency with respect to original Fs
order = 2;
[b_low, a_low] = butter(order, Wn_low, 'low');
ChannelLPF = filtfilt(b_low, a_low, Channel);

% Downsampling the Signal
downsample_factor = Fs_original / Fs_new;  % Factor by which to downsample
ChannelDownSample = downsample(ChannelLPF, downsample_factor);

% Highpass Filter
Wn_high = Fc_High / (Fs_new / 2);  % Normalized high-pass cutoff frequency
[b_high, a_high] = butter(order, Wn_high, 'high');  % 2nd-order Butterworth high-pass filter
ChannelOutput = (filtfilt(b_high, a_high, ChannelDownSample))*GainMultiplier;

% Creating Time Vector Based on Downsampled Data
timeVector_original = (0:1/Fs_new: ...
    (length(ChannelOutput)-1)/Fs_new)';

PreProcessedData.AllChannels = ChannelOutput;
PreProcessedData.Time = timeVector_original;
end
%% Function 1: Batch File Loader
function [FilePath, BaseNames] = FileLoader(myFolder)
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
        % User clicked Cancel
        return;
    end
end

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.txt'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
numFiles = length(theFiles);

FilePath = cell(numFiles, 1);
BaseNames = cell(numFiles, 1);

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = theFiles(k).folder;
    FilePath{k} = fullFileName;
    BaseNames{k} = baseFileName;
end
end
%% Things to do
% 1.) Granger Causality   X
% 2.) Adaptation from stimulation, fractional derivatives of signals at
% 0.25 fraction
% 3.) AM modulation, and look at phase shift of the derivative
% 4.) Changes in delay within each seizure event
% 5.) Check xcorr before a mirror might form
% 6.) Sample entropy for complexity
% 7.) Baseline analysis
% 8.) Seizure high frequency components non propagation
% 9.) Feature trends between seizures, regression model or machine learning
% 10.) For RMS plots, just pick a time point during baseline before
% injection and fixed time points post injection and plot that as an
% average for all samples
%% Function 0: Plotting Analysis Features
function [] = PlotData(DetectedEventFeatures, FilteredXCorr, UserInputs)



if UserInputs.Plots(4) == 1 %Cross Correlation Plots for Non Filtered Data
    bin_edges = -100:1:100;
    figure('name', 'Non Filtered', 'NumberTitle', 'off');
    subplot(2,3,1)
    histogram(DetectedEventFeatures.CrossCorrelation_Delays(:,1), bin_edges)
    title('Focus-CC(I)')

    subplot(2,3,2)
    histogram(DetectedEventFeatures.CrossCorrelation_Delays(:,2), bin_edges)
    title('Focus-CC(C)')

    subplot(2,3,3)
    histogram(DetectedEventFeatures.CrossCorrelation_Delays(:,3), bin_edges)
    title('Focus-Mirror')

    subplot(2,3,4)
    histogram(DetectedEventFeatures.CrossCorrelation_Delays(:,4), bin_edges)
    title('CC(I)-CC(C)')

    subplot(2,3,5)
    histogram(DetectedEventFeatures.CrossCorrelation_Delays(:,5), bin_edges)
    title('CC(I)-Mirror')

    subplot(2,3,6)
    histogram(DetectedEventFeatures.CrossCorrelation_Delays(:,6), bin_edges)
    title('CC(C)-Mirror')


    %Cross Correlation Plots for Filtered Data
    figure('name', 'Filtered', 'NumberTitle', 'off');
    subplot(2,3,1)
    histogram(FilteredXCorr.DetectedEvents(:,1), bin_edges)
    title('Focus-CC(I)')

    subplot(2,3,2)
    histogram(FilteredXCorr.DetectedEvents(:,2), bin_edges)
    title('Focus-CC(C)')

    subplot(2,3,3)
    histogram(FilteredXCorr.DetectedEvents(:,3), bin_edges)
    title('Focus-Mirror')

    subplot(2,3,4)
    histogram(FilteredXCorr.DetectedEvents(:,4), bin_edges)
    title('CC(I)-CC(C)')

    subplot(2,3,5)
    histogram(FilteredXCorr.DetectedEvents(:,5), bin_edges)
    title('CC(I)-Mirror')

    subplot(2,3,6)
    histogram(FilteredXCorr.DetectedEvents(:,6), bin_edges)
    title('CC(C)-Mirror')
else
end
end


%% Function 0: Plotting Analysis Features
function [] = PlotEvents(PreProcessedData, TeagerFeatures, TimeSeriesFeatures, DetectedEventFeatures, UserInputs)
% Checking which graphs to plot
PlotPreProcessed_Teager = UserInputs.Plots(1);
PlotTimeSeriesAnalysis = UserInputs.Plots(2);
PlotDetectedEventAnalysis = UserInputs.Plots(3);

% Preprocess and Teager Plot
if PlotPreProcessed_Teager == 1
    % Extract necessary data
    PreProcessedData_AllChannels = PreProcessedData.AllChannels;
    timeVector_original = PreProcessedData.Time;
    Teager = TeagerFeatures.TeagerAvg;
    Teager_Time = TeagerFeatures.TeagerTime;
    DetectedEvents = TeagerFeatures.DetectedEvents * 0.05; % Scale for visibility
    SLEThreshold = TeagerFeatures.SLEThresholdPlot;
    NoiseThreshold = TeagerFeatures.NoiseThresholdPlot;

    % Set axis limits based on 95% smallest values
    TeagerLimit = mean(maxk(abs(Teager), round(length(Teager) * 0.05)));
    ChannelLimit = mean(maxk(abs(PreProcessedData_AllChannels), round(size(PreProcessedData_AllChannels, 1) * 0.005)), 1);

    % Number of Channels and Time Window
    NumChannels = size(PreProcessedData_AllChannels, 2);
    timeWindow = 100;

    % Create figure for Preprocessed Data
    figure('Name', 'Preprocessed Data', 'NumberTitle', 'off');
    ax = gobjects(NumChannels + 1, 1);

    % Plot Preprocessed Data for Each Channel
    for ch = 1:NumChannels
        ax(ch) = subplot(NumChannels + 1, 1, ch);
        plot(ax(ch), timeVector_original, PreProcessedData_AllChannels(:, ch), 'b');
        hold(ax(ch), 'on');
        plot(ax(ch), Teager_Time, DetectedEvents, 'r'); % Detected events
        hold(ax(ch), 'off');
        SetCommonLabelsAndTitles(ax(ch), sprintf('Channel %d', ch), ChannelLimit(ch));
    end

    % Plot Teager Energy Features
    ax(NumChannels + 1) = subplot(NumChannels + 1, 1, NumChannels + 1);
    plot(ax(NumChannels + 1), Teager_Time, Teager, 'b', Teager_Time, SLEThreshold, 'g', Teager_Time, NoiseThreshold, 'r');
    SetCommonLabelsAndTitles(ax(NumChannels + 1), 'Event Detection Calibration', TeagerLimit); % Apply TeagerLimit
    legend('Teager Energy', 'Event Treshold', 'Noise Threshold');

    % Add super title and synchronize x-axes
    sgtitle('Preprocessed Data for All Channels');
    xlim(ax(1), [timeVector_original(1), timeVector_original(1) + timeWindow]);
    linkaxes(ax, 'x');

    % Slider for scrolling
    uicontrol('Style', 'slider', 'Min', 1, 'Max', timeVector_original(end) - timeWindow, ...
        'Value', 1, 'Units', 'normalized', 'Position', [0.1, 0.01, 0.8, 0.05], ...
        'Callback', @(src, ~) updateXLim(src, ax, timeWindow));
end

% Time Series Analysis Plots
if PlotTimeSeriesAnalysis == 1
    Time_TimeSeriesAnalysis = TimeSeriesFeatures.FeatureTime;
    RMS_TimeSeriesAnalysis = cell2mat(struct2cell(TimeSeriesFeatures.RMS)'); % Convert struct to matrix
    Teager_TimeSeriesAnalysis = cell2mat(struct2cell(TimeSeriesFeatures.Teager)');
    Entropy_TimeSeriesAnalysis = cell2mat(struct2cell(TimeSeriesFeatures.Entropy)');
    LineLength_TimeSeriesAnalysis = cell2mat(struct2cell(TimeSeriesFeatures.LineLength)');
    FrequencyBandPower_TimeSeriesAnalysis = TimeSeriesFeatures.FrequencyBandPower;

    % Frequency Bands Setup
    ChannelNames = fieldnames(TimeSeriesFeatures.FrequencyBandPower);
    [~, NumFrequencies] = size(TimeSeriesFeatures.FrequencyBandPower.Ch1);
    NumChannels = numel(ChannelNames);
    FrequencyBandsAre = UserInputs.TimeSeriesAnalysis{2};
    colorMap = lines(NumFrequencies);

    % Create Figures and Layouts
    fig_RMS = figure('Name', 'Time Series Analysis: RMS', 'NumberTitle', 'off');
    tl_RMS = tiledlayout(fig_RMS, NumChannels, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    fig_Teager = figure('Name', 'Time Series Analysis: Teager', 'NumberTitle', 'off');
    tl_Teager = tiledlayout(fig_Teager, NumChannels, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    fig_Entropy = figure('Name', 'Time Series Analysis: Shannon Entropy', 'NumberTitle', 'off');
    tl_Entropy = tiledlayout(fig_Entropy, NumChannels, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    fig_LineLength = figure('Name', 'Time Series Analysis: LineLength', 'NumberTitle', 'off');
    tl_LineLength = tiledlayout(fig_LineLength, NumChannels, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    fig_FreqBand = figure('Name', 'Time Series Analysis: Frequency Band Power', 'NumberTitle', 'off');
    tl_FreqBand = tiledlayout(fig_FreqBand, NumChannels, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Plot Data for Each Channel
    for ch = 1:NumChannels
        PlotChannelData(Time_TimeSeriesAnalysis, RMS_TimeSeriesAnalysis, Teager_TimeSeriesAnalysis, ...
            Entropy_TimeSeriesAnalysis, LineLength_TimeSeriesAnalysis, FrequencyBandPower_TimeSeriesAnalysis, ...
            NumFrequencies, FrequencyBandsAre, ChannelNames, colorMap, ch, tl_RMS, tl_Teager, tl_Entropy, tl_LineLength, tl_FreqBand);
    end
end

% Detected Event Plots
if PlotDetectedEventAnalysis == 1
    Time_DetectedEvents = DetectedEventFeatures.Time;
    CrossCorrelationRMax_DetectedEvents = DetectedEventFeatures.CrossCorrelation_RMax;
    CrossCorrelationDelays_DetectedEvents = DetectedEventFeatures.CrossCorrelation_Delays;
    CrossCorrelation_ChannelCombinations = DetectedEventFeatures.ChannelCombinations;
    LineLength_DetectedEvents = DetectedEventFeatures.LineLength;
    Entropy_DetectedEvents = DetectedEventFeatures.Entropy;
    RMS_DetectedEvents = DetectedEventFeatures.RMS;
    EventDuration_DetectedEvents = DetectedEventFeatures.EventDuration;

    % Cross-Correlation Plots
    PlotCrossCorrelationHistograms(CrossCorrelationRMax_DetectedEvents, CrossCorrelationDelays_DetectedEvents, ...
        CrossCorrelation_ChannelCombinations);

    % Line Graphs for Event Features
    figure('Name', 'Detected Event Analysis', 'NumberTitle', 'off');
    tiledLayout = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    PlotDetectedEventFeature(Time_DetectedEvents, LineLength_DetectedEvents, 'Line Length', tiledLayout);
    PlotDetectedEventFeature(Time_DetectedEvents, Entropy_DetectedEvents, 'Shannon Entropy', tiledLayout);
    PlotDetectedEventFeature(Time_DetectedEvents, RMS_DetectedEvents, 'RMS', tiledLayout);
    PlotDetectedEventFeature(Time_DetectedEvents, EventDuration_DetectedEvents, 'Event Duration', tiledLayout);
end
end
% Helper Functions for Plots
function PlotDetectedEventFeature(Time, FeatureData, FeatureName, tiledLayout)% Plots Detected Event Features
% Function to plot individual event features
nexttile(tiledLayout);
plot(Time, FeatureData, 'LineWidth', 1.5);
title(FeatureName);
xlabel('Time (s)');
ylabel(FeatureName);
grid on;
end
function PlotCrossCorrelationHistograms( ...%Plots Histograms for RMax and Delays
    CrossCorrelationRMax_DetectedEvents, CrossCorrelationDelays_DetectedEvents, CrossCorrelation_ChannelCombinations)
% Get the number of channel combinations
numCombinations = size(CrossCorrelationRMax_DetectedEvents, 2);

% Calculate the number of rows and columns for subplots
numRows = ceil(sqrt(numCombinations));  % Get the number of rows (rounded up)
numCols = ceil(numCombinations / numRows);  % Get the number of columns

% RMax Plots
figure('Name', 'Cross-Correlation RMax Analysis', 'NumberTitle', 'off');
for i = 1:numCombinations
    % Create subplot for RMax
    subplot(numRows, numCols, i);

    % Plot RMax histogram
    histogram(CrossCorrelationRMax_DetectedEvents(:, i), 20, 'FaceColor', 'b', 'EdgeColor', 'b', 'FaceAlpha', 0.7);

    % Get channel pair for the title
    channelPair = CrossCorrelation_ChannelCombinations(i, :);

    % Set title, labels, and grid
    title(sprintf('RMax Histogram for Channels %d & %d', channelPair(1), channelPair(2)));
    xlabel('RMax Value');
    ylabel('Frequency');
    grid on;
end

% Delay Plots
figure('Name', 'Cross-Correlation Delay Analysis', 'NumberTitle', 'off');
binEdges = -100:5:100;
for i = 1:numCombinations
    % Create subplot for Delay
    subplot(numRows, numCols, i);

    % Plot Delay histogram
    histogram(CrossCorrelationDelays_DetectedEvents(:, i), binEdges, 'FaceColor', 'r', 'EdgeColor', 'r', 'FaceAlpha', 0.7);

    % Get channel pair for the title
    channelPair = CrossCorrelation_ChannelCombinations(i, :);

    % Set title, labels, and grid
    title(sprintf('Delay Histogram for Channels %d & %d', channelPair(1), channelPair(2)));
    xlabel('Delay (ms)');
    ylabel('Frequency');
    grid on;
end
annotation('textbox', [0.5, 0.01, 0.4, 0.05], 'String', 'Footnote: Negative Delay Signifies the channel mentioned first is leading.', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
function PlotChannelData(Time, RMS, Teager, Entropy, LineLength, FrequencyBandPowers, ... % Plots Time Series Analysis
    NumFrequencies, FreqBandLimits, ChannelNames, colorMap, ch, ...
    tl_RMS, tl_Teager, tl_Entropy, tl_LineLength, tl_FreqBand)

% RMS subplot
nexttile(tl_RMS);
plot(Time, RMS(:, ch), 'r');
title(sprintf('Channel %d', ch));
ylabel('RMS');
grid on;

% Teager subplot
nexttile(tl_Teager);
plot(Time, Teager(:, ch), 'b');
title(sprintf('Channel %d', ch));
ylabel('Teager');
grid on;

% Shannon Entropy subplot
nexttile(tl_Entropy);
plot(Time, Entropy(:, ch), 'k');
title(sprintf('Channel %d', ch));
ylabel('Shannon Entropy');
grid on;

% Line Length subplot
nexttile(tl_LineLength);
plot(Time, LineLength(:, ch), 'r');
title(sprintf('Channel %d', ch));
ylabel('Line Length');
grid on;

% Frequency Band Power subplot
nexttile(tl_FreqBand);
for freq = 1:NumFrequencies %for each frequency band, we are already in another for loop for each channel indexed with ch
    FrequencyBandPower_IndividualChannel =  FrequencyBandPowers.(ChannelNames{ch});
    plot(Time, FrequencyBandPower_IndividualChannel(:, freq), 'Color', colorMap(freq, :));
    hold on
end
hold off;
title(sprintf('Channel %d', ch));
ylabel('Frequency Band Power');
grid on;

legendLabels = arrayfun(@(low, high) sprintf('%.1f - %.1f Hz', low, high), ...
    FreqBandLimits(:, 1), FreqBandLimits(:, 2), 'UniformOutput', false);
legend(legendLabels, 'Location', 'best');
end
function SetCommonLabelsAndTitles(ax, titleText, yLimit) %Sets Limits and Labels for Preprocessing and Teager Plots
xlabel(ax, 'Time (s)');
title(ax, titleText);
% Set the lower bound to 0 and upper bound to TeagerLimit for Teager plot only
if strcmp(titleText, 'Event Detection Calibration')
    ylim(ax, [0, yLimit]);  % Set lower bound to 0 and upper bound to TeagerLimit
    ylabel(ax, 'Teager Energy');
else
    ylim(ax, [-yLimit, yLimit]);  % Apply ChannelLimit for other plots
    ylabel(ax, 'Amplitude (mV)');
end
grid(ax, 'on');

end
function updateXLim(src, ax, timeWindow) %Scroll Bar for PreProcessing and Teager Plots
new_xlim = [src.Value, src.Value + timeWindow];
set(ax, 'XLim', new_xlim);
end
