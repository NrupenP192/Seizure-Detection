clc
clear 
close all

%% Loading file data
myFolder = 'C:\Users\Nrupen Pakalapati\OneDrive\Documents\MATLAB\During LFS';
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
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = theFiles(k).folder;
    slashString = '\';
    FilePath{k} = fullFileName;
    BaseNames{k} = baseFileName;
end

clearvars -except FilePath BaseNames
savedir = 'C:\Users\nrupe\Desktop\Matlab Files\In-Vivo Chronic\Matlab Code for analysis\Intan Files\Rat 1 Test';

%% Loading Data
for jj = 5: 5%length(BaseNames)
    filename = fullfile(FilePath(jj), BaseNames(jj));
    data = readmatrix(filename{1}); %reads two channels of data, double check which channel is CC and which is cortex
    CorpusCallossum = data(:,3);

    % User Inputs:
    Fs_original = 40000;  % Original sampling frequency (20,000 Hz)
    Fs_new = 10000;        % New sampling frequency (1,000 Hz) for downsampling
    Fc_low = 40;         % Lowpass filter cutoff frequency
    Fc_High = 0.25;          % Highpass filter cutoff frequency
    PadLength = 100;      % Padding either side of the teager window (in samples) with data.
    GainMultiplier = 10;   % Gain = 100; want things in mV, so 1000/100 = 10
    plotFigures = 0;        % 1 for plots, 0 for no plots
    %Teager Inputs
    window_size = 100; % 1000hz = 1s window
    NoiseMultiplier = 100000000;
    SLEMultiplier = 500;

    %% Preprocessing Data: Highpass Filter, Low Pass Filter, Downsample
    t_old = (0:1/Fs_original:(length(CorpusCallossum)-1)/Fs_original)'; %time vector
    % Design a 2nd order Butterworth low-pass filter and apply it
    Wn_low = Fc_low / (Fs_original / 2);  % Normalized cutoff frequency with respect to original Fs
    order = 2;
    % [b_low, a_low] = butter(order, Wn_low, 'low');
    % CorpusCallossum = filtfilt(b_low, a_low, CorpusCallossum);

    % Downsample the filtered signal
    downsample_factor = Fs_original / Fs_new;  % Factor by which to downsample
    CorpusCallossum = downsample(CorpusCallossum, downsample_factor);
    t = (0:1/Fs_new:(length(CorpusCallossum)-1)/Fs_new)'; %time vector

    % Highpass Filter
    Wn_high = Fc_High / (Fs_new / 2);  % Normalized high-pass cutoff frequency
    [b_high, a_high] = butter(order, Wn_high, 'high');  % 2nd-order Butterworth high-pass filter
    CorpusCallossum = (filtfilt(b_high, a_high, CorpusCallossum))*GainMultiplier;


    %% Teager Tresholding for stimulation artifacts
    % Define the window size and overlap for moving windows
    window_overlap = round(window_size/1.1111); % 90% overlap

    % Compute the windows
    windows = buffer((CorpusCallossum), window_size, window_overlap, 'nodelay');
    windowsT = buffer(t, window_size, window_overlap, 'nodelay');

    % Compute the Teager non-linear energy values
    %Teager Modifiers
    Tn = 2;
    %Teager discrete formula where n can be changed, and is vectorized: y = x(n+1:end-n).^2 - x(1:end-2n) . x(2*n+1:end); %general formula
    teager_values = (sum(windows(Tn+1:end-Tn, :).^2 - windows(1:end-2*Tn,:).*windows(2*Tn+1:end,:),1));
    teager_avg = teager_values./window_size;


    %teager_values = (sum(windows(2:end-1, :).^2 - windows(1:end-2, :).*windows(3:end, :), 1));%x^2 - x(n-1)*x(n+1), %time step = 0.1s, %f = 10hz
    time_axis = windowsT((window_size/2),:); %taking the midpoint of each window as the index for the teager value
    teager_avg = smoothdata(teager_avg,'movmean','SmoothingFactor',0.5,'SamplePoints',time_axis);


    % Calculating the Step Function for seizure detection Based on the Teager values
    minteager = mink(abs(teager_avg), round(length(teager_avg)*0.1)); % Takes the minimum 10% of teager values in the dataset and uses them to threshold
    avg_min = mean(minteager);
    maxteager = maxk(abs(teager_avg), round(length(teager_avg)*0.8));
    % SLEThreshold = avg_min*SLEMultiplier;
    SLEThreshold = mean(maxteager);
    NoiseThreshold = avg_min*NoiseMultiplier;

    %Threshold Plot
    indexSLE = find(teager_avg >= SLEThreshold & teager_avg <= NoiseThreshold);
    StepEnvelope = zeros(1, length(teager_avg));

    StepFunctionAmplitude = 5;
    StepEnvelope(indexSLE) = StepFunctionAmplitude;
    SLEStartCode = [0, StepFunctionAmplitude];
    SLEEndCode = [StepFunctionAmplitude, 0];

    %Creating Step Function for Event Marking
    [EVPStartTime, EVPEndTime, EVPStartIndex, EVPEndIndex] = findStartEndPoints(StepEnvelope, time_axis, SLEStartCode, SLEEndCode);
    Start = arrayfun(@(i) find(t == EVPStartTime(i)), 1:length(EVPStartTime));
    End = arrayfun(@(i) find(t == EVPEndTime(i)), 1:length(EVPEndTime));
    timeDiff = End - Start;
    minTime = min(timeDiff);
    halfTime = (minTime/3); %for centering everything around stim artifacts, assuming stim artifacts are similar
    
    %% Identifying evoked potentials in the window
    colors = lines((length(Start)/3)); % Get 5 distinct colors (you can modify this)
    num_colors = size(colors, 1); % Number of unique colors available
    initial_alpha = 1; % Fully opaque
    alpha_step = 1/length(Start); % How much alpha reduces per step

    figure(2)
    for c = 1: length(Start)
        EVPSnippet = CorpusCallossum(Start(c):End(c));
        %tSnippet = t(Start(c):End(c));
        
        [pks_artifact, locs_artifact] = findpeaks(smoothdata(-EVPSnippet,'movmean','SmoothingFactor',0.1), 'MinPeakProminence', 1);
        
        StimArtifact = locs_artifact(find(max(pks_artifact)));

        EVPSnippet = EVPSnippet(StimArtifact-halfTime : StimArtifact+halfTime);
        tSnippet = (0:1/Fs_new:(length(EVPSnippet)-1)/Fs_new)'; %time vector

        alpha_value = max(initial_alpha - (c-1) * alpha_step, 0.3); % Ensure alpha doesn't go below 0.1
        color_idx = mod(floor((c-1)/6),num_colors)+1;
        faded_color = (1 - alpha_value) * [1, 1, 1] + alpha_value * colors(color_idx, :);
        plot(tSnippet, EVPSnippet, 'k', 'linewidth', 1, 'color', faded_color);
        hold on
        %plot(tSnippet(locs), -pks, 'ro')
        ylim([-1.5 0.5]);
        dumdum = 1;
        
    end

    if plotFigures == 1;
        figure(jj);
        SLEThresholdPlot = SLEThreshold + zeros(1,length(t));
        NoiseThresholdPlot = NoiseThreshold + zeros(1,length(t));
        numPlots = 5;
        axesArray = gobjects(1, numPlots);

        % Generate each subplot individually with sample data
        xData = linspace(0, 100, 1000); % Define x-data range here
        ax1 = subplot(4, 1, 1);
        plot(ax1, t, CorpusCallossum);
        hold on
        plot (ax1, time_axis, StepEnvelope*0.1)
        ylim([-1 0.2]);


        ax2 = subplot(4, 1, 2);
        plot(ax2, time_axis, teager_avg);
        hold on
        plot(ax2, t, SLEThresholdPlot)
        hold on
        plot(ax2, t, NoiseThresholdPlot)
        ylim([0 0.01])

        ax3 = subplot(4, 1, 3);

        ax4 = subplot(4, 1, 4);

        % Link x-axes of all subplots
        linkaxes([ax1, ax2, ax3, ax4], 'x');

        % Set initial X-axis limits (adjust as needed)
        initialXLim = [0, 10];
        set([ax1, ax2, ax3, ax4], 'XLim', initialXLim);

        % Calculate the maximum scroll position based on xData range and initialXLim width
        xRange = t(end) - t(1);
        scrollMax = xRange - diff(initialXLim);

        % Create a slider to control x-axis scrolling
        scrollSlider = uicontrol('Style', 'slider', 'Min', xData(1), 'Max', scrollMax, ...
            'Value', initialXLim(1), 'Units', 'normalized', ...
            'Position', [0.1, 0.01, 0.8, 0.05]);

        % Callback function to update x-axis limits based on slider position
        scrollSlider.Callback = @(src, event) set([ax1, ax2, ax3, ax4], 'XLim', [src.Value, src.Value + diff(initialXLim)]);
    end


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