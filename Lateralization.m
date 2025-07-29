clc
clear
close all
data = load("Lateralization.txt");
class = load("Class.txt");

for i = 1:numel(class)
    currentValue = class(i);

    % Check the value and assign the corresponding string.
    if currentValue == 0
        outputCellArray{i} = 'Synchronized';
    elseif currentValue == 1
        outputCellArray{i} = 'Focus';
    elseif currentValue == 2
        outputCellArray{i} = 'Mirror';
    else
        % Handle unexpected values (optional: could throw an error or assign 'Unknown')
        warning('Unexpected value encountered: %d. Assigning "Unknown".', currentValue);
        outputCellArray{i} = 'Unknown';
    end
end

class = outputCellArray';

FocustoCCI = data(:,1);
FocustoCCC = data(:,2);
CCItoCCC = data(:,3);
CCItoMirror = data(:,4);
CCCtoMirror = data(:,5);

data2 = [FocustoCCI, FocustoCCC, CCItoMirror, CCCtoMirror];