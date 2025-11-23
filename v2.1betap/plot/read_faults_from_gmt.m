function faults = read_faults_from_gmt(inputFileName)

inputFileID = fopen(inputFileName, 'r');
textLine = fgetl(inputFileID);

faults=[];

while ischar(textLine)
    if strcmpi(textLine(1),'>')
        if size(faults,1) ~= 1
            faults(size(faults,1)+1,:) = [nan,nan];
        end
    elseif ~strcmpi(textLine(1),'#')
        ll = split(textLine);
        if size(faults,1) == 0
            faults = [str2double(ll{1}), str2double(ll{2})];
        else
        faults(size(faults,1) + 1,:) = [str2double(ll{1}), str2double(ll{2})];
        end
    end
    textLine = fgetl(inputFileID);
end

fclose(inputFileID);