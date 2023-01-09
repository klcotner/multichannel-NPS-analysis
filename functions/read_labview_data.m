%{
[rawtime, rawdata, samplerate] = read_labview_data(path_to_directory, sample_rate, sensitivity, [file_numbers], [new_filepath])

INPUTS
* path_to_directory = char pointing to the folder containing the acquired .txt files
    > directory must contain .txt files
    > data files must contain one data point on each line
    > file names ending in '-notes' will be ignored
* sample_rate = acquisition sample rate, in Hz
* sensitivity = preamp sensitivity during acquisition (A/V)
* file_numbers = (optional, array or bool) = [start_filenum, end_filenum]
    > if array, will only read in the file numbers start_filenum:end_filenum (inclusive)
        > this assumes all data file names end in '_' followed by a 5-digit number
        > this will not throw an error even if not all of the file numbers are valid
    > if false, will include all the data files
    > defaults to false
* new_filepath = (optional, bool or char)
    > if char, saves to new_filepath (must end in '.mat')
    > if true, saves to <inferred_folder>/<inferred_name>_input.mat
        > inferred_name is the name of the folder that contains the data files
        > inferred_folder is the parent folder of path_to_directory
    > if false, doesn't save
    > defaults to false

OUTPUTS
* rawtime = nx1 column vector of doubles (raw time values in seconds, starting from 0)
* rawdata = nx1 column vector of doubles (raw current readings in amps)
* samplerate = acquisition sample rate, in Hz

NOTES
>> note: ignores non-text files and files ending in '-notes.txt'
>> note: '../Util' must be on the search path (for check_filepath_mat, check_bool, and unused_filename)
%}

function [rawtime, rawdata, samplerate] = read_labview_data(path_to_directory, sample_rate, sensitivity, file_numbers, new_filepath)
    
    %% input parser
    p = inputParser;

    p.addRequired('path_to_directory', @(x) ischar(x) && isfolder(x));
    p.addRequired('sample_rate', @(x) isnumeric(x) && isscalar(x) && x>0);
    p.addRequired('sensitivity', @(x) isnumeric(x) && isscalar(x) && x>0);
    p.addOptional('file_numbers', false, @check_file_nums);
    p.addOptional('new_filepath', false, @(x) check_filepath_mat(x) || check_bool(x) );

    p.StructExpand = false;

    if nargin==3
        parse(p, path_to_directory, sample_rate, sensitivity);
        file_numbers = p.Results.file_numbers;
        new_filepath = p.Results.new_filepath;
    elseif nargin==4
        parse(p, path_to_directory, sample_rate, sensitivity, file_numbers);
        new_filepath = p.Results.new_filepath;
    else
        parse(p, path_to_directory, sample_rate, sensitivity, file_numbers, new_filepath);
    end
    
    %% parse files in the data directory
    
    disp('reading raw data');
    
    % get all .txt files in the directory
    infiles = dir(fullfile(path_to_directory, '*.txt'));
    if isempty(infiles)
        error('path_to_directory does not contain any .txt files');
    end
    
    % eliminate file names ending in '-notes'
    notesfiles = endsWith({infiles.name}, '-notes.txt');
    infiles = infiles(~notesfiles);
    if isempty(infiles)
        error('path_to_directory does not contain any valid data files');
    end
    
    % include only the desired range of file_numbers
    if file_numbers
        filenum_range = compose('_%05d.txt', file_numbers(1):file_numbers(2));
        inrangefiles = endsWith({infiles.name}, filenum_range);
        infiles = infiles(inrangefiles);
        if isempty(infiles)
            error('path_to_directory does not contain any data files in the range specified by file_numbers');
        end
    end
    
    %% read data & create output variables
    
    % concatenate data from all files into one column vector
    data_V = [];
    for ii=1:numel(infiles)
        filepath = fullfile(infiles(ii).folder, infiles(ii).name);
        filedata = importdata(filepath, ',', 0);
        data_V = [data_V; filedata];
    end
    
    % convert voltage to current
    rawdata = data_V * sensitivity;
    
    % create times vector
    numpoints = numel(rawdata);
    endtime = (1/sample_rate) * (numpoints-1);
    rawtime = (0:(1/sample_rate):endtime)';
    
    % finish
    samplerate = sample_rate; % different name for output
    disp('finished reading raw data');
    
    %% save parsed data, if desired

    if check_bool(new_filepath)
        dosave = new_filepath;
        
        % infer folder & filename for saving
        if ~endsWith(path_to_directory, filesep)
            path_to_directory = [path_to_directory, filesep];
        end
        pathsplit = split(path_to_directory, filesep);
        filestem = pathsplit{end-1};
        foldername = join(pathsplit(1:end-2), filesep);
        foldername = foldername{1};
        savename = fullfile(foldername, [filestem, '_input.mat']);
        
    else
        dosave = true;
        savename = new_filepath;
    end

    if dosave
        disp('saving raw data');
        [newsavename, saveflag] = unused_filename(savename);
        if saveflag
            warning([savename, ' already exists. suffix appended to prevent overwrite']);
        end
        
        save(newsavename, 'rawtime', 'rawdata', 'samplerate');
        disp(['raw data saved to ', newsavename]);
    end
    
end

%% helper function to check for valid values for `file_numbers`
function tf = check_file_nums(x)
    % file_numbers can be an array with 2 elements (positive & non-decreasing)
    if isnumeric(x) && numel(x)==2
        tf = all(x>0) && x(2)>=x(1);
    % file_numbers can have the value of false
    elseif check_bool(x)
        tf = ~x;
    % otherwise, x is an invalid value for file_numbers
    else
        tf = false;
    end
end

