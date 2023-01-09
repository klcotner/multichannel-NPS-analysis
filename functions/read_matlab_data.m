%{
[rawtime, rawdata, samplerate] = read_matlab_data(path_to_datafile, sensitivity, [index_range], [new_filepath])

INPUTS
* path_to_datafile = char pointing to the acquired .mat file
    > the file must contain "data" (2xn [seconds; voltage]) and "fs" (sample rate in Hz)
* sensitivity = preamp sensitivity during acquisition (A/V)
* index_range = (optional, array or bool) = [start_ix, end_ix]
    > if array, will only read in the column indices start_ix:end_ix (inclusive)
    > if false, will include all the data
    > defaults to false
* new_filepath = (optional, char or bool)
    > if char, saves to new_filepath (must end in '.mat')
    > if true, saves to <inferred_stem>_input.mat in the same directory as the data
    > if false, doesn't save
    > defaults to false

OUTPUTS
* rawtime = nx1 column vector of doubles (raw time values in seconds, starting from 0)
* rawdata = nx1 column vector of doubles (raw current readings in amps)
* samplerate = acquisition sample rate, in Hz

NOTES
>> note: checks that fs matches the inferred sample rate from the data vector
>> note: '../Util' must be on the search path (for check_filepath_mat, check_bool, and unused_filename)
%}

function [rawtime, rawdata, samplerate] = read_matlab_data(path_to_datafile, sensitivity, index_range, new_filepath)

    %% input parser
    p = inputParser;

    p.addRequired('path_to_datafile', @check_filepath_mat);
    p.addRequired('sensitivity', @(x) isnumeric(x) & isscalar(x) & x>0);
    p.addOptional('index_range', false, @check_index_range );
    p.addOptional('new_filepath', false, @(x) check_filepath_mat(x) || check_bool(x) );

    p.StructExpand = false;

    if nargin==2
        parse(p, path_to_datafile, sensitivity);
        index_range = p.Results.index_range;
        new_filepath = p.Results.new_filepath;
    elseif nargin==3
        parse(p, path_to_datafile, sensitivity, index_range);
        new_filepath = p.Results.new_filepath;
    else
        parse(p, path_to_datafile, sensitivity, index_range, new_filepath);
    end

    %% load data from .mat file and check it
    
    disp('reading raw data');
    
    m = matfile(fullfile(path_to_datafile));
    data = m.data;
    fs = m.fs;

    % check `data`
    if ~isnumeric(data)
        error('the variable `data` in the input file must be numeric');
    end
    if size(data,1) ~= 2
        error('the variable `data` in the input file must have size 2xn');
    end

    % check `fs`
    if ~isnumeric(fs)
        error('the variable `fs` in the input file must be numeric');
    end
    if ~isscalar(fs)
        error('the variable `fs` in the input file must be scalar');
    end
    if fs<=0
        error('the variable `fs` in the input file must be positive');
    end

    % check whether `fs` matches the inferred sample rate from `data`
    timediff = diff(data(1,:));
    interval = 1/fs;
    tolerance = (interval) * 1e-5;
    if any( abs(timediff - interval) > tolerance )
        error('data time points must be spaced at intervals equal to 1/fs');
    end

    %% parse input data

    % include only the columns from index_range
    if index_range
        if any(index_range > size(data,2)) % check that values are valid
            error('`index_range` exceeds size of `data`');
        else
            data = data(:, index_range(1):index_range(2) );
        end
    end

    % split into column vectors for time & current
    rawtime = data(1,:)'; % sec
    rawdata = data(2,:)' * sensitivity; % A
    
    % ensure that rawtime starts at 0sec
    if rawtime(1) ~= 0
        numpoints = numel(rawtime);
        endtime = (1/fs) * (numpoints-1);
        rawtime = (0:(1/fs):endtime)';
    end

    samplerate = fs;    
    disp('finished reading raw data');

    %% save parsed data, if desired

    if check_bool(new_filepath)
        dosave = new_filepath;
        [folder, filename, ~] = fileparts(path_to_datafile);
        savename = fullfile(folder, [filename, '_input.mat']);
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

%% helper function to check for valid values for `index_range`
function tf = check_index_range(x)
    % index_range can be an array with 2 elements (positive & increasing)
    if isnumeric(x) && numel(x)==2
        tf = all(x>0) && x(2)>x(1);
    % index_range can have the value of false
    elseif check_bool(x)
        tf = ~x;
    % otherwise, x is an invalid value for index_range
    else
        tf = false;
    end
end
