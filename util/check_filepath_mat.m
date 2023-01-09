% helper function to check whether x is a valid path to a .mat file
function tf = check_filepath_mat(x)
    % x must be a char/str ending in .mat, with a filename & valid folder
    if ischar(x) || isstring(x)
        [folder, filename, ext] = fileparts(x);
        tf = strcmpi(ext,'.mat') && ~isempty(filename) && (isfolder(folder) || isempty(folder));
    % if x isn't a char/str, it isn't a valid filepath
    else
        tf = false;
    end
end