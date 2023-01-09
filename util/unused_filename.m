% if oldname exists, append a number to the base filename
% newname = the potentially altered filename
% nameflag = whether or not the filename had to be altered
% NB: preserves filepath as well

function [newname, nameflag] = unused_filename(oldname)

    [filepath, basename, ext] = fileparts(oldname);
    
    iter = 0;
    while isfile(oldname)
        iter = iter + 1;
        oldname = fullfile(filepath, [basename, '_', num2str(iter), ext]);
    end
    
    newname = fullfile(oldname);
    nameflag = iter>0;
    
end
