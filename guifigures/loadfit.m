function fit = loadfit(fn)
%FIT = LOADFIT(FN) loads the file FN and returns the variable named fit or
%   fitg. Other variables in the file are assigned in the caller workspace.
% Jacob Kowalewski 2010 SPTrack v4

loadstruct = load(fn,'-mat');
vars = fieldnames(loadstruct);
for i = 1:length(vars)
    switch vars{i}
        case 'fitg'
            fit = loadstruct.fitg;
        case 'fit'
            fit = loadstruct.fit;
        otherwise
            assignin('caller',vars{i},loadstruct.(vars{i}));
    end
end
end


