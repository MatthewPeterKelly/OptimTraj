function NewOpt = mergeDefaults(Opt, Default)
% Opt = mergeDefaults(Opt, Default)
%
% This function merges the stuct Opt and the struct Default. If a field is
% present in both structs, the Opt takes precedence. If a field is found in
% Default but not Opt, then it is added. If a field is found in Opt but not
% in Default, then a warning is printed.
%

% Set default, and override if found in Opt
NewOpt = mergeFields(Opt, Default);

% Check if there were any unused fields in Opt
checkFields(Opt, NewOpt);

end


function NewOpt = mergeFields(Opt, Default)
%
% This function merges the stuct Opt and the struct Default. If a field is
% present in both structs, the Opt takes precedence. If a field is found in
% Default but not Opt, then it is added.
%

names = fieldnames(Default);
for i=1:length(names)
    if isstruct(Default.(names{i}))
        if ~isfield(Opt,names{i})
            Opt.(names{i}) = [];
        end
        NewOpt.(names{i}) = mergeFields(Opt.(names{i}), Default.(names{i}));
    else
        if isfield(Opt,names{i})
            NewOpt.(names{i}) = Opt.(names{i});
        else
            NewOpt.(names{i}) = Default.(names{i});
        end
    end
end

end



function checkFields(Opt, NewOpt)
%
% This function checks for any fields that are present in Opt, but not in
% NewOpt, and prints a warning if so.
%

if isstruct(Opt)
    names = fieldnames(Opt);
    for i=1:length(names)
        if isstruct(Opt.(names{i}))
            checkFields(Opt.(names{i}), NewOpt.(names{i}));
        else
            if ~isfield(NewOpt,names{i})
                disp(['Warning:  field ''' names{i} ''' in options struct is unused']);
            end
        end
    end
end
end