function ops = gitStatus(ops)
% function ops = gitStatus(ops)
% 
% Retrieve & store info about current state of git repo(s)
% - enable for use by setting flag:
%   [ops.useGit] == true;
% - if enabled, will track kilosort git repo by default
% - see code comments for how to add your own repos for tracking
%
% 2021-03-08  T.Czuba  Wrote it to inherit git tracking aspects from PLDAPS
% 

if getOr(ops, 'useGit', 0)
    if ~isfield(ops,'git')
        ops.git.kilosort = struct('status',[]);
    end
    
    fn = fieldnames(ops.git);
    fn = ['kilosort', fn(~strcmp(fn{:},'kilosort'))];
    % Each fieldname under [ops.git] should correspond to a separate git repo to be tracked
    % - [.kilosort]  will be tracked by default
    % - follow format to add additional repos (e.g. lab-specific kilosort config, data preprocessing/conversion, utilities, etc)

    for i = 1:length(fn)
        try
            % only execute git check once
            if ~isfield(ops.git.(fn{i}),'status') || isempty(ops.git.(fn{i}).status)
                % Retrieve info about repo path, source, status, version/commit, & diff
                if ~isfield(ops.git.(fn{i}), 'mainFxn') || isempty(ops.git.(fn{i}).mainFxn)
                    % if no main function provided, assume repo name is a unique function in the repo's base directory
                    % - e.g. 'kilosort.m'
                    ops.git.(fn{i}).mainFxn = [fn{i},'.m'];
                end
                ops.git.(fn{i}).basePath = fileparts(which(ops.git.(fn{i}).mainFxn));
                thisBase = ops.git.(fn{i}).basePath; % shorthand
                
                ops.git.(fn{i}).status      = git.git(['-C ' thisBase ' status']);
                ops.git.(fn{i}).remote      = git.git(['-C ' thisBase ' remote -v']);
                ops.git.(fn{i}).branch      = git.git(['-C ' thisBase ' symbolic-ref --short HEAD']);
                ops.git.(fn{i}).revision    = git.git(['-C ' thisBase ' rev-parse HEAD']);
                ops.git.(fn{i}).diff        = git.git(['-C ' thisBase ' diff']);
            end
        catch
            % separate error field allows everything up to error to carry through
            errString = sprintf('%s repo not found, or inaccessible.',fn{i})
            ops.git.(fn{i}).error      = errString;
        end
    end
end
