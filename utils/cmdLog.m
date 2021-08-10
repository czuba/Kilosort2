function cmdLog(inStr, varargin)
% 
% Standardize Kilosort command line logging
% 

timeStr = datestr(now, 'HH:MM:SS');
if ~isempty(varargin)
    tocSec = varargin{1};
    if tocSec<120
        timeStr = sprintf('%s (%2.1fsec)', timeStr, tocSec);
    else
        timeStr = sprintf('%s (%2.1fmin)', timeStr, tocSec/60);
    end
end

fprintf('  %s:\t%s\n', timeStr, inStr);