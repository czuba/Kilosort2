% ksGUI_updatePars.m
% 
% script to update GUI pars with recommend settings for stereo Uprobes (32 chan, 50/100 um geometry)
% 

fprintf(['\n',repmat('=-',1,20),'\n']);

%% Retrieve gui settings
ks = get(figure(1029321), 'UserData'); % use standard kilosort [magic] figure number to select proper handle
ops = ks.ops
ops0 = ops; % backup initial ops settings

fprintf([repmat('--',1,20),'\n']);


%% Custom fields or flags
% Parallel memmap raw data for whitening/loading
% - requires use of get_whitening_matrix_faster.m w/in preprocessDataSub.m
ops.useMemMapping   = 1;

ops.fig = 2; % 1==standard plots, 2=extra debug plots (more verbose, but marginally slower)

% time range for spike sorting analysis
% - def = [0 inf];
% ops.trange = [0, inf];

% drift correction params (see datashift2.m)
%       Probe recommendation:  ops.nblocks = 1;  ops.integerShifts = 1;
% [.nblocks] type of data shifting (0 = none, 1 = rigid, 2 = nonrigid)
ops.nblocks = 0; % non-rigid only really relevant for mmmany channels or probe length is long relative to brain (i.e. rodents)

% flag to round [data]shifts to nearest electrode spacing integer
% - ALWAYS use integerShifts
ops.integerShifts = 1;

% preselect target batch for drift alignment
% - if  <1, will be batch nearest targBatch% of total batches
% - if >=1, will be direct index to batch#
% - default = 2/3;
% ops.targBatch = 0.3;

%%



% Randomize batch order during learning
% - provides more stable/effective set of learned templates across entire file
ops.learnRand = 1;

% clip template updating to a minimum number of contributing spikes
% - helps prevent inversions (due to subtle/irregular noise being injected into dWU0 output of mexMPnu8.cu)
% 20 spike cutoff works well for 10 sec batch
ops.clipMin = 20;
ops.clipMinFit = .8;  % can survive clipping if median accounts for at least this much variance (ratio of vexp./amp)

% Apply detailed ccg analysis function to remove double-counted spike clusters
% - Orig from Bondy fork, but integrated into standard around kilosort 2(.5)
% - this is useful feature, but actually makes manual curation somewhat more challenging,
%   because strong ccg peak is informative for merge decisions
% - best left disabled for probes (even hopes that it would allow threshold cutoff to be less errorprone didn't work out)
ops.rmDuplicates    = 0;

% Post-hoc split clusters by:  1==template projections, 2==amplitudes, 0==don't split
% - amplitude splits seem reasonably trustworthy (...template splits suceptible to oddities of templates (e.g. inversions))
% ops.splitClustersBy = 2;    % (relatively safe & effective, but can mask problems w/sort parameters & fitting)
ops.splitClustersBy = 0;    % 0 recommended for full assessment of what sorting is doing

% standard cutoff can be overly aggressive
% - best left disabled for probes
ops.applyCutoff = 0;

% Git repo status & diff
ops.useGit = 1;

% add this kilosort_utils repo to git version tracking
% - moved into ks25 ./configFiles dir
% - retained here for example on how to add other repos to git tracking functionality
% % ops.git.kilosort_utils.mainFxn = 'ksGUI_updatePars.m';


%% Apply changes
ops.fshigh = 300; % map system has hardware high pass filters at 300

% make waveform length independent of sampling rate
% ops.nt0                 = ceil( 0.002 * ops.fs); % width of waveform templates (makes consistent N-ms either side of spike, regardless of sampling rate)
% ops.nt0                 = ops.nt0 + mod(ops.nt0-1, 2); % ensure oddity (...forum something about aiding spike trough alignment)
% ops.nt0 = 61;

% when Kilosort does CAR, is on batch time blocks, which mitigates risk of steps injected btwn demeaned segments
% -- NOTE: flag is "CAR", but uses common median referencing [CMR]
ops.CAR                 = 1;
ops.useStableMode       = 1;
ops.reorder             = 0; % this should always be disabled (==0); reordering time to fix probe drift was no good
%    if ~ops.reorder, fprintf(2,'\n!!!  NOTICE:  !!!\tKilosort batch reordering has been disabled in settings [ops.reorder==0]\n\n'); end


% try to address scaling discrepancy btwn data and template/whitened by increasing scaleproc; see get_whitening_matrix.m (def=200)
ops.scaleproc = 200;%600;    % [def=200]  % == 600 equates filtered/unfiltered amps of Plexon Omniplex recording in gui, BUT unsure if/when too high a value might start clipping underlying voltage sampling resolution

ops.throw_out_channels = 0; % NO! confounds source identity; never throw out chans during sort
ops.minfr_goodchannels = 0; % minimum firing rate on a "good" channel (0 to skip); always disable, see above   (def=0.1)

% [minFR] prevents errant 'units' with just a few detected spikes from proliferating
% - implementation is a little dicey with longer batch durations (5-10 sec) and/or shorter files (<1hr)
% - ...particularly when nbatches<=100, or templates that are added late in learning phase
ops.minFR = 0.02; 

% try non-zero min rate cutoff
% - clip truly useless templates, but don't drop less active ones (esp with randomized batch order during learning)
ops.minFR = ops.minFR / max([ops.learnRand*2,1]);


% threshold used when establishing baseline templates from raw data
% - standard codebase tends to [frustratingly] overwrite this param, but working to straighten out those instances
ops.spkTh = -6;     % [def= -6]

ops.ThPre = 8;      % [def= 8]

% splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
ops.AUCsplit = 0.9; %0.95; % ks2 default=0.9;  ks3 default=0.8;

% how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
ops.lam = 10;  % ks3 default is 20; previously 10...    (TBC: this has always been totally cryptic)

% threshold on projections (like in Kilosort1, can be different for last pass like [10 4])  (def=[10,4])
ops.Th = [10 4]; %[10 4];  (TBC: [8 4] better for awake nhp, but still clipping)
%               % eMouse benchmark config uses [6 2], but this output is garbage (pinned at max units, even high amp/snr units are over split)

% ~!~ KILOSORT 3.0 ~!~
% - appears to be sea change in how ops.Th is used, default is now [9 9]
% ops.Th = [8  4]; % no, vvery few units with [9 9]...must be throwing out lots of viable signal on oversampled Neuropixels


%% Stereo-probe specific adjustments (standard geom: 50um within, 100um between stereopairs)
% spatial constant in um for computing residual variance of spike     (def=30)
ops.sigmaMask = 70;  % 50-70 better for 50/100um intra/inter-trode spacing; else no spread across channels
ops.whiteningRange = 32; % use all chanels available

% Spike detection: .loc_range & .long_range are used isolate threshold peaks in:
%       [nSamples in time, nChannels];
% - BUT relevant uprobe 'channel' units are very different from nanopixel spacing
%       nChannels==3 will include lateral & longitudinally adjacent channels of stereo probe
%       nChannels==5 will include adjacent stereopairs (**but b/c spacing asymmetry of channel indices, this can include channel up to 300 microns away...:-/ )
ops.loc_range   = [5, 2];   % def=[5, 4]
ops.long_range  = [30, 4];  % def=[30, 6]


%% ---!!!--- NEW for kilosort 2.5,  "datashift" ---!!!---
ops.sig = 20; % [20] "spatial smoothness constant for registration"  (TBC: unclear how related to ops.sigmaMask ...try expanding slightly for uprobe)
% - cant be 0; used by standalone_detector.m & 0 will wipe out multi-scale gaussian spread of generic templates across channels
% So looks like this param (or the often hardcoded .sig param) is used when applying datashift drift corrections in increments smaller than [y] sampling of recording sites
% - this effectively blurrs shifted data traces into/across adjacent channels
% - maybe doing so flies with high res sampling of neuropixels, but abruptly jacks data quality/signal on more coarsely sampled devices (e.g. uprobes)
% - DONT expand, like we do for sigmaMask
% - Maybe minimize this as much as possible (==10), esp with [ops.integerShifts==1]  (see custom flags above)


%% Update temporal/batch parameters
% [.nTEMP] number of initial templates to extract from threshold crossings
% - These form the basis for any new templates added during learning & the initial PCA dimensions
% - If undefined, 6 is the usual number of templates, but more seems generally non-detrimental & likely helpful
ops.nTEMP = 12;

% number of samples to average over (annealed from first to second value)     (def=[20,400])
ops.momentum = [60 600]; % should really take batch size into account, but try expanding first [40 400];
% - this doesn't appear to be active anymore during actual extraction (trackAndSort.m)
% - ...it needs to be for template shapes to be dynamically updated (e.g.  rez.WU)
% - AH! another feature that was dropped completely between Kilosort 2.0 & 2.5/3.0!  ...now revived in czuba's [ks25] branch

ops.nfilt_factor        = 4; % (def=4) max number of clusters per ['good'] channel (even temporary ones)    

% TBC version define batches in seconds of data, not abstract bit chunks
batchSec                = 8;  % define batch number in seconds of data     (TBC: 10 seems good for 1-2 hr files and/or 32 channels)

% samples of symmetrical buffer for whitening and spike detection
% - implementation FINALLY FIXED in ks25!! (...was broken since kilosort 2.0)
ops.ntbuff              =  ceil(2*ops.fs/64)*64;%  ceil(batchSec/4*ops.fs/64)*64;%  64; % 64*300; % (def=64) 

% buffer size in samples
ops.NT                  = ceil(batchSec*ops.fs/32)*32; % convert to 32 count increments of samples

% sample from batches more sparsely (in certain circumstances/analyses)
batchSkips              = ceil(60/batchSec); % do high-level assessments at least once every minute of data
ops.nskip               = batchSkips;  % 1; % how many batches to skip for determining spike PCs
ops.nSkipCov            = batchSkips;  % 1; % compute whitening matrix from every N-th batch


%% Apply to updates to gui params object
fprintf('New ops settings:\n')
disp(ops)

ks.ops = ops;

% Update GUI parameter values
ks.H.settings.setFsEdt.String = num2str(ks.ops.fs);
ks.H.settings.setTrangeEdt.String = num2str(ks.ops.trange);
% !! CAUTION: jacked GUI parameter name correspondence for >=ks2.5
if any(contains(ks.H.settings.setMinfrTxt.String, 'blocks', 'ignoreCase',1))
    ks.H.settings.setMinfrEdt.String = num2str(ks.ops.nblocks); 
else
    ks.H.settings.setMinfrEdt.String = num2str(ks.ops.minfr_goodchannels);
end
ks.H.settings.setThEdt.String = num2str(ks.ops.Th);
ks.H.settings.setLambdaEdt.String = num2str(ks.ops.lam);
ks.H.settings.setCcsplitEdt.String = num2str(ks.ops.AUCsplit);

clear ops;

% set figDir in base workspace
figDir = fullfile(ks.ops.saveDir,'figs');

% move focus to command window
commandwindow;
% ks.updateFileSettings
fprintf(['\nDone.\tNew settings have been applied to kilosort GUI object [ks].\n',repmat('=-',1,20),'\n']);

% make entry in gui log
try
    ks.updateFileSettings;
    ks.log('Advanced Kilosort [& GUI] params updated.');
end
