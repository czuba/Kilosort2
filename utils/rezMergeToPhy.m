function rezMergeToPhy(rez1, rez2, savePath)
% merge two kilosort rez structs into one,
% save all requisit output files for loading merged dataset into Phy
%
% ---
%   W.I.P.:: does not handle template or feature projections, which
%   are typically excluded from rez.mat save struct, and should really
%   be recomputed based on merged content (e.g. template similarity 
%   & feature projections of coherent clusters from each rez session)
% ---
% 2021-06-xx  TBC  Hacked together based on standard rezToPhy.m
%                  
% 

%% st3 content:
% % % % From learnAndSolve8b >> runTemplates >> trackAndSort.m
% % %     st3(irange,1) = double(st); % spike times
% % %     st3(irange,2) = double(id0+1); % spike clusters (1-indexing)
% % %     st3(irange,3) = double(x0); % template amplitudes
% % %     st3(irange,4) = double(vexp); % residual variance of this spike
% % %     st3(irange,5) = ibatch; % batch from which this spike was found

%% collect & sort spike vars in local workspace

if ~exist(savePath,'dir')
    mkdir(savePath);
elseif strcmp(rez1.ops.saveDir, savePath) || strcmp(rez2.ops.saveDir, savePath)
    error('Merged destination directory cannot be the same as either of the input rez structs.');
elseif exist(savePath,'dir')
    savePath = uigetdir(savePath, 'Destination exists, please confirm');
end
    
    
% add index for each rez struct
rez1.rid = 1;
rez2.rid = 2;

% combine rez structs
rez(1) = rez1;
rez(2) = rez2;
ntemps = cumsum([0, arrayfun(@(x) length(x.mu), rez)]);


%% clear input rez vars (excess memory overhead)
clear rez1 rez2

%     % already done when rez.mat created 
%     % spikeTimes will be in samples, not seconds
%     rez.W = gather(single(rez.Wphy));
%     rez.U = gather(single(rez.U));
%     rez.mu = gather(single(rez.mu));


%     % already done when rez.mat created 
%     [~, isort]   = sort(rez.st3(:,1), 'ascend');
%     rez.st3      = rez.st3(isort, :);
%     rez.cProj    = rez.cProj(isort, :);
%     rez.cProjPC  = rez.cProjPC(isort, :, :);

%% clear existing/conflicting files from destination
fs = dir(fullfile(savePath, '*.npy'));
for i = 1:length(fs)
   delete(fullfile(savePath, fs(i).name));
end
if exist(fullfile(savePath, '.phy'), 'dir')
    rmdir(fullfile(savePath, '.phy'), 's');
end


%% Compile params from input rez structs
spikeTimes = cell2mat(arrayfun(@(x) uint64(x.st3(:,1)), rez, 'uni',0)');
% spikeTimes = cell2mat(spikeTimes'); % concatenate

[spikeTimes, ii] = sort(spikeTimes);

% - add offset to template indices of second rez struct to ensure ids are unique
% - offset must match with index of concatenated template shapes as well
spikeTemplates = cell2mat(arrayfun(@(x) uint32(x.st3(:,2) + ntemps(x.rid)), rez, 'uni',0)');
spikeTemplates = spikeTemplates(ii);
% NO:  st3(:,5) is really batch#, not cluster# (!??...KS1 holdover?)
% if size(rez.st3,2)>4
%     spikeClusters = uint32(1+rez.st3(:,5));
% end

... unused:     spikeBatch = uint32(rez.st3(:,5)); 

amplitudes = cell2mat(arrayfun(@(x) x.st3(:,3), rez, 'uni',0)');
amplitudes = amplitudes(ii);
% Calc amplitudes to reflect temporal variations in waveform templates
isgood = cell2mat(arrayfun(@(x) x.good, rez, 'uni',0)');

estContam = cell2mat(arrayfun(@(x) x.est_contam_rate, rez, 'uni',0)');

% the following fields MUST BE IDENTICAL for both rez structs
Nchan = rez(1).ops.Nchan;

xcoords     = rez(1).xcoords(:);
ycoords     = rez(1).ycoords(:);
chanMap     = rez(1).ops.chanMap(:);
chanMap0ind = chanMap - 1;

nt0 = size(rez(1).W,1);


U = arrayfun(@(x) x.U, rez,'uni',0);
U = cat(2, U{:}); % must do two step for multi dimensional
W = arrayfun(@(x) x.W, rez,'uni',0);
W = cat(2, W{:});

% total number of templates
Nfilt = ntemps(end);% size(W,2);

templates = zeros(Nchan, nt0, Nfilt, 'single');
for iNN = 1:size(templates,3)
   templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))';
end
templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
templatesInds = repmat((0:size(templates,3)-1), size(templates,1), 1); % we include all channels so this is trivial

% cProj & cProjPC fields are typically excluded from rez.mat save
% because can balloon file into gigs of data
%     if isfield(rez, 'cProj') && all(arrayfun(@(x) ~isempty(x.cProj), rez))
%         templateFeatures = rez.cProj;
%         templateFeatureInds = uint32(rez.iNeigh);
%         pcFeatures = rez.cProjPC;
%         pcFeatureInds = uint32(rez.iNeighPC);
%     end

% Here things get tricky
% whiteningMatrix = rez.Wrot/rez.ops.scaleproc;
% whiteningMatrix = eye(size(rez.Wrot)) / rez.ops.scaleproc;
% ...luckily, this isn't actually using the whitening matrix, just undoing the scaleproc
whiteningMatrix = eye(size(rez(1).Wrot)) / rez(1).ops.scaleproc;
whiteningMatrixInv = whiteningMatrix^-1;


%% This section should all 'just work' on the concatenated data
% here we compute the amplitude of every template...

% unwhiten all the templates
tempsUnW = zeros(size(templates));
for t = 1:size(templates,1)
    tempsUnW(t,:,:) = squeeze(templates(t,:,:))*whiteningMatrixInv;
end

% The amplitude on each channel is the positive peak minus the negative
tempChanAmps = squeeze(max(tempsUnW,[],2))-squeeze(min(tempsUnW,[],2));

% The template amplitude is the amplitude of its largest channel
tempAmpsUnscaled = max(tempChanAmps,[],2);

% assign all spikes the amplitude of their template multiplied by their
% scaling amplitudes
spikeAmps = tempAmpsUnscaled(spikeTemplates).*amplitudes;

% take the average of all spike amps to get actual template amps (since
% tempScalingAmps are equal mean for all templates)
ta = clusterAverage(spikeTemplates, spikeAmps);
tids = unique(spikeTemplates);
tempAmps = zeros(ntemps(end),1);       % zeros(numel(rez.mu),1);
tempAmps(tids) = ta; % because ta only has entries for templates that had at least one spike
tempAmps = tempAmps';   % gain is fixed
%     gain = getOr(rez.ops, 'gain', 1);
%     tempAmps = gain*tempAmps'; % for consistency, make first dimension template number

if ~isempty(savePath)
    fileID = fopen(fullfile(savePath, 'cluster_KSLabel.tsv'),'w');
    fprintf(fileID, 'cluster_id%sKSLabel', char(9));
    fprintf(fileID, char([13 10]));
    
    fileIDCP = fopen(fullfile(savePath, 'cluster_ContamPct.tsv'),'w');
    fprintf(fileIDCP, 'cluster_id%sContamPct', char(9));
    fprintf(fileIDCP, char([13 10]));
    
    fileIDA = fopen(fullfile(savePath, 'cluster_Amplitude.tsv'),'w');
    fprintf(fileIDA, 'cluster_id%sAmplitude', char(9));
    fprintf(fileIDA, char([13 10]));
        
    for j = 1:length(isgood)
        if isgood(j)
            fprintf(fileID, '%d%sgood', j-1, char(9));
        else
            fprintf(fileID, '%d%smua', j-1, char(9));
        end
        fprintf(fileID, char([13 10]));
        
        if isfield(rez, 'est_contam_rate')
            fprintf(fileIDCP, '%d%s%.1f', j-1, char(9), estContam(j)*100);
            fprintf(fileIDCP, char([13 10]));
        end
        
        fprintf(fileIDA, '%d%s%.1f', j-1, char(9), tempAmps(j));
        fprintf(fileIDA, char([13 10]));
        
    end
    fclose(fileID);
    fclose(fileIDCP);
    fclose(fileIDA);
    
    
    writeNPY(spikeTimes,    fullfile(savePath, 'spike_times.npy'));
    writeNPY(uint32(spikeTemplates-1),  fullfile(savePath, 'spike_templates.npy')); % -1 for zero indexing

    writeNPY(amplitudes,    fullfile(savePath, 'amplitudes.npy'));
    writeNPY(templates,     fullfile(savePath, 'templates.npy'));
    writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));

    %chanMap0ind = int32(chanMap0ind);
    chanMap0ind = int32([1:Nchan]-1);
    writeNPY(chanMap0ind,   fullfile(savePath, 'channel_map.npy'));
    writeNPY([xcoords ycoords],     fullfile(savePath, 'channel_positions.npy'));
    
    
%     % excluded from rez merge....
%     writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
%     writeNPY(templateFeatureInds'-1, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
%     writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
%     writeNPY(pcFeatureInds'-1, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing


    writeNPY(whiteningMatrix,       fullfile(savePath, 'whitening_mat.npy'));
    writeNPY(whiteningMatrixInv,    fullfile(savePath, 'whitening_mat_inv.npy'));

    if isfield(rez, 'simScore')
        % similarTemplates = cell2mat(arrayfun(@(x) x.simScore, rez, 'uni',0)');
        similarTemplates = zeros(ntemps(end));
        sims = arrayfun(@(x) x.simScore, rez, 'uni',0);
        for i = 1:length(sims)
            nt = size(sims{i},1);
            ii = (1:nt)+ntemps(i);
            similarTemplates(ii,ii) = sims{i};
        end
        writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
    end


    % Duplicate "KSLabel" as "group", a special metadata ID for Phy, so that
    % filtering works as expected in the cluster view
    KSLabelFilename = fullfile(savePath, 'cluster_KSLabel.tsv');
    copyfile(KSLabelFilename, fullfile(savePath, 'cluster_group.tsv'));

    %     % if raw/binary data file location is not same as save destination,
    %     % attempt to create symlink to raw file
    %     if ~strcmpi( fileparts(rez.ops.fbinary), rez.ops.saveDir)
    %         fprintf(2, ['\n\tWARNING: raw data directory and save output data directory are distinct locations.'...
    %             '\n\tAttempt to create symlink to raw data in save output directory...']);
    %         try
    %             [~, fname, ext] = fileparts(rez.ops.fbinary);
    %             [err, msg] = system( sprintf('ln -sv %s %s', rez.ops.fbinary, fullfile(rez.ops.saveDir, [fname ext]) ));
    %             if ~err
    %                 fprintf('successful!\n\t%s\n',msg)
    %             else
    %                 % Note: symlinks won't work on certain file systems (needs extended attributes; not Fat32)
    %                 fprintf(2, 'failed.\n\t%s',msg)
    %                 fprintf(['\n\t%s','\n\t>>','\n\t%s'...
    %                     '\n\tA copy of raw data may need to be added to output directory before starting Phy\n\n'], rez.ops.fbinary,rez.ops.saveDir);
    %             end
    %         end
    %     end
    
    %make params file
    if ~exist(fullfile(savePath,'params.py'),'file')
        fid = fopen(fullfile(savePath,'params.py'), 'w');

        % use relative path name for preprocessed data file in params.py
        % - defaults to preprocessed file of last rez struct
        % - assuming they're in order
        % - **** AND that the preprocessed data file was created with the [ks25] branch
        %   - which include everything in the preprocessed file from t0 to tend
        [~, fname, ext] = fileparts(rez(end).ops.fproc);
        copyfile(rez(end).ops.fproc, fullfile(savePath, [fname,ext]));
        fprintf(fid, 'dat_path = ''%s''\n', fullfile('.',[fname,ext]));
        fprintf(fid,'n_channels_dat = %i\n', Nchan);
        fprintf(fid,'dtype = ''int16''\n');
        fprintf(fid,'offset = 0\n');
        if mod(rez(1).ops.fs,1)
            fprintf(fid,'sample_rate = %i\n', rez(1).ops.fs);
        else
            fprintf(fid,'sample_rate = %i.\n', rez(1).ops.fs);
        end
        fprintf(fid,'hp_filtered = True\n');
        fprintf(fid,'template_scaling = 5.0\n');
        
        fclose(fid);
    end
end

end %main function
