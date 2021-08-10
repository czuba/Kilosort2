function rez = learnAndSolve8b(rez)
tic

ops = rez.ops;
ops.fig = getOr(ops, 'fig', 1);

NrankPC = 6;  
Nrank = 3;
rng('default'); rng(1);

[wTEMP, wPCA]    = extractTemplatesfromSnippets(rez, NrankPC);
% wPCA    = extractPCfromSnippets(rez, Nrank);
% wPCA(:,1) = - wPCA(:,1) * sign(wPCA(20,1));


wPCA = gpuArray(wPCA(:, 1:Nrank));
wTEMP = gpuArray(wTEMP);

wPCAd = double(wPCA);

ops.wPCA = gather(wPCA);
ops.wTEMP = gather(wTEMP);
rez.ops = ops;

NchanNear   = min(ops.Nchan, 32);
Nnearest    = min(ops.Nchan, 32);

sigmaMask  = ops.sigmaMask;


% ops.spkTh = -6; % why am I overwriting this here?
% I dont know...maybe just don't  --TBC


nt0 = ops.nt0;
nt0min  = rez.ops.nt0min; 

nBatches  = rez.temp.Nbatch;
NT  	= ops.NT;
batchstart = 0:NT:NT*nBatches;
Nfilt 	= ops.Nfilt;

Nchan 	= ops.Nchan;

% [iC, mask] = getClosestChannels(rez, sigmaMask, NchanNear);
[iC, mask, C2C] = getClosestChannels(rez, sigmaMask, NchanNear);

isortbatches = rez.iorig(:);
nhalf = floor(nBatches/2);   % nhalf = ceil(nBatches/2);


% Setup double pass sequence through batches of time
% (**NOTE relative times in this sweep are in resorted 'time'**)
% 1st Pass:  [midpoint-to-end, midpoint-to-start]
% 2nd Pass:  [midpoint-to-start, midpoint-to-end]
% % % pass1 = [(nhalf+1):nBatches, nBatches:-1:1];
pass1 = [1:1:nBatches, (nBatches-1):-1:nhalf+1];%[(nhalf+1):nBatches, nBatches:-1:1];
% --NOTE: this  

% sweep in opposite directions for second pass
% % % pass2 = [1:1:nBatches];%]circshift(pass1, nhalf+1);
pass2 = [nhalf:-1:1, nhalf+1:1:nBatches];% [1:1:nBatches];
irounds = [pass1, pass2];
passId = [ones(size(pass1)), 2*ones(size(pass2))];

% % git orig code
% ischedule = [nhalf:nBatches nBatches:-1:nhalf];
% i1 = [(nhalf-1):-1:1];
% i2 = [nhalf:nBatches];
% irounds = cat(2, ischedule, i1, i2);

niter   = numel(irounds);

% replace "niter-nBatches" with direct count of 1st pass instances
nPass1batches = sum(passId==1);

% % % if irounds(niter - nBatches)~=nhalf
% % %     error('mismatch between number of batches');
% % % end

%
flag_final = 0;
flag_resort      = 1;

t0 = ceil(rez.ops.trange(1) * ops.fs);

nInnerIter  = 60;


pmi = exp(-1./linspace(ops.momentum(1), ops.momentum(2), nPass1batches));

Nsum = 7; % how many channels to extend out the waveform in mexgetspikes
Params     = double([NT Nfilt ops.Th(1) nInnerIter nt0 Nnearest ...
    Nrank ops.lam pmi(1) Nchan NchanNear ops.nt0min 1 Nsum NrankPC ops.Th(1)]);

W0 = permute(double(wPCA), [1 3 2]);

iList = int32(gpuArray(zeros(Nnearest, Nfilt)));

nsp = gpuArray.zeros(0,1, 'double');

Params(13) = 0;

[Ka, Kb] = getKernels(ops, 10, 1);

p1 = .95; % decay of nsp estimate

fprintf('Time %3.0fs. Optimizing templates ...\n', toc)

T0 = tic;

fid = fopen(ops.fproc, 'r');

ntot = 0;
ndrop = zeros(1,2);

m0 = ops.minFR * ops.NT/ops.fs;

% 1st pass == ibatch<nPass1batches
% 2nd pass == ibatch>nPass1batches

for ibatch = 1:niter
    if size(nsp,2)>1
        nsp = nsp(:);
        fprintf(2, '!');    %  keyboard
    end
    
    %     k = irounds(ibatch);
    korder = irounds(ibatch);
    k = isortbatches(korder);
    % which pass are we on?
    passNo = passId(ibatch);
    
    % 2nd pass == ibatch>nPass1batches
    %    Drift from midpoint to tEnd, jump back to midpoint (here), then
    %    works from midpoint to t0
    
    if passNo>1 && korder==nhalf+1 %(korder==nhalf || korder==nhalf+1)    %ibatch>nPass1batches && korder==nhalf
        [W, dWU] = revertW(rez);
        fprintf('reverted back to middle timepoint [%d, %d]\n', ibatch, korder)
%         keyboard
    end

    % 1st pass:  Apply variable annealing strength (average of last n-spikes)
    if passNo==1    %ibatch<=nPass1batches
        Params(9) = pmi(ibatch);
        pm = pmi(ibatch) * gpuArray.ones(Nfilt, 1, 'double');
    end

    % dat load \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    offset = 2 * ops.Nchan*batchstart(k);
    fseek(fid, offset, 'bof');
    dat = fread(fid, [NT ops.Nchan], '*int16');
    dataRAW = single(gpuArray(dat))/ ops.scaleproc;

    % Initial pass only
    if ibatch==1            
        [dWU, cmap] = mexGetSpikes2(Params, dataRAW, wTEMP, iC-1);        
%         dWU = mexGetSpikes(Params, dataRAW, wPCA);
        dWU = double(dWU);
        dWU = reshape(wPCAd * (wPCAd' * dWU(:,:)), size(dWU));
        
        
        W = W0(:,ones(1,size(dWU,3)),:);
        Nfilt = size(W,2);
        nsp(1:Nfilt) = m0;        
        Params(2) = Nfilt;
    end
    
    if flag_resort
        [~, iW] = max(abs(dWU(nt0min, :, :)), [], 2);
        iW = int32(squeeze(iW));

        [iW, isort] = sort(iW);
        W = W(:,isort, :);
        dWU = dWU(:,:,isort);
        nsp = nsp(isort);        
    end

    % decompose dWU by svd of time and space (61 by 61)
    [W, U, mu] = mexSVDsmall2(Params, dWU, W, iC-1, iW-1, Ka, Kb);
  
    % this needs to change
    [UtU, maskU] = getMeUtU(iW, iC, mask, Nnearest, Nchan);
   
    
    [st0, id0, x0, featW, dWU0, drez, nsp0, featPC, vexp] = ...
        mexMPnu8(Params, dataRAW, single(U), single(W), single(mu), iC-1, iW-1, UtU, iList-1, ...
        wPCA);
   
    fexp = exp(double(nsp0).*log(pm(1:Nfilt)));
    fexp = reshape(fexp, 1,1,[]);
    nsp = nsp(:) * p1 + (1-p1) * double(nsp0);
    dWU = dWU .* fexp + (1-fexp) .* (dWU0./reshape(max(1, double(nsp0)), 1,1, []));
    
    % \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    % last batch of this pass: Lock in settings
    %  & adjust flags for next pass
    if passNo~=passId(min(ibatch+1,end)) % ibatch==nPass1batches
        flag_resort   = 0;
        flag_final = 1;
        
        % final clean up
        [W, U, dWU, mu, nsp, ndrop] = ...
            triageTemplates2(ops, iW, C2C, W, U, dWU, mu, nsp, ndrop);
        fprintf('Params = \n\t%s\n',sprintf('%04.3f\t',Params))
        
        Nfilt = size(W,2);
        Params(2) = Nfilt;

        [WtW, iList] = getMeWtW(single(W), single(U), Nnearest);

        [~, iW] = max(abs(dWU(nt0min, :, :)), [], 2);
        iW = int32(squeeze(iW));

        % extract ALL features on the last pass
        Params(13) = 2;

        % different threshold on last pass?
        Params(3) = ops.Th(end);

        rez = memorizeW(rez, W, dWU, U, mu);
        fprintf('memorized middle timepoint [%d, %d]\n', ibatch, korder)
        fprintf('Params = \n\t%s\n',sprintf('%04.3f\t',Params))

    %end 
    
    elseif passNo==1 %ibatch<nPass1batches %-50

        %         if ibatch==nhalf
        %             fprintf('Params = \n\t%s\n',sprintf('%04.3f\t',Params))
        %             Params([3]) = ops.Th(end);
        %             fprintf('Params = \n\t%s\n',sprintf('%04.3f\t',Params))
        %         end

        % 1st pass stuff
        if rem(ibatch, 5)==1
            % this drops templates
            [W, U, dWU, mu, nsp, ndrop] = ...
                triageTemplates2(ops, iW, C2C, W, U, dWU, mu, nsp, ndrop);
        end
        Nfilt = size(W,2);
        Params(2) = Nfilt;

        % this adds templates        
%         dWU0 = mexGetSpikes(Params, drez, wPCA);
        [dWU0,cmap] = mexGetSpikes2(Params, drez, wTEMP, iC-1);
        
        if size(dWU0,3)>0    
            dWU0 = double(dWU0);
            dWU0 = reshape(wPCAd * (wPCAd' * dWU0(:,:)), size(dWU0));            
            dWU = cat(3, dWU, dWU0);

            W(:,Nfilt + [1:size(dWU0,3)],:) = W0(:,ones(1,size(dWU0,3)),:);

            nsp(Nfilt + [1:size(dWU0,3)]) = ops.minFR * NT/ops.fs;
            mu(Nfilt + [1:size(dWU0,3)])  = 10;            

            Nfilt = min(ops.Nfilt, size(W,2));
            Params(2) = Nfilt;

            W   = W(:, 1:Nfilt, :);
            dWU = dWU(:, :, 1:Nfilt);
            nsp = nsp(1:Nfilt);
            mu  = mu(1:Nfilt);            
        end

    %end

    elseif passNo==2 %ibatch>nPass1batches        
        
        %         if korder==nhalf+5
        %             fprintf('Params = \n\t%s\n',sprintf('%04.3f\t',Params))
        %             Params([3,16]) = ops.Th(end);
        %             fprintf('Params = \n\t%s\n',sprintf('%04.3f\t',Params))
        %         end

        % 2nd pass stuff
        rez.WA(:,:,:,k) = gather(W);
        rez.UA(:,:,:,k) = gather(U);
        rez.muA(:,k) = gather(mu);
        
        ioffset         = ops.ntbuff;
        if k==1
            ioffset         = 0;
        end
        toff = nt0min + t0 -ioffset + (NT-ops.ntbuff)*(k-1);        
        
        st = toff + double(st0);
        irange = ntot + [1:numel(x0)];
        
        if ntot+numel(x0)>size(st3,1)
           fW(:, 2*size(st3,1))    = 0;
           fWpc(:,:,2*size(st3,1)) = 0;
           st3(2*size(st3,1), 1)   = 0;
        end
        
        st3(irange,1) = double(st);
        st3(irange,2) = double(id0+1);
        st3(irange,3) = double(x0);        
        st3(irange,4) = double(vexp);
        st3(irange,5) = korder;

        fW(:, irange) = gather(featW);

        fWpc(:, :, irange) = gather(featPC);

        ntot = ntot + numel(x0);
    end
    
% %     if ibatch>1 && passNo==1 &&  korder==nhalf+1
% %         rez = memorizeW(rez, W, dWU, U, mu);
% %         fprintf('memorized middle timepoint [%d, %d]\n', ibatch, korder)
% %     end

    % clean up for next pass type
    if passNo~=passId(min(ibatch+1,end))  %ibatch==nPass1batches        
        st3 = zeros(1e7, 5);
        rez.WA = zeros(nt0, Nfilt, Nrank,nBatches,  'single');
        rez.UA = zeros(Nchan, Nfilt, Nrank,nBatches,  'single');
        rez.muA = zeros(Nfilt, nBatches,  'single');
        
        fW  = zeros(Nnearest, 1e7, 'single');
        fWpc = zeros(NchanNear, Nrank, 1e7, 'single');
    end

    doRefresh = (ibatch<10) || (ibatch<20 && rem(ibatch, 10)==1) || (ibatch<500 && rem(ibatch, 20)==1) || (rem(ibatch, 100)==1);
    if doRefresh    %rem(ibatch, 100)==1
        tstring = datestr(toc(T0)/86400,'HH:MM:SS'); % sec to days to string
        fprintf('  %s, %d / %d batches, %d units, nspks: %2.4f, mu: %2.4f, nst0: %d, merges: %2.4f, %2.4f\n', ...
            tstring, ibatch, niter, Nfilt, sum(nsp), median(mu), numel(st0), ndrop)
        
        if ops.fig
            
            if ibatch==1
                figHand = figure;
            else
                figure(figHand);
            end
            
            make_fig(W, U, mu, nsp);
       end
    end
end

fclose(fid);


st3 = st3(1:ntot, :);
fW = fW(:, 1:ntot);
fWpc = fWpc(:,:, 1:ntot);

% ntot

% [~, isort] = sort(st3(:,1), 'ascend');
% fW = fW(:, isort);
% fWpc = fWpc(:,:,isort);
% st3 = st3(isort, :);

rez.st3 = st3;
rez.st2 = st3;

rez.simScore = gather(max(WtW, [], 3));

rez.cProj    = fW';
rez.iNeigh   = gather(iList);

rez.ops = ops;

rez.nsp = nsp;

% nNeighPC        = size(fWpc,1);
rez.cProjPC     = permute(fWpc, [3 2 1]); %zeros(size(st3,1), 3, nNeighPC, 'single');

% [~, iNch]       = sort(abs(rez.U(:,:,1)), 1, 'descend');
% maskPC          = zeros(Nchan, Nfilt, 'single');
rez.iNeighPC    = gather(iC(:, iW));


nKeep = 20; % how many PCs to keep
rez.W_a = zeros(nt0 * Nrank, nKeep, Nfilt, 'single');
rez.W_b = zeros(nBatches, nKeep, Nfilt, 'single');
rez.U_a = zeros(Nchan* Nrank, nKeep, Nfilt, 'single');
rez.U_b = zeros(nBatches, nKeep, Nfilt, 'single');
for j = 1:Nfilt    
    WA = reshape(rez.WA(:, j, :, :), [], nBatches);
    WA = gpuArray(WA);
    [A, B, C] = svdecon(WA);
    rez.W_a(:,:,j) = gather(A(:, 1:nKeep) * B(1:nKeep, 1:nKeep));
    rez.W_b(:,:,j) = gather(C(:, 1:nKeep));
    
    UA = reshape(rez.UA(:, j, :, :), [], nBatches);
    UA = gpuArray(UA);
    [A, B, C] = svdecon(UA);
    rez.U_a(:,:,j) = gather(A(:, 1:nKeep) * B(1:nKeep, 1:nKeep));
    rez.U_b(:,:,j) = gather(C(:, 1:nKeep));
end

fprintf('Time %3.0fs. Finished compressing time-varying templates \n',toc)
%%



