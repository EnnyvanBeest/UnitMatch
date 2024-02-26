function [corrWCCA, W, pUall] = computeNatImCorr(spikeData,cc2keep)
    %% Projects neural responses in lower-dimensional space
    % This script will perform CCA across days and use the CCA weights of
    % each neurons along the first "cc2keep" components to compute a
    % similarity score (correlation) across sessions.
    %
    % Parameters:
    % -------------------
    % spikeData: cell
    %   Data array with binned PSTHs for all clusters (of stim ID x 
    % time x neurons x repeats)
    % cc2keep (optional): int
    %   Number of Canonical components to keep
    %
    % Returns:
    % -------------------
    % corrWCCA: cell
    %   Contains all the correlation values of CCA weights of the neurons 
    %   across days.
    % W: cell
    %   Contains the CCA weights

    if ~exist('cc2keep','var')
        cc2keep = 1:75; % arbitrary number
    end
    
    % perform CCA
    Xwall = [];
    idxall = [];
    clear U S V
    for k = 1:numel(spikeData)
        if ~isempty(spikeData{k})
            s = size(spikeData{k});
            X = reshape(nanmean(spikeData{k},4),[s(1)*s(2),s(3)]);
            % Xz = zscore(X);
            [U{k},S{k},V{k}] = svd(X, 'econ');
            Xw = U{k}*V{k}';
            Xwall = cat(2,Xwall, Xw);
            idxall = [idxall, ones(1,s(3))*k];
        end
    end
    [Uall,~,~] = svd(Xwall-mean(Xwall), 'econ');

    cc2keep = cc2keep(cc2keep<=size(Uall,2));
    
    % get weight of each neuron on each day, for each component
    pUall = pinv(Uall(:,cc2keep));
    for k = 1:numel(spikeData)
        if ~isempty(spikeData{k})
            s = size(spikeData{k});
            X = reshape(nanmean(spikeData{k},4),[s(1)*s(2),s(3)]);
            W{k} = (pUall*X)'; % directly recompute a linear regression of Uall onto the data
        else 
            W{k} = [];
        end
    end

    %% look at correlation of all cells with all the other cells, to see if any is stable...

    corrWCCA = cell(numel(spikeData),numel(spikeData));
    for d1 = 1:numel(spikeData)
        for d2 = 1:numel(spikeData)
            if ~isempty(W{d1}) && ~isempty(W{d2})
                corrWCCA{d1,d2} = corr(W{d1}(:,cc2keep)',W{d2}(:,cc2keep)');
            else
                corrWCCA{d1,d2} = nan;
            end
        end
    end
    
end
