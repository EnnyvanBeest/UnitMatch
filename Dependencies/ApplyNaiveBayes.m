function [label,probability, Performance]=ApplyNaiveBayes(Tbl,Parameterkernels,label,Prior)
    % Input:
    % Tbl: table with all parameters (columns), rows for observations
    % Parameterkernels: Probability distribution for every parameter, for every
    % condition
    % label: either label per observation OR the options for classes. 
    % Prior: if empty, use 0.5 0.5
    
    % Output:
    % label given Naive Bayes rule
    % Probability of it being a certain class
    % performance, if label is a label per observation
    
    % Calculate posterior probability given Paremeterkernels (=probability
    % distributions)
    %% User input. Better not changed unless you do it everywhere
    nPairs = size(Tbl,1);
    stepsize = 1./size(Parameterkernels,1);

    %% Initialize
    Edges = 0:stepsize:1;
    ScoreVector = Edges(1)+stepsize/2:stepsize:Edges(end)-stepsize/2;
    
    % Input:
    if length(label)==nPairs
        Cond = unique(label);
        GivePerformance=1;
    elseif length(label)==length(unique(label))
        Cond = unique(label);
        GivePerformance=0;
    else
        disp('I don''t understand this input')
        keyboard
    end
    
    
    %% Bayes Theorem
    %% p(Ck) Prior
    if nargin<4 || ~exist('Prior','var')
        Prior = [0.5 0.5]; %equal assumptions instead
    end

    %% p(X|Ck) Likelihood
    % Find index for every observation in scorevector to assign correct
    % probability
    chunkSize = 1e6;
    nChunk = 1;
    minidx = nan(size(Tbl));
    while nChunk*chunkSize < nPairs + chunkSize
        idx = (nChunk-1)*chunkSize+1:min((nChunk-1)*chunkSize+chunkSize,nPairs);
        x1 = repmat(table2array(Tbl(idx,:)), [1 1 numel(ScoreVector)]);
        x2 = repmat(permute(ScoreVector, [1 3 2]), [size(x1,1), size(x1,2), 1]);
        [~,minidx(idx,:)] = min(abs(x1-x2),[],3);
        clear x1 x2
        nChunk = nChunk + 1;
    end

    % likelihood is the product of all probabilities given a certain condition
    likelihood = nan(nPairs,length(Cond));
    for Ck = 1:length(Cond)
        % Probability of observation given a given class
        tmpp = nan(size(minidx));
        for yy = 1:size(minidx,2)
            tmpp(:,yy) = Parameterkernels(minidx(:,yy),yy,Ck);
        end
        % Naive Bayes assumes independence: product of the probabilities
        likelihood(:,Ck) = prod(tmpp,2);
    end
    
    %% Probability
    probability = nan(nPairs,2);
    for Ck=1:length(Cond)
        probability(:,Ck) = (Prior(Ck).*likelihood(:,Ck))./nansum((Prior.*likelihood),2);
    end
    if GivePerformance
        originallabel = label;
    end
    label = false(size(probability,1),1);
    label(probability(:,2)>0.5) = 1;
    
    if GivePerformance
        Performance = nan(1,length(Cond));
        for Ck=1:length(Cond)
            Performance(Ck) = sum(originallabel==Cond(Ck)&label==Cond(Ck))./sum(originallabel==Cond(Ck));
        end    
    else
        Performance = nan;
    end
end
