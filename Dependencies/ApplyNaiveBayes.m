function [label,probability, Performance]=ApplyNaiveBayes(Tbl,Parameterkernels,label)
% Input:
% Tbl: table with all parameters (columns), rows for observations
% Parameterkernels: Probability distribution for every parameter, for every
% condition
% label: either label per observation OR the options for classes. 

% Output:
% label given Naive Bayes rule
% Probability of it being a certain class
% performance, if label is a label per observation

% Calculate posterior probability given Paremeterkernels (=probability
% distributions)
%% User input. Better not changed unless you do it everywhere
stepsize = 0.01;

%% Initialize
Edges = [0:stepsize:1];
ScoreVector = Edges(1)+stepsize/2:stepsize:Edges(end)-stepsize/2;

% Input:
if length(label)==size(Tbl,1)
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
Prior = nan(1,length(Cond));
for Ck = 1:length(Cond)
    Prior(Ck) = sum(label==Cond(Ck))/length(label);
end

%% p(X|Ck) Likelihood
% Find index for every observation in scorevector to assign correct
% probability
[mindif,minidx] = arrayfun(@(X) min(abs(X-ScoreVector)),table2array(Tbl),'Uni',0);
minidx = cell2mat(minidx);

% likelihood is the product of all probabilities given a certain condition
likelihood = nan(size(Tbl,1),length(Cond));
for Ck = 1:length(Cond)
    % Probability of observation given a given class
    tmpp = arrayfun(@(Y) cell2mat(arrayfun(@(X) Parameterkernels(minidx(X,Y),Y,Ck),1:size(minidx,1),'Uni',0)),1:size(minidx,2),'Uni',0);
    tmpp = cat(1,tmpp{:});
    % Naive Bayes assumes independence: product of the probabilities
    likelihood(:,Ck) = prod(tmpp,1);
end

%% evidence p(x)
Evidence = nansum((Prior.*likelihood),2);

%% Probability
probability = nan(size(Tbl,1),2);
for Ck=1:length(Cond)
    probability(:,Ck) = (Prior(Ck).*likelihood(:,Ck))./Evidence;
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
