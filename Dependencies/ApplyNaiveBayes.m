function [label,probability, Performance]=ApplyNaiveBayes(Tbl,Parameterkernels,label)
% Calculate posterior probability given Paremeterkernels (=probability
% distributions)
stepsize = 0.01;
Edges = [0:stepsize:1];
ScoreVector = Edges(1)+stepsize/2:stepsize:Edges(end)-stepsize/2;

% Bayes Theorem
% p(Ck) is irrelevant
% p(X|Ck)
% Find index in scorevector
[mindif,minidx] = arrayfun(@(X) min(abs(X-ScoreVector)),table2array(Tbl),'Uni',0);
minidx = cell2mat(minidx);

% likelihood is the product of all probabilities given a certain condition
likelihood = nan(size(Tbl,1),2);
for Ck = 1:2
    % Probability of observation given a given class
    tmpp = arrayfun(@(Y) cell2mat(arrayfun(@(X) Parameterkernels(minidx(X,Y),Y,Ck),1:size(minidx,1),'Uni',0)),1:size(minidx,2),'Uni',0);
    tmpp = cat(1,tmpp{:});
    % Naive Bayes assumes independence: product of the probabilities
    likelihood(:,Ck) = prod(tmpp,1);
end

probability = nan(size(Tbl,1),2);
for Ck=1:2
    probability(:,Ck) = likelihood(:,Ck)./(nansum(likelihood,2));
end
if nargin>2
    originallabel = label;
end
label = false(size(probability,1),1);
label(probability(:,2)>0.5) = 1;

if nargin>2
    Performance = [sum(originallabel==1&label==1)./sum(originallabel) sum(originallabel==0&label==0)./sum(originallabel==0)];
else
    Performance = nan;
end

end
