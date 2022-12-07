function [InvPerformance, PerformanceRank, PerformanceLabels] = TwoFactorPerformance(Parameters,Predictors,Rankscore,DrawThis,Labels)
% Function to calculate the performance of the matching scores, based on
% parameter weights and predictors. Labels are assumed labels from PyKS,
% CCFingerprint is the cross-correlation fingerprint
% Parameters: threshold  and Weights for predictors (in that order!)

if nargin>4
    % Labels 'known'
    knownlabel = 1;
else
    knownlabel = 0;
end
% Calculate total matching score basedon predictors, weighted by their
% weight
TotalMatchingScore=arrayfun(@(X) Predictors(:,:,X).*Parameters(X+1),1:size(Predictors,3),'Uni',0);
TotalMatchingScore = nansum(cat(3,TotalMatchingScore{:}),3);

% Find identified matches, dependent on the threshold
IdentifiedMatches = TotalMatchingScore>Parameters(1);
RankBin = Rankscore==1;
% PerformanceCorrelations
PerformanceRank = (sum(abs(IdentifiedMatches(:)-RankBin(:))<1))...
    ./length(IdentifiedMatches(:));
% Punish off-diagonal matches unless rank is 1, and punish diagonal matches
% with low ranks
InvPerformance = sum(abs(diag(IdentifiedMatches)-diag(Rankscore))) + 0.5*sum(sum((triu(IdentifiedMatches,1)-triu(Rankscore,1))>0)); %

% PerformanceRank = sum((IdentifiedMatches(:) == 1 & Rankscore(:)==1))...
%     ./sum(IdentifiedMatches(:) | Rankscore(:)); 

% PerformanceCorrelations = corr(,TotalMatchingScore(:)); %Variance explained by cross-correlation fingerprints

% PerformanceLabels
if knownlabel
    PerformanceLabels= sum(IdentifiedMatches(:) == Labels(:))...
    ./length(IdentifiedMatches(:)); 
    %     PerformanceLabels= sum((IdentifiedMatches(:) == 1 & Labels(:)==1))...
    %     ./sum(IdentifiedMatches(:));
    InvPerformance = sum(abs(diag(IdentifiedMatches)-diag(Labels))) + 0.5*sum(sum(abs(triu(IdentifiedMatches,1)-triu(Labels,1)))); %Off diagonal counts more
    nplot = 3;
else
    nplot=2;
end


if DrawThis
    figure('name','TwoFactorPerformance')
    subplot(1,nplot,1)
    imagesc(IdentifiedMatches)
    xlabel('Unit_i'); ylabel('Unit_j'); title('TMS Matches'); colorbar; colormap(flipud(gray)); makepretty

    subplot(1,nplot,2)
    imagesc(RankBin)
    xlabel('Unit_i'); ylabel('Unit_j'); title(['Rank==1 Matches ' num2str(round(PerformanceRank*1000)./10) '%']); colorbar; colormap(flipud(gray)); makepretty

    if knownlabel
        subplot(1,nplot,3)
        imagesc(Labels)
        xlabel('Unit_i'); ylabel('Unit_j'); title(['PyKS Matches ' num2str(round(PerformanceLabels*1000)./10) '%']); colorbar; colormap(flipud(gray)); makepretty


    end

end
end